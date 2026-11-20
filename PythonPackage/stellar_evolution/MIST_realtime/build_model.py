#!/usr/bin/env python3
"""Define class for building and using MIST models in POET."""

import os as _os
from os import path, listdir, makedirs, chmod, remove, rename
from abc import ABC, abstractmethod
from tempfile import TemporaryDirectory
import glob as _glob
import re as _re
import shutil
import logging

try:
    from spython.main import Client as SingularityClient
except ModuleNotFoundError:
    SingularityClient = None

try:
    import docker
except ModuleNotFoundError:
    docker = None
from mesa_reader import MesaData
import numpy
import pandas

from stellar_evolution.library_interface import library, MESAInterpolator

_workdir_template = path.join(path.dirname(__file__), 'mesa_workdir_template')


def mist2_workdir(base_dir, mass, feh, afe=0.0, vvcrit=0.0):
    """Return path to the pre-existing MIST2 workdir for the given parameters.

    Args:
        base_dir(str):  Root of the MIST2 data tree (contains
            ``nonrotating_workdirs/``).
        mass(float):    Stellar mass in M_sun.
        feh(float):     [Fe/H] metallicity.
        afe(float):     [alpha/Fe] enhancement (default 0.0).
        vvcrit(float):  v/v_crit rotation fraction (default 0.0).

    Returns:
        str: Absolute path to the matching MIST2 workdir.
    """
    feh_sign = 'm' if feh < 0 else 'p'
    feh_str = f'{abs(round(feh * 100)):03d}'
    afe_sign = 'm' if afe < 0 else 'p'
    afe_str = str(abs(round(afe * 10)))
    mass_str = f'{round(mass * 100):05d}M_dir'
    feh_dir = f'feh_{feh_sign}{feh_str}_afe_{afe_sign}{afe_str}_vvcrit{vvcrit:.1f}'
    return path.join(base_dir, 'nonrotating_workdirs', feh_dir, mass_str)


_FEH_DIR_PAT = _re.compile(
    r'^feh_([mp])(\d+)_afe_([mp])(\d+)_vvcrit([\d.]+)$'
)
_MASS_DIR_PAT = _re.compile(r'^(\d+)M(?:_VLM)?_dir$')


def _parse_mist2_grid(base_dir):
    """Return list of (feh, afe, vvcrit, mass, path) for all non-VLM MIST2 workdirs."""
    entries = []
    nwr = path.join(base_dir, 'nonrotating_workdirs')
    for feh_dir in listdir(nwr):
        m = _FEH_DIR_PAT.match(feh_dir)
        if not m:
            continue
        feh_val = (-1 if m.group(1) == 'm' else 1) * int(m.group(2)) / 100.0
        afe_val = (-1 if m.group(3) == 'm' else 1) * int(m.group(4)) / 10.0
        vvcrit_val = float(m.group(5))
        feh_path = path.join(nwr, feh_dir)
        for mass_dir in listdir(feh_path):
            mm = _MASS_DIR_PAT.match(mass_dir)
            if not mm:
                continue
            entries.append((feh_val, afe_val, vvcrit_val,
                            int(mm.group(1)) / 100.0,
                            path.join(feh_path, mass_dir)))
    return entries


def find_closest_mist2_workdir(base_dir, mass, feh, afe=0.0, vvcrit=0.0):
    """Find the MIST2 workdir with the closest composition then closest mass.

    Minimises (feh, afe) Euclidean distance first (composition and opacity
    are fixed per grid point), then picks the closest mass at that point.

    Args:
        base_dir(str):  Root of the MIST2 data tree.
        mass(float):    Target stellar mass in M_sun.
        feh(float):     Target [Fe/H].
        afe(float):     Target [alpha/Fe] (default 0.0).
        vvcrit(float):  v/v_crit rotation (default 0.0).

    Returns:
        str: Path to the best-matching MIST2 workdir.
    """
    entries = _parse_mist2_grid(base_dir)
    if not entries:
        raise RuntimeError(f'No MIST2 workdirs found under {base_dir!r}')
    best_fa = min(
        {(e[0], e[1]) for e in entries},
        key=lambda fa: (fa[0] - feh) ** 2 + (fa[1] - afe) ** 2,
    )
    candidates = [e for e in entries
                  if e[0] == best_fa[0] and e[1] == best_fa[1]]
    return min(candidates, key=lambda e: abs(e[3] - mass))[4]


def _read_grid_xyz(base_dir, afe=0.0, vvcrit=0.0):
    """Return {feh: (X, Y, Z)} from MIST2 grid for the nearest afe value."""
    entries = _parse_mist2_grid(base_dir)
    afe_vals = sorted({e[1] for e in entries if e[2] == vvcrit})
    if not afe_vals:
        raise RuntimeError(f'No MIST2 workdirs with vvcrit={vvcrit} under {base_dir!r}')
    closest_afe = min(afe_vals, key=lambda a: abs(a - afe))
    result = {}
    for e in entries:
        if e[1] != closest_afe or e[2] != vvcrit:
            continue
        gfeh = e[0]
        if gfeh in result:
            continue
        xyz_file = path.join(e[4], 'input_XYZ')
        if not path.exists(xyz_file):
            continue
        with open(xyz_file) as f:
            vals = [float(x) for x in f.read().split()]
        if len(vals) == 3:
            result[gfeh] = tuple(vals)
    return result


def _interpolate_xyz(grid_xyz, feh):
    """Linearly interpolate (X, Y, Z) to target feh from MIST2 grid data."""
    fehs = sorted(grid_xyz.keys())
    if len(fehs) < 2:
        return grid_xyz[fehs[0]]
    if feh <= fehs[0]:
        f1, f2 = fehs[0], fehs[1]
    elif feh >= fehs[-1]:
        f1, f2 = fehs[-2], fehs[-1]
    else:
        for i in range(len(fehs) - 1):
            if fehs[i] <= feh <= fehs[i + 1]:
                f1, f2 = fehs[i], fehs[i + 1]
                break
    t = (feh - f1) / (f2 - f1)
    x1, y1, z1 = grid_xyz[f1]
    x2, y2, z2 = grid_xyz[f2]
    return (x1 + t * (x2 - x1), y1 + t * (y2 - y1), z1 + t * (z2 - z1))


def _update_initial_abundances(mesa_workdir, feh, base_dir, afe=0.0, vvcrit=0.0):
    """
    Rewrite input_XYZ and input_initial_xa.data for the target (feh, afe).

    H species (h1, h2) are scaled by X_target/X_source preserving h2/h1 ratio.
    He species (he3, he4) are scaled by Y_target/Y_source preserving he3/he4 ratio.
    All metal species are scaled by Z_target/Z_source.

    Returns Z_target for use as Zbase in the inlist.
    """
    grid_xyz = _read_grid_xyz(base_dir, afe, vvcrit)
    X_target, Y_target, Z_target = _interpolate_xyz(grid_xyz, feh)

    xyz_path = path.join(mesa_workdir, 'input_XYZ')
    if not path.exists(xyz_path):
        return Z_target

    with open(xyz_path) as f:
        vals = [float(x) for x in f.read().split()]
    X_source, Y_source, Z_source = vals

    with open(xyz_path, 'w') as f:
        f.write(f'    {X_target:.10E}\n')
        f.write(f'    {Y_target:.10E}\n')
        f.write(f'    {Z_target:.10E}\n')

    xa_path = path.join(mesa_workdir, 'input_initial_xa.data')
    if not path.exists(xa_path):
        return Z_target

    x_scale = X_target / X_source if X_source > 0 else 1.0
    y_scale = Y_target / Y_source if Y_source > 0 else 1.0
    z_scale = Z_target / Z_source if Z_source > 0 else 1.0
    h_species  = {'h1', 'h2'}
    he_species = {'he3', 'he4'}

    lines_out = []
    with open(xa_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('!'):
                lines_out.append(line)
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                sp = parts[0].lower()
                val = float(parts[1])
                if sp in h_species:
                    new_val = val * x_scale
                elif sp in he_species:
                    new_val = val * y_scale
                else:
                    new_val = val * z_scale
                lines_out.append(f'{parts[0]:<12}{new_val:.10E}\n')
            else:
                lines_out.append(line)

    with open(xa_path, 'w') as f:
        f.writelines(lines_out)

    return Z_target


def _patch_inlist_project(inlist_path, mass, feh, max_age,
                           zbase=None,
                           when_to_stop_rtol=0, when_to_stop_atol=1000,
                           profile_interval=None,
                           varcontrol_target=3e-4,
                           delta_lgTeff_limit=0.005,
                           delta_lgTeff_hard_limit=0.01,
                           delta_lgL_limit=0.02,
                           delta_lgL_hard_limit=0.05):
    """Patch a MIST2 inlist_project in-place for a POET run.

    Changes made:
      - initial_mass set to target mass
      - star_history_name set to 'history.data'
      - extras_rpar(4) set to feh (needed by f_ov_below_nonburn in extras_startup)
      - extras_lpar(2) commented out (post-AGB burning suppression, not needed)
      - gamma_center_limit commented out (POET uses max_age / H1-depletion stop)
      - Zbase set to zbase (target Z for Type2 opacities)
      - write_profiles_flag updated
      - max_age, stopping criteria, and precision params appended before final /
    """
    with open(inlist_path) as f:
        lines = f.readlines()

    result = []
    for line in lines:
        bare = line.strip()
        if _re.match(r'\s*initial_mass\s*=', line):
            result.append(f'      initial_mass = {mass}d0\n')
        elif _re.match(r'\s*star_history_name\s*=', line):
            result.append("      star_history_name = 'history.data'\n")
        elif _re.match(r'\s*extras_rpar\s*\(\s*4\s*\)\s*=', line):
            result.append(f'      extras_rpar(4) = {feh}d0 ! [Fe/H]\n')
        elif _re.match(r'\s*extras_lpar\s*\(\s*2\s*\)\s*=', line):
            result.append('      !' + bare + '\n')
        elif _re.match(r'\s*gamma_center_limit\s*=', line):
            result.append('      !' + bare + '\n')
        elif zbase is not None and _re.match(r'\s*Zbase\s*=', line):
            result.append(f'      Zbase = {zbase!r}\n')
        elif _re.match(r'\s*write_profiles_flag\s*=', line):
            flag = 'true' if bool(profile_interval) else 'false'
            result.append(f'      write_profiles_flag = .{flag}.\n')
        else:
            result.append(line)

    # Insert POET stopping/precision params before the final closing '/'
    last_slash = next(
        (i for i in range(len(result) - 1, -1, -1)
         if result[i].strip() == '/'),
        None
    )
    if last_slash is not None:
        inserts = [
            f'      max_age = {max_age * 1e9!r}\n',
            f'      when_to_stop_rtol = {when_to_stop_rtol!r}\n',
            f'      when_to_stop_atol = {when_to_stop_atol!r}\n',
            f'      varcontrol_target = {varcontrol_target!r}\n',
            f'      delta_lgTeff_limit = {delta_lgTeff_limit!r}\n',
            f'      delta_lgTeff_hard_limit = {delta_lgTeff_hard_limit!r}\n',
            f'      delta_lgL_limit = {delta_lgL_limit!r}\n',
            f'      delta_lgL_hard_limit = {delta_lgL_hard_limit!r}\n',
        ]
        if profile_interval:
            inserts.append(f'      profile_interval = {int(profile_interval):d}\n')
        result[last_slash:last_slash] = inserts

    with open(inlist_path, 'w') as f:
        f.writelines(result)


class TemporaryMESAWorkDirectory(TemporaryDirectory):
    """Wrap around python temporary directories to add MESA contents."""

    @classmethod
    def set_workdir_parent(cls, workdir_parent):
        """Set directory under which MESA work directories should be created."""

        cls._workdir_parent = workdir_parent


    def __init__(self, *args, **kwargs):
        """Set the parent directory as needed for MESA."""

        if len(args) >= 3 or 'dir' in kwargs:
            raise RuntimeError('Argument `dir` must not be specified when '
                               'creating temporary MESA work directories.')
        super().__init__(*args, dir=self._workdir_parent, **kwargs)


def convert_history_to_track(initial_mass,
                             history_fname,
                             track_fname,
                             min_age=None,
                             zone_separated='normal'):
    """Convert a MESA history file to a track ready to interplate.

    Args:
        zone_separated(str):    Controls how the POET two-zone (conv/rad)
            inertias are assigned from the three MESA zones (I_cc, I_rs, I_cs):

            * ``'normal'``: dynamic per-timestep assignment based on which zone
              dominates — if I_cc > I_cs (core-dominant):
              I_conv = I_rs + I_cs, I_rad = I_cc;
              otherwise (shell-dominant): I_conv = I_cs, I_rad = I_cc + I_rs.

            * ``'rc_env'``: radiative core treated as envelope —
              I_conv = I_rs + I_cs, I_rad = I_cc (always).

            * ``'rc_core'``: radiative core treated as core —
              I_conv = I_cs, I_rad = I_cc + I_rs (always).
    """

    history = MesaData(history_fname)

    if min_age:
        selected_points = (history.star_age >= min_age * 1e9)
        num_points = selected_points.sum()
    else:
        selected_points = slice(None)
        num_points = history.star_age.size

    i_cc = history.core_inertia[selected_points]
    i_rs = history.rad_inertia[selected_points]
    i_cs = history.env_inertia[selected_points]
    r_core = history.rcore[selected_points]
    r_shell = history.rrad[selected_points]

    m_core_orig = history.mcore[selected_points]
    m_star = history.star_mass[selected_points]

    if zone_separated == 'rc_env':
        i_conv = i_rs + i_cs
        i_rad  = i_cc
        r_tachocline = r_core
        m_rad = numpy.maximum(m_core_orig, 1e-3)
        m_conv = m_star - m_rad
    elif zone_separated == 'rc_core':
        i_conv = i_cs
        i_rad  = i_cc + i_rs
        r_tachocline = r_shell
        m_rad_mesa = history.mrad[selected_points]
        m_rad = numpy.maximum(m_core_orig, m_rad_mesa)
        fully_radiative = (i_rad > 0) & (i_conv == 0)
        m_rad[fully_radiative] = m_star[fully_radiative]
        m_conv = m_star - m_rad
    else:
        # 'normal': dynamic per-timestep assignment — no permanent lock-in.
        # I_cc > I_cs: core-dominant  → I_conv = I_rs + I_cs, I_rad = I_cc
        # I_cc ≤ I_cs: shell-dominant → I_conv = I_cs,        I_rad = I_cc + I_rs
        n = len(i_cc)
        i_rad  = numpy.empty(n)
        i_conv = numpy.empty(n)
        r_tachocline = numpy.empty(n)

        core_mask  = i_cc > i_cs
        shell_mask = ~core_mask

        i_conv[core_mask]  = i_rs[core_mask]  + i_cs[core_mask]
        i_rad[core_mask]   = i_cc[core_mask]
        r_tachocline[core_mask] = r_core[core_mask]

        i_conv[shell_mask] = i_cs[shell_mask]
        i_rad[shell_mask]  = i_cc[shell_mask] + i_rs[shell_mask]
        r_tachocline[shell_mask] = r_shell[shell_mask]

        m_rad = m_core_orig
        m_conv = m_star - m_rad

    pandas.DataFrame(
        dict(age=history.star_age[selected_points],
             R_star=10.0**history.log_R[selected_points],
             L_star=10.0**history.log_L[selected_points],
             M_star=m_star,
             R_tachocline=r_tachocline,
             M_ini=numpy.full(num_points, initial_mass),
             M_conv=m_conv,
             M_rad=m_rad,
             I_conv=numpy.maximum(i_conv, 1e-4),
             I_rad=numpy.maximum(i_rad, 1e-4),
             center_h1=history.center_h1[selected_points])
    ).to_csv(
        track_fname,
        na_rep='NaN',
        index=False
    )


def _setup_mesa_workdir(mesa_workdir, substitutions):
    """Use given substitutions and templates to create MESA working dir."""


    ignore_workdir_template_files = ['README.rst', '.gitignore']

    for entry in listdir(_workdir_template):
        source = path.join(_workdir_template, entry)
        destination = path.join(mesa_workdir, entry)
        if entry.endswith('_template'):
            destination = destination[:-len('_template')]
            assert not path.exists(destination)
            with open(destination, 'w') as inlist:
                inlist.write(
                    open(source).read().format_map(substitutions)
                )
            shutil.copymode(source, destination)
        elif entry not in ignore_workdir_template_files:
            logging.getLogger(__name__).debug('Copying %s -> %s',
                                              source,
                                              destination)
            try:
                shutil.copyfile(source, destination)
                shutil.copymode(source, destination)
            except OSError:
                shutil.copytree(source, destination)


class MISTTrackMakerBase(ABC):
    """Class for generating and processing MIST stellar evolution tracks."""

    @classmethod
    def _prepare_run(cls,
                     mass,
                     feh,
                     max_age,
                     mesa_workdir,
                     *,
                     profile_interval=None,
                     num_threads=4,
                     afe=0.0,
                     mist2_base_dir=None,
                     **config):
        """
        Generate required inputs to run MESA for specified stellar properties.

        Args:
            mass(float):    The mass of the star to simulate.

            feh(float):    The [Fe/H] of the star to simulate.

            max_age(float):    The last age to include in the generated
                evolution track.

            mesa_workdir(str):    The MESA working directory to configure.

            afe(float):     [alpha/Fe] enhancement (default 0.0, MIST2 only).

            mist2_base_dir(str or None):    Root of the MIST2 data tree.  When
                set, the closest MIST2 workdir is used as the base: its files
                are copied, POET's run_star_extras is substituted, the inlist
                is patched, and a .bashrc-sourcing mk/rn pair is generated.
                When None, the legacy POET template (MIST 1.0 physics) is used.

            config:    See ``mesa_config`` argument to `run_mesa()`

            profile_interval(int or None):    Whether to generate profiles and
                with what frequency (number of models between profiles).

            num_threads(int):    The maximum number of threads to allow MESA to
                use.
        """

        if mist2_base_dir is not None:
            # === MIST2-based approach ===
            source_workdir = find_closest_mist2_workdir(
                mist2_base_dir, mass, feh, afe
            )
            logging.getLogger(__name__).debug(
                'Using MIST2 workdir %s for M=%s feh=%s afe=%s',
                source_workdir, mass, feh, afe
            )

            # Copy the entire MIST2 workdir into mesa_workdir.
            for entry in listdir(source_workdir):
                src = path.join(source_workdir, entry)
                dst = path.join(mesa_workdir, entry)
                try:
                    shutil.copytree(src, dst)
                except NotADirectoryError:
                    shutil.copy2(src, dst)

            # Replace MIST2's run_star_extras.f90 with POET's version.
            # Keep the .f90 extension so work_standard_makefile finds it.
            # POET's file is free-format F90 despite its .f extension.
            # Also remove use mlt_def — in r7503 those constants live in const_def.
            f90_dst = path.join(mesa_workdir, 'src', 'run_star_extras.f90')
            if path.exists(f90_dst):
                remove(f90_dst)
            shutil.copy2(
                path.join(_workdir_template, 'src', 'run_star_extras.f'),
                f90_dst
            )

            # Update input_XYZ and input_initial_xa.data for target (feh, afe).
            # Returns Z_target for Zbase in the inlist.
            z_target = _update_initial_abundances(
                mesa_workdir, feh, mist2_base_dir, afe=afe
            )

            # Patch inlist_project for POET (mass, stopping criteria, etc.).
            patch_kwargs = dict(
                varcontrol_target=3e-4,
                delta_lgTeff_limit=0.005,
                delta_lgTeff_hard_limit=0.01,
                delta_lgL_limit=0.02,
                delta_lgL_hard_limit=0.05,
            )
            patch_kwargs.update({k.lower(): v for k, v in config.items()})
            _patch_inlist_project(
                path.join(mesa_workdir, 'inlist_project'),
                mass=mass, feh=feh, max_age=max_age,
                zbase=z_target,
                when_to_stop_rtol=0, when_to_stop_atol=1000,
                profile_interval=profile_interval,
                **patch_kwargs
            )
            # Rename to 'inlist' so ./star (invoked by rn) reads it by default.
            rename(
                path.join(mesa_workdir, 'inlist_project'),
                path.join(mesa_workdir, 'inlist')
            )

            # Generate rn from POET's rn_template (sources .bashrc for env).
            rn_tmpl = path.join(_workdir_template, 'rn_template')
            with open(rn_tmpl) as _f:
                rn_content = _f.read().format_map({'NUM_THREADS': num_threads})
            rn_path = path.join(mesa_workdir, 'rn')
            with open(rn_path, 'w') as _f:
                _f.write(rn_content)
            shutil.copymode(rn_tmpl, rn_path)

            # Generate an mk wrapper that sources .bashrc before building
            # (the MIST2 mk script runs bare make without MESASDK in env).
            mk_path = path.join(mesa_workdir, 'mk')
            with open(mk_path, 'w') as _f:
                _f.write(
                    '#!/bin/bash -il\n\n'
                    'export HOME=/home/user\n\n'
                    'source /home/user/.bashrc\n\n'
                    'function check_okay {\n'
                    '    if [ $? -ne 0 ]; then\n'
                    '        echo\n'
                    '        echo "FAILED"\n'
                    '        echo\n'
                    '        exit 1\n'
                    '    fi\n'
                    '}\n\n'
                    'cd make; make; check_okay\n'
                )
            shutil.copymode(rn_tmpl, mk_path)

            # Make the entire workdir writable by the container user.
            chmod(mesa_workdir, 0o777)
            for _root, _dirs, _files in _os.walk(mesa_workdir):
                for _d in _dirs:
                    chmod(path.join(_root, _d), 0o777)
                for _f in _files:
                    chmod(path.join(_root, _f), 0o777)

        else:
            # === Legacy POET template approach (MIST 1.0) ===
            initial_z = library.z_from_feh(feh)
            substitutions = dict(
                INITIAL_MASS=mass,
                INITIAL_Z=initial_z,
                MAX_AGE=max_age * 1e9,
                WHEN_TO_STOP_RTOL=0,
                WHEN_TO_STOP_ATOL=1000,
                ENABLE_PROFILES=('true' if bool(profile_interval) else 'false'),
                PROFILE_INTERVAL=(profile_interval or 1),
                VARCONTROL_TARGET=3e-4,
                DELTA_LGTEFF_LIMIT=0.005,
                DELTA_LGTEFF_HARD_LIMIT=0.01,
                DELTA_LGL_LIMIT=0.02,
                DELTA_LGL_HARD_LIMIT=0.05
            )
            for k in config:
                substitutions[k.upper()] = config[k]
            h2_h1_ratio = 2e-5
            he3_he4_ratio = 1.66e-4
            initial_he = (
                library.constants.Yprimordial
                +
                (library.constants.Yprotosun - library.constants.Yprimordial)
                /
                library.constants.Zprotosun
                *
                initial_z
            )
            initial_h = 1.0 - initial_he - initial_z
            h1 = initial_h / (1.0 + h2_h1_ratio)
            h2 = h1 * h2_h1_ratio
            he4 = initial_he / (1.0 + he3_he4_ratio)
            he3 = he4 * he3_he4_ratio
            substitutions['COMPOSITION_BLOCK'] = (
                "      set_uniform_initial_composition = .true.\n"
                f"      initial_h1 = {h1!r}\n"
                f"      initial_h2 = {h2!r}\n"
                f"      initial_he3 = {he3!r}\n"
                f"      initial_he4 = {he4!r}\n"
                f"      !for initial Z = {initial_z!r}\n\n"
                "      initial_zfracs = 6 ! AGSS09_zfracs"
            )
            substitutions['NUM_THREADS'] = num_threads
            _setup_mesa_workdir(mesa_workdir, substitutions)


    @abstractmethod
    def _execute_container_command(self, command, container_workdir):
        """Run command in the MESA container from the given container dir."""


    def _execute_mesa_command(self, command, mesa_workdir):
        """Execute one of the MESA commands in the given working directory."""

        container_workdir = path.join(
            self._container_workdir_parent,
            path.relpath(mesa_workdir, self._mesa_workdir_parent)
        )
        logging.getLogger(__name__).debug(
            'Starting MESA command %s in %s -> %s working directory in docker',
            repr(command),
            repr(mesa_workdir),
            repr(container_workdir)
        )
        exit_code, output = self._execute_container_command(
            command,
            container_workdir
        )
        if exit_code:
            raise RuntimeError(
                'MESA command {!r} in {!r} -> {!r} working directory failed '
                'with return code {:d} and output:\n{!s}'.format(
                    command,
                    mesa_workdir,
                    container_workdir,
                    exit_code,
                    output
                )
            )



    def _ensure_mesa_executable(self):
        """Check if MESA star executable exists and compile it if not."""

        expected_fname = path.join(_workdir_template, 'star')
        if path.exists(expected_fname):
            return

        with TemporaryMESAWorkDirectory() as mesa_workdir:
            chmod(mesa_workdir, 0o777)
            for entry in listdir(_workdir_template):
                src = path.join(_workdir_template, entry)
                dst = path.join(mesa_workdir, entry)
                if path.isdir(src):
                    shutil.copytree(src, dst)
                else:
                    shutil.copy2(src, dst)
            # Recursively make everything world-writable so the Docker
            # container (uid=1000) can write .mod0 and object files.
            for root, dirs, files in _os.walk(mesa_workdir):
                for d in dirs:
                    chmod(_os.path.join(root, d), 0o777)
                for f in files:
                    chmod(_os.path.join(root, f), 0o777)
            self._execute_mesa_command('./mk', mesa_workdir)
            logging.getLogger(__name__).debug(
                'After compiling, %s contents are: %s',
                repr(mesa_workdir),
                repr(listdir(mesa_workdir))
            )
            created_fname = path.join(mesa_workdir, 'star')
            assert path.exists(created_fname)
            shutil.copyfile(created_fname, expected_fname)
            shutil.copystat(created_fname, expected_fname)


    def __init__(self, mesa_workdir_parent, container_workdir_parent):
        """Prepare to generate MIST tracks."""

        self._mesa_workdir_parent = mesa_workdir_parent
        self._container_workdir_parent = container_workdir_parent
        TemporaryMESAWorkDirectory.set_workdir_parent(mesa_workdir_parent)
        self._ensure_mesa_executable()


    def run_mesa(self, mass, feh, max_age, mesa_workdir, *, afe=0.0, mist2_base_dir=None, **mesa_config):
        """
        Use MESA to generate a MIST stellar evolution track in given work dir.

        Args:
            mass(float):    The mass of the star to simulate.

            feh(float):    The [Fe/H] of the star to simulate.

            max_age(float):    The last age to include in the generated
                evolution track.

            mesa_workdir(str):    An initialized (per `_init_mesa_workdir()`)
                MESA working directory to generate history and possibly profiles
                in.

            mesa_config:    Controls the following::

                * the precision with which stellar evolution is calculated. See
                  MESA documentation for description of parameters. If
                  unspecified, the following precision defaults are used and no
                  profiles are generated:

                    * when_to_stop_rtol = 0

                    * when_to_stop_atol = 1e-6

                    * varcontrol_target = 3d-4

                    * delta_lgTeff_limit = 0.005

                    * delta_lgTeff_hard_limit = 0.01

                    * delta_lgL_limit = 0.02

                    * delta_lgL_hard_limit = 0.05


                * profile_interval: Profile is output when step number is
                  divisible by this value.

                * num_threads: maximum number of threads to use (see MESA use of
                  OMP_NUM_THREADS environment variable).

        Returns:
            None
        """

        self._prepare_run(mass, feh, max_age, mesa_workdir,
                          afe=afe, mist2_base_dir=mist2_base_dir,
                          **mesa_config)
        if mist2_base_dir is not None:
            self._execute_mesa_command('./mk', mesa_workdir)
        self._execute_mesa_command('./rn', mesa_workdir)


    def create_interpolator(self,
                            mass,
                            feh,
                            max_age,
                            mesa_workdir,
                            *,
                            min_age=None,
                            zone_separated='normal',
                            interpolation_config=None,
                            force_rerun=False,
                            afe=0.0,
                            mist2_base_dir=None,
                            **mesa_config):
        """
        Create interpolator from existing history or creating one if needed.

        Args:
            mass(float):    The mass of the star to create an interpolator for.

            feh(float):    The [Fe/H] of the star to create an interpolator for.

            max_age(float):    The latest age at which the interpolator can be
                evaluated.

            mesa_workdir(str):    Path to directory either containing a MESA run
                or which will be used for a fresh MESA run. If a history file is
                found, it is assumed to apply for the given star.

            min_age(float):    If specified, entries from the MESA history file
                for ages less than this in Gyr are not used in the interplation.

            interpolation_config(dict):    See
            `MESAInterpolator.get_create_interpolator_config()` for full
            description.

            mesa_config:    See `run_mesa()`.

        Returns:
            MESAInterpolator:
                The interpolator allowing evaluating POET relevant quantities
                for the given star.
        """

        history_fname = path.join(mesa_workdir, 'LOGS', 'history.data')
        track_dir = path.join(mesa_workdir, 'TRACK')
        track_fname = path.join(
            track_dir,
            'M{!r}_Z{!r}.csv'.format(mass, library.z_from_feh(feh))
        )
        if force_rerun or not path.exists(track_fname):
            if force_rerun or not path.exists(history_fname):
                self.run_mesa(mass, feh, max_age, mesa_workdir,
                              afe=afe, mist2_base_dir=mist2_base_dir,
                              **mesa_config)
                assert path.exists(history_fname)
            if not path.exists(track_dir):
                makedirs(track_dir)
            convert_history_to_track(mass, history_fname, track_fname, min_age,
                                     zone_separated)

        if interpolation_config is None:
            interpolation_config = dict()
        return MESAInterpolator(
            mesa_dir=track_dir,
            **MESAInterpolator.get_create_interpolator_config(
                **interpolation_config
            ),
            num_threads=interpolation_config.get('num_threads', 1)
        )


class DockerMISTTrackMaker(MISTTrackMakerBase):
    """Generate MIST tracks using mesa through docker."""

    def _execute_container_command(self, command, container_workdir):
        """Execute one of the MESA commands in the given working directory."""

        sdk_lib = '/home/user/mesasdk/lib:/home/user/mesasdk/lib64'
        return self._mesa.exec_run(
            command,
            workdir=container_workdir,
            environment={'LD_LIBRARY_PATH': sdk_lib}
        )


    def __init__(self,
                 container_name='mesa-r7503',
                 container_workdir_parent='/home/user/mnt'):
        """
        Start the container with MESA v7503.

        Args:
            container_name(str):    The name of the docker container where MESA
                was compiled.
        """

        self._mesa = docker.from_env().containers.get(container_name)
        container_mounts = (
            docker.APIClient().inspect_container(container_name)['Mounts']
        )
        mesa_workdir_parent = None

        for mount in container_mounts:
            if mount['Destination'] == container_workdir_parent:
                mesa_workdir_parent = mount['Source']

        if mesa_workdir_parent is None:
            raise RuntimeError('No mount point found for MESA working '
                               'directories in docker container '
                               +
                               container_name)

        if self._mesa.status != 'running':
            self._mesa.start()
        super().__init__(mesa_workdir_parent, container_workdir_parent)


class SingularityMISTTrackMaker(MISTTrackMakerBase):
    """Generate MIST tracks using mesa through singularity."""

    def _execute_container_command(self, command, container_workdir):
        """Execute one of the MESA commands in the given working directory."""

        result = SingularityClient.execute(self._mesa,
                                           command,
                                           options=['--pwd', container_workdir],
                                           return_result=True)
        return result['return_code'], result['message']


    def __init__(
            self,
            mesa_dir='/work/05392/kpenev/shared/mesa/mesa-r7503',
            container_fname='/work/05392/kpenev/shared/mesa/nudome_14.0.sif',
            mesa_workdir_parent=(
                '/work/05392/kpenev/shared/mesa/MIST_workdir_root'
            ),
            container_workdir_parent='/home/user/mnt'
    ):
        """Start the container instance to run MESA v7503 in."""

        self._mesa = SingularityClient.instance(
            container_fname,
            options=[
                '--cleanenv',
                '--hostname', 'mesa-r7503',
                '--no-home',
                '--bind', ','.join([
                    mesa_workdir_parent + ':' + container_workdir_parent,
                    mesa_dir + ':/home/user/mesa'
                ])
            ]
        )
        super().__init__(mesa_workdir_parent, container_workdir_parent)

MISTTrackMaker = (SingularityMISTTrackMaker if SingularityClient is not None
                  else DockerMISTTrackMaker)
