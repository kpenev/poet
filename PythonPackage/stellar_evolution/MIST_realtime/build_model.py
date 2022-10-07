#!/usr/bin/env python3
"""Define class for building and using MIST models in POET."""

from os import path, listdir, makedirs
from abc import ABC, abstractmethod
from tempfile import TemporaryDirectory
import shutil
import logging

try:
    import spython
except ModuleNotFoundError:
    import docker

from stellar_evolution.library_interface import library

_workdir_template = path.join(path.dirname(__file__), 'mesa_workdir_template')

def _init_mesa_workdir(workdir_name):
    """Add all non-template entries of _workdir_template to given directory."""

    for entry in listdir(_workdir_template):
        if (
                not entry.endswith('_template')
                and
                entry not in ['README.rst', '.gitignore']
        ):
            orig = path.join(_workdir_template, entry)
            copy = path.join(workdir_name, entry)
            try:
                shutil.copyfile(orig, copy)
                shutil.copymode(orig, copy)
            except IsADirectoryError:
                shutil.copytree(orig, copy)


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


    def __enter__(self):
        """Make temp directory and initialize it per `_init_mesa_workdir()."""

        workdir_name = super().__enter__()
        _init_mesa_workdir(workdir_name)
        return workdir_name


class MISTTrackMakerBase(ABC):
    """Class for generating and processing MIST stellar evolution tracks."""

    @staticmethod
    def _generate_inlist(mass,
                         feh,
                         max_age,
                         destination,
                         *,
                         profile_interval=None,
                         **config):
        """
        Generate the inlist to run MESA with per current stellar properties.

        Args:
            mass(float):    The mass of the star to simulate.

            feh(float):    The [Fe/H] of the star to simulate.

            max_age(float):    The last age to include in the generated
                evolution track.

            destination(str):    The filename of the inlist to generate.

            config:    Controls the precision with which stellar evolution is
                calculated. See MESA documentation for description of
                parameters. If unspecified, the following defaults are used:

                    * when_to_stop_rtol = 0

                    * when_to_stop_atol = 1e-6

                    * varcontrol_target = 3d-4

                    * delta_lgTeff_limit = 0.005

                    * delta_lgTeff_hard_limit = 0.01

                    * delta_lgL_limit = 0.02

                    * delta_lgL_hard_limit = 0.05


            profile_interval(int or None):    Whether to generate profiles and
                with what frequnecy (number of models between profiles).

        Returns:
            str:
                Name of the inlist file generated.
        """

        assert not path.exists(destination)

        initial_z = library.z_from_feh(feh)
        inlist_substitutions = dict(
            INITIAL_MASS=mass,
            INITIAL_Z=initial_z,
            MAX_AGE=max_age * 1e9,
            WHEN_TO_STOP_RTOL=0,
            WHEN_TO_STOP_ATOL=1.0,
            ENABLE_PROFILES=('true' if bool(profile_interval) else 'false'),
            PROFILE_INTERVAL=(profile_interval or 1),
            VARCONTROL_TARGET=3e-4,
            DELTA_LGTEFF_LIMIT=0.005,
            DELTA_LGTEFF_HARD_LIMIT=0.01,
            DELTA_LGL_LIMIT=0.02,
            DELTA_LGL_HARD_LIMIT=0.05
        )
        for k in config:
            inlist_substitutions[k.upper()] = config[k]

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

        inlist_substitutions['INITIAL_H1'] = initial_h / (1.0 + h2_h1_ratio)
        inlist_substitutions['INITIAL_H2'] = (
            inlist_substitutions['INITIAL_H1']
            *
            h2_h1_ratio
        )
        inlist_substitutions['INITIAL_HE4'] = initial_he / (1.0 + he3_he4_ratio)
        inlist_substitutions['INITIAL_HE3'] = (
            inlist_substitutions['INITIAL_HE4']
            *
            he3_he4_ratio
        )

        with open(destination, 'w') as inlist:
            inlist.write(
                open(
                    path.join(_workdir_template, 'inlist_template')
                ).read(
                ).format_map(
                    inlist_substitutions
                )
            )


    @abstractmethod
    def _execute_mesa_command(self, command, mesa_workdir):
        """Execute one of the MESA commands in the given working directory."""


    def _ensure_mesa_executable(self):
        """Check if MESA star executable exists and compile it if not."""

        expected_fname = path.join(_workdir_template, 'star')
        if path.exists(expected_fname):
            return

        with TemporaryMESAWorkDirectory() as mesa_workdir:
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


    def __init__(self):
        """Prepare to generate MIST tracks."""

        self._ensure_mesa_executable()


    def run_mesa(self, mass, feh, max_age, mesa_workdir, **mesa_config):
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

            mesa_config:    Additional configuration for how to run MESA (see
                `generate_inlist()` for allowed configurations

        Returns:
            None
        """

        self._generate_inlist(mass,
                              feh,
                              max_age,
                              path.join(mesa_workdir, 'inlist'),
                              **mesa_config)
        self._execute_mesa_command('./rn', mesa_workdir)


class DockerMISTTrackMaker(MISTTrackMakerBase):
    """Generate MIST tracks using mesa through docker."""

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
        exit_code, output = self._mesa.exec_run(command,
                                                workdir=container_workdir)
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
            docker.APIClient().inspect_container('mesa-r7503')['Mounts']
        )
        self._mesa_workdir_parent = None
        self._container_workdir_parent = container_workdir_parent
        for mount in container_mounts:
            if mount['Destination'] == container_workdir_parent:
                self._mesa_workdir_parent = mount['Source']

        if self._mesa_workdir_parent is None:
            raise RuntimeError('No mount point found for MESA working '
                               'directories in docker container '
                               +
                               container_name)

        TemporaryMESAWorkDirectory.set_workdir_parent(
            self._mesa_workdir_parent
        )

        if self._mesa.status != 'running':
            self._mesa.start()
        super().__init__()


class SingularityMISTTrackMaker(MISTTrackMakerBase):
    """Generate MIST tracks using mesa through singularity."""

    def _execute_mesa_command(self, command, mesa_workdir):
        raise NotImplementedError()


MISTTrackMaker = (SingularityMISTTrackMaker if 'spython' in globals()
                  else DockerMISTTrackMaker)
