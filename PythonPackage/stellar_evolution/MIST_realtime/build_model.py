#!/usr/bin/env python3
"""Define class for building and using MIST models in POET."""

from os import path, listdir, makedirs
from abc import ABC, abstractmethod
from tempfile import TemporaryDirectory
import shutil
import logging

try:
    from spython.main import Client as SingularityClient
except ModuleNotFoundError:
    import docker
from mesa_reader import MesaData
import pandas

from stellar_evolution.library_interface import library, MESAInterpolator

_workdir_template = path.join(path.dirname(__file__), 'mesa_workdir_template')


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
                             min_age=None):
    """Convert a MESA history file to a track ready to interplate."""

    history = MesaData(history_fname)

    if min_age:
        selected_points = (history.star_age >= min_age)
    else:
        selected_points = slice(None)

    pandas.DataFrame(
        dict(age=history.star_age[selected_points],
             R_star=10.0**history.log_R[selected_points],
             L_star=10.0**history.log_L[selected_points],
             M_star=history.star_mass[selected_points],
             R_tachocline=history.rcore[selected_points],
             M_ini=initial_mass[selected_points],
             M_conv=(history.star_mass - history.mcore)[selected_points],
             M_rad=history.mcore[selected_points],
             I_conv=history.env_inertia[selected_points],
             I_rad=history.core_inertia[selected_points])
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
        elif entry not in ignore_workdir_template_files:
            try:
                shutil.copyfile(source, destination)
                shutil.copymode(source, destination)
            except OSError:
                shutil.copytree(source, destination)


class MISTTrackMakerBase(ABC):
    """Class for generating and processing MIST stellar evolution tracks."""

    #TODO: Implement more general templating of any file in mesa work dir and
    #use to allow specifying number of threads to use

    @classmethod
    def _prepare_run(cls,
                     mass,
                     feh,
                     max_age,
                     mesa_workdir,
                     *,
                     profile_interval=None,
                     num_threads=1,
                     **config):
        """
        Generate required inputs to run MESA for specified stellar properties.

        Args:
            mass(float):    The mass of the star to simulate.

            feh(float):    The [Fe/H] of the star to simulate.

            max_age(float):    The last age to include in the generated
                evolution track.

            mesa_workdir(str):    The MESA working directory to configure.

            config:    See ``mesa_config`` argument to `run_mesa()`

            profile_interval(int or None):    Whether to generate profiles and
                with what frequnecy (number of models between profiles).

            num_threads(int):    The maximum number of threads to allow MESA to
                use.

        Returns:
            str:
                Name of the inlist file generated.
        """

        initial_z = library.z_from_feh(feh)
        substitutions = dict(
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

        substitutions['INITIAL_H1'] = initial_h / (1.0 + h2_h1_ratio)
        substitutions['INITIAL_H2'] = (
            substitutions['INITIAL_H1']
            *
            h2_h1_ratio
        )
        substitutions['INITIAL_HE4'] = initial_he / (1.0 + he3_he4_ratio)
        substitutions['INITIAL_HE3'] = (
            substitutions['INITIAL_HE4']
            *
            he3_he4_ratio
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

        self._prepare_run(mass,
                          feh,
                          max_age,
                          path.join(mesa_workdir, 'inlist'),
                          **mesa_config)
        self._execute_mesa_command('./rn', mesa_workdir)


    def create_interpolator(self,
                            mass,
                            feh,
                            max_age,
                            mesa_workdir,
                            *,
                            min_age=None,
                            interpolation_config=None,
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
        if not path.exists(track_fname):
            if not path.exists(history_fname):
                self.run_mesa(mass, feh, max_age, mesa_workdir, **mesa_config)
                assert path.exists(history_fname)
            if not path.exists(track_dir):
                makedirs(track_dir)
            convert_history_to_track(mass, history_fname, track_fname, min_age)

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

        return self._mesa.exec_run(command, workdir=container_workdir)


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


MISTTrackMaker = (SingularityMISTTrackMaker if 'SingularityClient' in globals()
                  else DockerMISTTrackMaker)
