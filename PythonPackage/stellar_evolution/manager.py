"""Define class for managing many stellar evolution interpolations."""

import os
import os.path
import shutil
from uuid import uuid4 as get_uuid
from decimal import Decimal
import ctypes
from math import isnan
import sys

import numpy
from sqlalchemy import exists, and_, func, create_engine

sys.path.append('..')

#Need to add POET package to module search path before importing
#pylint: disable=wrong-import-position
from stellar_evolution.change_variables import VarChangingInterpolator
from stellar_evolution.library_interface import\
    library_track_fname_rex,\
    library_track_fname,\
    library
from stellar_evolution import Session
from stellar_evolution.basic_utils import db_session_scope, tempdir_scope
from basic_utils import Structure

from .manager_data_model import\
    DataModelBase,\
    SerializedInterpolator,\
    InterpolationParameters,\
    Quantity,\
    Track,\
    ModelSuite

from .managed_interpolator import\
    ManagedInterpolator,\
    verify_checksum,\
    checksum_filename
#pylint: enable=wrong-import-position

class StellarEvolutionManager:
    """
    Class for managing a collection of stellar evolution inteprolations.
    """

    @staticmethod
    def _get_decimal(value):
        """Prepare a value for db storage as a limited precision decimal."""

        return Decimal(value).quantize(Decimal('0.0001'))

    @staticmethod
    def _define_evolution_quantities(db_session):
        """
        Define the quantities tracked by VarChangingInterpolator instances.

        Args:
            db_session:    The currently active database session.

        Returns:
            None
        """

        db_quantities = [
            Quantity(id=id, name=name)
            for name, id in VarChangingInterpolator.quantity_ids.items()
        ]
        db_session.add_all(db_quantities)

    def _initialize_database(self, db_engine, db_session):
        """
        Ensure all database tables exist and contain at least required data.

        Args:
            db_session:    An sqlalchemy session, used to update the database.

        Returns:
            None
        """

        for table in DataModelBase.metadata.sorted_tables:
            if not db_engine.has_table(table.name):
                table.create(db_engine)
                if table.name == 'quantities':
                    self._define_evolution_quantities(db_session)
                elif table.name == 'model_suites':
                    db_session.add(ModelSuite(name='MESA'))

    def _get_db_config(self, db_session):
        """Read some configuration from the database."""

        self._quantities = [
            Structure(id=q.id, name=q.name)
            for q in db_session.query(Quantity).all()
        ]

        self._suites = db_session.query(ModelSuite).all()
        self._new_track_id = db_session.query(Track).count() + 1

    #Not simplifiable without loss of intuitive use.
    #pylint: disable=too-many-arguments
    def _add_track(self,
                   track_fname,
                   mass,
                   feh,
                   model_suite,
                   db_session):
        """Add a track to the database."""

        db_track = Track(id=self._new_track_id,
                         filename=os.path.abspath(track_fname),
                         mass=mass,
                         feh=feh,
                         checksum=checksum_filename(track_fname),
                         suite=db_session.query(
                             ModelSuite
                         ).filter_by(
                             name=model_suite
                         ).one())
        db_session.add(db_track)
        self._new_track_id += 1
    #pylint: enable=too-many-arguments

    @staticmethod
    def _track_grid_from_files(track_fnames,
                               db_session,
                               model_suite=None):
        """
        Organize track files in a mass - [Fe/H] grid and verify checksums.

        In addition, if model_suite is not None, all tracks are verified to
        belong to this suite.

        Args:
            track_fnames:    See get_interpolator track_fnames argument.

        Returns:
            dict:
                keys - the masses of tracks; and values - further dictionaries
                with keys the [Fe/H] of tracks and values the filename of each
                track.
        """

        track_grid = dict()
        for fname in track_fnames:
            absolute_fname = os.path.abspath(fname)
            db_track = db_session.query(Track).filter_by(
                filename=absolute_fname
            )

            if model_suite is not None:
                assert db_track.suite.name == model_suite
            verify_checksum(fname, db_track.checksum, 'track')

            mass_key = db_track.mass
            feh_key = db_track.feh
            if mass_key not in track_grid:
                track_grid[mass_key] = dict()
            assert feh_key not in track_grid[mass_key]
            track_grid[mass_key][feh_key] = (fname, db_track.id)
        return track_grid

    def _track_grid_from_grid(self,
                              mass_list,
                              feh_list,
                              model_suite,
                              db_session):
        """
        Return a mass - [Fe/H] grid with filenames and checksums.

        Fails if multiple tracks are registered for some (mass, [Fe/H],
        model suite) combination.

        Args:
            mass_list:    The masses for which to include tracks.

            feh_list:    The [Fe/H values for which to include tracks.

            model_suite:    The software suite from whose tracks to choose.

            db_session:    A database session to submit queries to.

        Returns:
            dict:
                See _track_grid_from_files()
        """

        if mass_list is None:
            mass_list = [
                record[0]
                for record in db_session.query(Track.mass).filter(
                    Track.model_suite_id == ModelSuite.id,
                    ModelSuite.name == model_suite
                ).distinct().order_by(Track.mass).all()
            ]
        if feh_list is None:
            feh_list = [
                record[0]
                for record in db_session.query(Track.feh).filter(
                    Track.model_suite_id == ModelSuite.id,
                    ModelSuite.name == model_suite
                ).distinct().order_by(Track.feh).all()
            ]
        track_grid = {m: dict() for m in mass_list}
        for mass in mass_list:
            for feh in feh_list:
                db_track = db_session.query(Track).filter(
                    Track.mass == self._get_decimal(mass),
                    Track.feh == self._get_decimal(feh),
                    Track.model_suite_id == ModelSuite.id,
                    ModelSuite.name == model_suite
                ).one()

                verify_checksum(db_track.filename,
                                db_track.checksum,
                                'track')

                track_grid[mass][feh] = (db_track.filename, db_track.id)
        return track_grid

    def _find_existing_interpolator(self,
                                    *,
                                    track_grid,
                                    nodes,
                                    smoothing,
                                    db_session,
                                    vs_log_age,
                                    log_quantity):
        """
        Return the specified interpolation if already exists, otherwise None.

        Args:
            track_grid:    See result of _track_grid_from_files()

            nodes:    see get_interpolator.

            smoothing:    see get_interpolator.

            db_session:    The currently active database session.

            vs_log_age:    see get_interpolator.

            log_quantity:    see get_interpolator.

        Returns:
            VarChangingInterpolator:
                Pre-serialized interpolation matching the given arguments if
                one is found in the interpolation archive. If no
                pre-serialized interpolation exists, returns None.
        """


        num_tracks = len(track_grid) * len(next(iter(track_grid.values())))

        track_counts = db_session.query(
            SerializedInterpolator.id,
            func.count('*').label('num_tracks')
        ).join(
            SerializedInterpolator.tracks
        ).group_by(
            SerializedInterpolator.id
        ).subquery()

        match_config = (
            [
                exists().where(
                    and_(
                        InterpolationParameters.interpolator_id
                        ==
                        SerializedInterpolator.id
                        ,
                        InterpolationParameters.quantity_id == quantity.id
                        ,
                        InterpolationParameters.smoothing
                        ==
                        (None if isnan(smoothing[quantity.name])
                         else self._get_decimal(smoothing[quantity.name]))
                        ,
                        InterpolationParameters.nodes
                        ==
                        nodes[quantity.name],
                        InterpolationParameters.vs_log_age
                        ==
                        vs_log_age[quantity.name],
                        InterpolationParameters.log_quantity
                        ==
                        log_quantity[quantity.name]
                    )
                )
                for quantity in self._quantities
            ]
            +
            [
                SerializedInterpolator.tracks.any(Track.id == track_id)
                for mass, mass_row in track_grid.items()
                for feh, (track_fname, track_id) in mass_row.items()
            ]
            +
            [track_counts.c.num_tracks == num_tracks]
        )

        result = db_session.query(
            SerializedInterpolator
        ).join(
            track_counts,
            SerializedInterpolator.id == track_counts.c.id
        ).filter(
            *match_config
        ).one_or_none()
        if result is None:
            return result

        return ManagedInterpolator(
            db_interpolator=result,
            serialization_path=self._serialization_path,
            db_session=db_session
        )

    #TODO: think about breaking this function up.
    #pylint: disable=too-many-locals
    def _create_new_interpolator(self,
                                 *,
                                 track_grid,
                                 nodes,
                                 smoothing,
                                 db_session,
                                 vs_log_age,
                                 log_quantity,
                                 num_threads,
                                 name=None):
        """
        Generate the specified interpolation and add it to the archive.

        Args:
            track_grid:    See result of _track_grid_from_files()

            nodes:    see get_interpolator().

            smoothing:    see get_interpolator().

            vs_log_age:    see get_interpolator().

            log_quantity:    see get_interpolator().

            num_threads:    The number of simultaneous threads to use when
                constructing the interpolation.

            db_session:    The database query session to use.

            name:    The name to assign to the new interpolator. If None, the
                UUID used to form the filename is used.

        Returns:
            VarChangingInterpolator:
                Created from scratch based on the given arguments.
        """

        interp_str = str(get_uuid())
        interp_fname = os.path.join(self._serialization_path,
                                    interp_str)

        db_interpolator = SerializedInterpolator(
            id=(db_session.query(SerializedInterpolator).count() + 1),
            name=(name or interp_str),
            filename=interp_str
        )

        if (
                db_session.query(
                    SerializedInterpolator
                ).filter_by(name=name).count()
        ):
            raise ValueError('Interpolator named %s already exists, with a '
                             'different configuration than the one being '
                             'constructed!' % repr(name))

        with tempdir_scope() as track_dir:
            for mass, mass_row in track_grid.items():
                for feh, (track_fname, track_id) in mass_row.items():
                    track = db_session.query(
                        Track
                    ).filter_by(
                        id=track_id
                    ).one()
                    db_interpolator.tracks.append(track)
                    shutil.copy(
                        track_fname,
                        os.path.join(track_dir,
                                     library_track_fname(mass, feh))
                    )
            interp_smoothing = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype=ctypes.c_double
            )
            interp_nodes = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype=ctypes.c_int
            )
            interp_vs_log_age = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype=ctypes.c_bool
            )
            interp_log_quantity = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype=ctypes.c_bool
            )
            for q_name, q_index in \
                    VarChangingInterpolator.quantity_ids.items():
                interp_smoothing[q_index] = smoothing[q_name]
                interp_nodes[q_index] = nodes[q_name]
                interp_vs_log_age[q_index] = vs_log_age[q_name]
                interp_log_quantity[q_index] = log_quantity[q_name]

            db_interpolator.parameters = [
                InterpolationParameters(quantity_id=q.id,
                                        nodes=nodes[q.name],
                                        smoothing=smoothing[q.name],
                                        vs_log_age=vs_log_age[q.name],
                                        log_quantity=log_quantity[q.name],
                                        interpolator=db_interpolator)
                for q in self._quantities
            ]

            actual_interpolator = ManagedInterpolator(
                db_interpolator=db_interpolator,
                serialization_path=self._serialization_path,
                db_session=db_session,
                mesa_dir=track_dir,
                smoothing=interp_smoothing,
                nodes=interp_nodes,
                vs_log_age=interp_vs_log_age,
                log_quantity=interp_log_quantity,
                num_threads=num_threads
            )

        actual_interpolator.save(interp_fname)
        db_interpolator.checksum = checksum_filename(interp_fname)
        db_session.add(db_interpolator)
        db_session.add_all(db_interpolator.parameters)

        return actual_interpolator
    #pylint: enable=too-many-locals

    def __init__(self, serialization_path):
        """
        Create a manager storing serialized interpolators in the given path.

        Args:
            serialization_path:    The path where to store
                serialized interpolators.

        Returns:
            None
        """

        if not os.path.exists(serialization_path):
            os.makedirs(serialization_path)
        db_engine = create_engine(
            'sqlite:///'
            +
            os.path.join(serialization_path, 'serialized.sqlite'),
            echo=False
        )
        Session.configure(bind=db_engine)
        self._serialization_path = serialization_path
        with db_session_scope() as db_session:
            self._initialize_database(db_engine, db_session)
        with db_session_scope() as db_session:
            self._get_db_config(db_session)

    def get_interpolator(
            self,
            *,
            nodes=VarChangingInterpolator.default_nodes,
            smoothing=VarChangingInterpolator.default_smoothing,
            vs_log_age=VarChangingInterpolator.default_vs_log_age,
            log_quantity=VarChangingInterpolator.default_log_quantity,
            track_fnames=None,
            masses=None,
            feh=None,
            model_suite='MESA',
            new_interp_name=None,
            num_threads=1
    ):
        """
        Return a stellar evolution interpolator with the given configuration.

        All tracks that the interpolator should be based on must be
        pre-registered with the manager. Two ways are supported for
        identifying tracks: as a list of filenames or as a mass-[Fe/H]
        grid combined with a suite. The first case always works, while the
        second requires that the set of identified tracks is unique, i.e. for
        none of the mass - [Fe/H] combinations there are two or more
        tracks registered for the given suite.

        Args:
            nodes:    The number of nodes to use for the age interpolation of
                each quantity of each track. Should be a dictionary with keys
                VarChangingInterpolator.quantity_list. See the POET code
                StellarEvolution::Interpolator::create_from() documentation for
                a description of what this actually means.

            smoothing:    The amount of smoothing to use for the age
                interpolation of each quantity of each track. Should be a
                dictionary with keys VarChangingInterpolator.quantity_list. See
                the POET code StellarEvolution::Interpolator::create_from()
                documentation for a description of what this actually means.

            vs_log_age:    Use log(age) instead of age as the independent
                argument for the intperpolation? Should be a dictionary with
                keys VarChangingInterpolator.quantity_list.

            log_quantity:    Interpolate log(quantity) instead of quantity?
                Should be a dictionary with keys
                VarChangingInterpolator.quantity_list.

            track_fnames:    A list of files containing stellar evolution tracks
                the interpolator should be based on.

            masses:    A list of the stellar masses to include in the
                interpolation. Unique tracks with those masses and all selected
                [Fe/H] (see next argument) must already be registered
                with the database for the given suite. If None, all track
                masses from the given suite are used.

            feh:    A list of the stellar [Fe/H] values to include in the
                interpolation. If None, all track [Fe/H] from the
                given suite are used.

            model_suite:    The software suite used to generate the stellar
                evolution tracks. May be omitted if tracks are specified by
                filename, but must be supplied if using masses and [Fe/H].

            new_interp_name:    Name to assign to the a newly generated
                interolator. Ignored if an interpolator matching all other
                arguments already exists. If not specified, and no interpolator
                exists matching the remining arguments, a new interpolator is
                not generated.

            num_threads:    If a new interpolator is created this many
                simultaneous interpolation threads are used.

        Returns:
            VarChangingInterpolator:
                Configured per the arguments supplied or None if
                no existing interpolator is found and creating a new one is
                forbidden (see new_interp_name argument).
        """

        with db_session_scope() as db_session:
            if track_fnames is None:
                track_grid = self._track_grid_from_grid(masses,
                                                        feh,
                                                        model_suite,
                                                        db_session)
            else:
                track_grid = self._track_grid_from_files(track_fnames,
                                                         db_session,
                                                         model_suite)

            result = self._find_existing_interpolator(
                track_grid=track_grid,
                nodes=nodes,
                smoothing=smoothing,
                vs_log_age=vs_log_age,
                log_quantity=log_quantity,
                db_session=db_session
            )
            if result is not None:
                return result
            if new_interp_name is None:
                return None

            return self._create_new_interpolator(
                track_grid=track_grid,
                nodes=nodes,
                smoothing=smoothing,
                vs_log_age=vs_log_age,
                log_quantity=log_quantity,
                db_session=db_session,
                name=new_interp_name,
                num_threads=num_threads
            )

    def get_interpolator_by_name(self, name):
        """Return the interpolator with the given name."""

        with db_session_scope() as db_session:
            return ManagedInterpolator(
                db_interpolator=db_session.query(
                    SerializedInterpolator
                ).filter_by(
                    name=name
                ).one(),
                serialization_path=self._serialization_path,
                db_session=db_session
            )

    @staticmethod
    def list_interpolator_names():
        """Return a list of all intorpolator names."""

        with db_session_scope() as db_session:
            return [record[0] for record in
                    db_session.query(SerializedInterpolator.name).all()]

    @staticmethod
    def list_suites():
        """Return a list of all software suites with available tracks."""

        with db_session_scope() as db_session:
            return [record[0] for record in
                    db_session.query(ModelSuite.name).all()]

    @staticmethod
    def get_suite_tracks(model_suite='MESA'):
        """Return all tracks from a given suite."""

        with db_session_scope() as db_session:
            return db_session.query(Track).filter(
                and_(Track.model_suite_id == ModelSuite.id,
                     ModelSuite.name == model_suite)
            )

    def register_track(self,
                       track_fname,
                       mass,
                       feh,
                       model_suite='MESA'):
        """Register a track for use in creating interpolators."""

        with db_session_scope() as db_session:
            self._add_track(track_fname,
                            self._get_decimal(mass),
                            self._get_decimal(feh),
                            model_suite,
                            db_session)

    def register_track_collection(self,
                                  track_fnames,
                                  fname_rex=library_track_fname_rex,
                                  model_suite='MESA'):
        """
        Add a collection of tracks with [Fe/H] and M* encoded in filename.

        Args:
            -track_fnames:
                The filenames of the tracks to add.
            - fname_rex:
                A regular expression defining groups named 'MASS' and
                either 'Z' or 'FeH' used to parse the filename for the
                stellar mass and metallicity each track applies to.
            - model_suite:
                The software suite used to generate the stellar evolution
                tracks.
        """

        for fname in track_fnames:
            parsed_fname = fname_rex.match(
                os.path.basename(fname)
            ).groupdict()
            mass = self._get_decimal(parsed_fname['MASS'])
            if 'Z' in parsed_fname:
                feh = self._get_decimal(
                    library.feh_from_z(float(parsed_fname['Z']))
                )
            else:
                feh = self._get_decimal(parsed_fname['FeH'])

            with db_session_scope() as db_session:
                self._add_track(fname,
                                mass,
                                feh,
                                model_suite,
                                db_session)
