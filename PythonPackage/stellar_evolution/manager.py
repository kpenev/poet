"""Define class for managing many stellar evolution interpolations."""

import sys
sys.path.append('..')

from stellar_evolution.change_variables import VarChangingInterpolator
from stellar_evolution.library_interface import\
    library_track_fname_rex,\
    library_track_fname,\
    library
from .manager_data_model import\
    DataModelBase,\
    SerializedInterpolator,\
    InterpolationParameters,\
    Quantity,\
    Track,\
    ModelSuite,\
    VarchangeGrid,\
    VarchangeDependentValue,\
    VarchangeAgeNode,\
    VarchangeMassNode,\
    VarchangeMetallicityNode,\
    VarchangeDependentVariable
from stellar_evolution import Session
from stellar_evolution.basic_utils import db_session_scope, tempdir_scope
from basic_utils import Structure
from sqlalchemy import exists, and_, literal, func, create_engine
import hashlib
import re
import os
import os.path
import shutil
import numpy
from uuid import uuid4 as get_uuid
from decimal import Decimal
import ctypes
from glob import glob
from math import isnan

def checksum_filename(fname) :
    """Return a str checksum of the file with the given name."""

    assert(os.path.exists(fname))
    with open(fname, 'rb') as opened_file :
        return hashlib.sha1(opened_file.read()).hexdigest()

def verify_checksum(filename, checksum, what) :
    """
    Check if the given file has the expected checksum.

    Args:
        - filename:
            The name of the file whose checksum to verify.
        - checksum:
            The expected value of the checksum.
        - what: 
            What is being verified (only used if error message if checksums
            do not match).

    Returns: None
    """

    if checksum != checksum_filename(filename) :
        raise IOError(
            '%s with filename %s registered, with a different checksum!'
            %
            (what.title(), repr(filename))
        )

class ManagedInterpolator(VarChangingInterpolator) :
    """Add properties describing the configuration of an interpolator."""

    def _new_var_change_grid(self,
                             grid_name,
                             metallicities,
                             masses,
                             ages,
                             db_session) :
        """
        Create a new grid with the given nodes and register it with the DB.

        Args:
            - grid_name:
                The name to assign to the new grid in the database.

            - metallicities:
                The [Fe/H] values at which to tabulate the dependent
                variables.

            - masses:
                The stellar masses at which to tabulate the dependent
                variables.

            - ages:
                The ages (in Gyrs) at which to tabulate the dependent
                variables.

            - db_session:
                A database session to submit queries to.

        Returns: None
        """

        self._varchange_grid_name = grid_name
        self._define_var_change_grid(metallicities, masses, ages)

        self._grid_db_id = db_session.query(VarchangeGrid).count() + 1
        db_grid = VarchangeGrid(
            id = self._grid_db_id,
            name = grid_name,
            interpolator_id = self._db_id
        )
        db_grid.metallicity_nodes = [
            VarchangeMetallicityNode(index = index, value = value)
            for index, value in enumerate(metallicities)
        ]
        db_grid.mass_nodes = [
            VarchangeMassNode(index = index, value = value)
            for index, value in enumerate(masses)
        ]
        db_grid.age_nodes = [
            VarchangeAgeNode(index = index, value = value)
            for index, value in enumerate(ages)
        ]

        db_session.add(db_grid)
        db_session.add_all(db_grid.metallicity_nodes)
        db_session.add_all(db_grid.mass_nodes)
        db_session.add_all(db_grid.age_nodes)

    def _variable_db_id(self, variable, db_session, must_exist = True) :
        """
        Return the ID of the given varibale in the database.

        Args:
            - variable:
                The name of the variable whose ID to return.

            - must_exst:
                If False, and the variable is not yet in the
                varchange_dependent_variables table, a new entry is added.
                Otherwise, an exception is raised if it is not there.

        Returns:
            - variable_db_id:
                The ID of the variable in the database.
        """

        variable_db_id = db_session.query(
            VarchangeDependentVariable.id
        ).filter_by(name = variable).all()
        if must_exist or len(variable_db_id) > 0 :
            return variable_db_id[0][0]
        else :
            variable_db_id = (
                db_session.query(VarchangeDependentVariable).count() + 1
            )
            db_session.add(VarchangeDependentVariable(id = variable_db_id,
                                                      name = variable))
            return variable_db_id


    def _read_variable_from_db(self, variable, db_session) :
        """
        Read the given variable's grid values from the DB.

        Args:
            - variable_name:
                The name of the variable to add.

            - db_session:
                A database session for queries.

        Returns: 
            None, but has the same effect as calling
            VarChangingInterpolator._add_grid_variable, but finishes much
            faster.
        """

        variable_db_id = self._variable_db_id(variable, db_session)
        setattr(
            self.grid,
            variable,
            numpy.empty((self.grid.masses.size,
                         self.grid.ages.size,
                         self.grid.metallicities.size),
                        dtype = bool if variable == 'weights' else float)
        )
        grid_var = getattr(self.grid, variable)
        for metallicity_i, mass_i, age_i, value in db_session.query(
                VarchangeDependentValue.metallicity_node_index,
                VarchangeDependentValue.mass_node_index,
                VarchangeDependentValue.age_node_index,
                VarchangeDependentValue.value
        ).filter_by(
            variable_id = variable_db_id,
            grid_id = self._grid_db_id
        ) :
            grid_var[mass_i, age_i, metallicity_i] = value

    def _add_variable_to_db(self, variable, db_session) :
        """Add pre-calculated node values of a variable to DB."""

        variable_db_id = self._variable_db_id(variable, db_session, False)
        grid_var = getattr(self.grid, variable)
        db_session.add_all(
            (
                VarchangeDependentValue(
                    variable_id = variable_db_id,
                    grid_id = self._grid_db_id,
                    metallicity_node_index = metallicity_i,
                    mass_node_index = mass_i,
                    age_node_index = age_i,
                    value = grid_var[mass_i, age_i, metallicity_i]
                )
                for metallicity_i in range(self.grid.metallicities.size)
                for mass_i in range(self.grid.masses.size)
                for age_i in range(self.grid.ages.size)
            )
        )

    def _add_grid_variable(self, variable) :
        """
        Prepares to use another dependent variable to change from.

        Args: see VarChangingInterpolator._add_grid_variable.

        Returns: None
        """

        with db_session_scope() as db_session :
            variable_db_id = self._variable_db_id(variable,
                                                  db_session,
                                                  False)
            if (
                    db_session.query(
                        VarchangeDependentValue
                    ).filter_by(
                        variable_id = variable_db_id,
                        grid_id = self._grid_db_id
                    ).count() > 0
            ) :
                self._read_variable_from_db(variable, db_session)
                if not self._defined_weights :
                    self._read_variable_from_db('weights', db_session)
                    self._defined_weights = True
            else :
                new_weights = not self._defined_weights
                super()._add_grid_variable(variable)
                assert(self._defined_weights)
                self._add_variable_to_db(variable, db_session)
                if new_weights : self._add_variable_to_db('weights',
                                                          db_session)

    def _set_var_change_grid(self, grid_name, db_session) :
        """
        Read a varchange grid from the DB and set the interpolator to use it.

        See VarChangingInterpolator._define_var_change_grid for newly created
        members of self.

        Args:
            - grid_name:
                The name of the grid in the database to read.

            - db_session:
                A database session for queries.

        Returns:
            True if a grid with the given name exists, False otherwise.
        """

        grid_db_id = db_session.query(
            VarchangeGrid.id
        ).filter_by(
            name = grid_name,
            interpolator_id = self._db_id
        ).all()
        if len(grid_db_id) == 0 :
            return False
        else :
            self._grid_db_id = grid_db_id[0][0]
            self._define_var_change_grid(
                metallicities = numpy.array(
                    db_session.query(
                        VarchangeMetallicityNode.value
                    ).filter_by(
                        grid_id = self._grid_db_id,
                    ).order_by(
                        VarchangeMetallicityNode.index
                    ).all()
                ).flatten(),
                masses = numpy.array(
                    db_session.query(
                        VarchangeMassNode.value
                    ).filter_by(
                        grid_id = self._grid_db_id,
                    ).order_by(
                        VarchangeMassNode.index
                    ).all()
                ).flatten(),
                ages = numpy.array(
                    db_session.query(
                        VarchangeAgeNode.value
                    ).filter_by(
                        grid_id = self._grid_db_id,
                    ).order_by(
                        VarchangeAgeNode.index
                    ).all()
                ).flatten()
            )
            return True

    def __init__(self,
                 db_interpolator,
                 serialization_path,
                 db_session,
                 **kwargs) :
        """
        Create VarChangingInterpolator and add properties describing config.

        Defines the following properties containing the information from
        db_interpolator:
            - name:
                The human readable name of the interpolator

            - _db_id:
                The ID of the interpolator in the database.

            - filename:
                The filename from which the interpolator was read.

            - nodes:
                A dictionary indexed by quantity giving the number of
                interpolation nodes used.

            - smoothing:
                Same as nodes but for the smoothing arguments.

            - suite:
                The software suite used to generate the tracks on which the
                interpolator is based.

            - track_masses:
                List of stellar masses on whose tracks the interpolation is
                based.

            - track_metallicities:
                List of stellar metallicities on whose tracks the
                interpolation is based.

        Args:
            - db_interpolator:
                SerializedInterpolator instance from which to initialize
                self.
            - serialization_path:
                The directory where serialized interpolators are stored.

        Keyword only arguments:
            If not an empty dictionary, the underlying interpolator is
            constructed using those instead of the serialized filename.

        Returns: None
        """

        interpolator_fname = os.path.join(serialization_path,
                                          db_interpolator.filename)

        if db_interpolator.checksum is not None :
            verify_checksum(interpolator_fname,
                            db_interpolator.checksum,
                            'interpolator')

        self.name = db_interpolator.name
        self._db_id = db_interpolator.id
        self.filename = db_interpolator.filename
        self.smoothing = dict()
        self.nodes = dict()
        self.vs_log_age = dict()
        self.log_quantity = dict()
        suite = {track.suite.name for track in db_interpolator.tracks}
        assert(len(suite) == 1)
        self.suite = suite.pop()
        for param in db_interpolator.parameters :
            quantity = self.quantity_names[param.quantity_id]
            self.smoothing[quantity] = param.smoothing or float('nan')
            self.nodes[quantity] = param.nodes
            self.vs_log_age[quantity] = param.vs_log_age
            self.log_quantity[quantity] = param.log_quantity
        self.track_masses = sorted(
            {track.mass  for track in db_interpolator.tracks}
        )
        self.track_metallicities = sorted(
            {track.metallicity for track in db_interpolator.tracks}
        )
        if kwargs :
            super().__init__(grid_metallicities = numpy.array([]),
                             grid_masses = numpy.array([]),
                             grid_ages = numpy.array([]),
                             **kwargs)
        else :
            super().__init__(grid_metallicities = numpy.array([]),
                             grid_masses = numpy.array([]),
                             grid_ages = numpy.array([]),
                             interpolator_fname = interpolator_fname)
        if not self._set_var_change_grid('default', db_session) :
            self._new_var_change_grid(
                'default',
                metallicities = numpy.linspace(
                    float(self.track_metallicities[0]) * 0.99,
                    float(self.track_metallicities[-1]) * 0.99,
                    3 * len(self.track_metallicities)
                ),
                masses = numpy.linspace(
                    float(self.track_masses[0]),
                    float(self.track_masses[-1]),
                    10 * len(self.track_masses)
                ),
                ages = numpy.linspace(1e-2, 13.71, 412),
                db_session = db_session
            )

    def __str__(self) :
        """Human readable representation of the interpolator."""

        return (self.name
                +
                '['
                +
                ', '.join(['%s(n: %d, s: %g)' % (quantity,
                                                 self.nodes[quantity],
                                                 self.smoothing[quantity])
                           for quantity in self.quantity_list])
                +
                '] masses: ['
                +
                ', '.join([str(m) for m in self.track_masses])
                +
                '], metallicites: ['
                +
                ', '.join([str(feh) for feh in self.track_metallicities])
                +
                ']')

class StellarEvolutionManager :
    """
    Class for managing a collection of stellar evolution inteprolations.
    """

    def _get_decimal(self, value) :
        """Prepare a value for db storage as a limited precision decimal."""

        return Decimal(value).quantize(Decimal('0.0001'))

    def _define_evolution_quantities(self, db_session) :
        """
        Define the quantities tracked by VarChangingInterpolator instances.

        Args:
            - db_session: The currently active database session.

        Returns: None
        """

        db_quantities = [
            Quantity(id = id, name = name)
            for name, id in VarChangingInterpolator.quantity_ids.items()
        ]
        db_session.add_all(db_quantities)

    def _initialize_database(self, db_engine, db_session) :
        """
        Ensure all database tables exist and contain at least required data.

        Args:
            - db_session: An sqlalchemy session, used to update the database.

        Returns: None
        """

        for table in DataModelBase.metadata.sorted_tables :
            if not db_engine.has_table(table.name) :
                table.create(db_engine)
                if table.name == 'quantities' :
                    self._define_evolution_quantities(db_session)
                elif table.name == 'model_suites' :
                    db_session.add(ModelSuite(name = 'MESA'))

    def _get_db_config(self, db_session) :
        """Read some configuration from the database."""

        self._quantities = [
            Structure(id = q.id, name = q.name)
            for q in db_session.query(Quantity).all()
        ]

        self._suites = db_session.query(ModelSuite).all()
        self._new_track_id = db_session.query(Track).count() + 1

    def _add_track(self,
                   track_fname,
                   mass,
                   metallicity,
                   model_suite,
                   db_session) :
        """Add a track to the database."""

        db_track = Track(id = self._new_track_id,
                         filename = os.path.abspath(track_fname),
                         mass = mass,
                         metallicity = metallicity,
                         checksum = checksum_filename(track_fname),
                         suite = db_session.query(
                             ModelSuite
                         ).filter_by(
                             name = model_suite
                         ).one())
        db_session.add(db_track)
        self._new_track_id += 1

    def _track_grid_from_files(self,
                               track_fnames,
                               db_session,
                               model_suite = None) :
        """
        Organize track files in a mass - [Fe/H] grid and verify checksums.

        In addition, if model_suite is not None, all tracks are verified to
        belong to this suite.

        Args:
            - track_fnames:
                See get_interpolator track_fnames argument.

        Returns:
            A dictionary with keys the masses of tracks and values
            further dictionaries with keys the metallicities of tracks
            and values the filename of each track.
        """

        track_grid = dict()
        for fname in track_fnames :
            absolute_fname = os.path.abspath(fname)
            db_track = db_session.query(Track).filter_by(
                filename = absolute_fname
            )

            if model_suite is not None :
                assert(db_track.suite.name == model_suite)
            verify_checksum(fname, db_track.checksum, 'track')

            mass_key = db_track.mass
            metallicity_key = db_track.metallicity
            if mass_key not in track_grid : track_grid[mass_key] = dict()
            assert(metallicity_key not in track_grid[mass_key])
            track_grid[mass_key][metallicity_key] = (fname, db_track.id)
        return track_grid

    def _track_grid_from_grid(self,
                              masses,
                              metallicities,
                              model_suite,
                              db_session) :
        """
        Return a mass - [Fe/H] grid with filenames and checksums.

        Fails if multiple tracks are registered for some (mass, metallicity,
        model suite) combination.

        Args:
            - masses:
                The masses for which to include tracks.
            - metallicities:
                The metallicities for which to include tracks.
            - model_suite:
                The software suite from whose tracks to choose.
            - db_session:
                A database session to submit queries to.

        Returns:
            The structure as _track_grid_from_files.
        """

        if masses is None :
            masses = [
                record[0]
                for record in db_session.query(Track.mass).filter(
                    Track.model_suite_id == ModelSuite.id,
                    ModelSuite.name == model_suite
                ).distinct().order_by(Track.mass).all()
            ]
        if metallicities is None :
            metallicities = [
                record[0]
                for record in db_session.query(Track.metallicity).filter(
                    Track.model_suite_id == ModelSuite.id,
                    ModelSuite.name == model_suite
                ).distinct().order_by(Track.metallicity).all()
            ]
        track_grid = {m: dict() for m in masses}
        for m in masses :
            for feh in metallicities :
                db_track = db_session.query(Track).filter(
                    Track.mass == self._get_decimal(m),
                    Track.metallicity == self._get_decimal(feh),
                    Track.model_suite_id == ModelSuite.id,
                    ModelSuite.name == model_suite
                ).one()

                verify_checksum(db_track.filename,
                                db_track.checksum,
                                'track')

                track_grid[m][feh] = (db_track.filename, db_track.id)
        return track_grid

    def _find_existing_interpolator(self,
                                    track_grid,
                                    nodes,
                                    smoothing,
                                    vs_log_age,
                                    log_quantity,
                                    db_session) :
        """
        Return the specified interpolation if already exists, otherwise None.

        Args:
            - track_grid:
                See result of _track_grid_from_files()
            - nodes:
                see get_interpolator.
            - smoothing:
                see get_interpolator.
            - vs_log_age:
                see get_interpolator.
            - log_quantity:
                see get_interpolator.
            - db_session:
                The currently active database session.

        Returns:
            An instance of VarChangingInterpolator created from a
            pre-serialized interpolation matching the given arguments if one
            is found in the interpolation archive. If no pre-serialized
            interpolation exists, returns None.
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
                for metallicity, (track_fname, track_id) in mass_row.items()
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
        if result is None : 
            return result
        else :
            return ManagedInterpolator(
                db_interpolator = result,
                serialization_path = self._serialization_path,
                db_session = db_session
            )

    def _create_new_interpolator(self,
                                 track_grid,
                                 nodes,
                                 smoothing,
                                 vs_log_age,
                                 log_quantity,
                                 num_threads,
                                 db_session,
                                 name = None) :
        """
        Generate the specified interpolation and add it to the archive.

        Args:
            - track_grid:
                See result of _track_grid_from_files()
            - nodes:
                see get_interpolator().
            - smoothing:
                see get_interpolator().
            - vs_log_age:
                 see get_interpolator().
            - log_quantity:
                 see get_interpolator().
            - num_threads:
                The number of simultaneous threads to use when constructing
                the interpolation.
            - db_session:
                The database query session to use.
            - name:
                The name to assign to the new interpolator. If None, the UUID
                used to form the filename is used.

        Returns: 
            An instance of VarChangingInterpolator created from scratch
            based on the given arguments.
        """

        interp_str = str(get_uuid())
        interp_fname = os.path.join(self._serialization_path, 
                                    interp_str)

        db_interpolator = SerializedInterpolator(
            id = db_session.query(SerializedInterpolator).count() + 1,
            name = name or interp_str,
            filename = interp_str
        )

        if (
                db_session.query(
                    SerializedInterpolator
                ).filter_by(name = name).count()
        ) : 
            raise ValueError('Interpolator named %s already exists, with a '
                             'different configuration than the one being '
                             'constructed!' % repr(name))

        with tempdir_scope() as track_dir :
            for mass, mass_row in track_grid.items() :
                for metallicity, (track_fname, track_id) in mass_row.items():
                    track = db_session.query(
                        Track
                    ).filter_by(
                        id = track_id
                    ).one()
                    db_interpolator.tracks.append(track)
                    shutil.copy(
                        track_fname, 
                        os.path.join(track_dir,
                                     library_track_fname(mass, metallicity))
                    )
            interp_smoothing = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype = ctypes.c_double
            )
            interp_nodes = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype = ctypes.c_int
            )
            interp_vs_log_age = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype = ctypes.c_bool
            )
            interp_log_quantity = numpy.empty(
                len(VarChangingInterpolator.quantity_list),
                dtype = ctypes.c_bool
            )
            for q_name, q_index in \
                    VarChangingInterpolator.quantity_ids.items() :
                interp_smoothing[q_index] = smoothing[q_name]
                interp_nodes[q_index] = nodes[q_name]
                interp_vs_log_age[q_index] = vs_log_age[q_name]
                interp_log_quantity[q_index] = log_quantity[q_name]

            db_interpolator.parameters = [
                InterpolationParameters(quantity_id = q.id,
                                        nodes = nodes[q.name],
                                        smoothing = smoothing[q.name],
                                        vs_log_age = vs_log_age[q.name],
                                        log_quantity = log_quantity[q.name],
                                        interpolator = db_interpolator)
                for q in self._quantities
            ]

            actual_interpolator = ManagedInterpolator(
                db_interpolator = db_interpolator,
                serialization_path = self._serialization_path,
                db_session = db_session,
                mesa_dir = track_dir,
                smoothing = interp_smoothing,
                nodes = interp_nodes,
                vs_log_age = interp_vs_log_age,
                log_quantity = interp_log_quantity,
                num_threads = num_threads
            )

        actual_interpolator.save(interp_fname)
        db_interpolator.checksum = checksum_filename(interp_fname)
        db_session.add(db_interpolator)
        db_session.add_all(db_interpolator.parameters)

        return actual_interpolator

    def __init__(self, serialization_path) :
        """
        Create a manager storing serialized interpolators in the given path.

        Args:
            - serialization_path: The path where to store serialized
                                  interpolators.

        Returns: None.
        """

        if not os.path.exists(serialization_path) :
            os.makedirs(serialization_path)
        db_engine = create_engine(
            'sqlite:///'
            + 
            os.path.join(serialization_path, 'serialized.sqlite'),
            echo = True
        )
        Session.configure(bind = db_engine)
        self._serialization_path = serialization_path
        with db_session_scope() as db_session :
            self._initialize_database(db_engine, db_session)
        with db_session_scope() as db_session :
            self._get_db_config(db_session)

    def get_interpolator(
            self,
            nodes = VarChangingInterpolator.default_nodes,
            smoothing = VarChangingInterpolator.default_smoothing,
            vs_log_age = VarChangingInterpolator.default_vs_log_age,
            log_quantity = VarChangingInterpolator.default_log_quantity,
            track_fnames = None,
            masses = None,
            metallicities = None,
            model_suite = 'MESA',
            new_interp_name = None,
            num_threads = 1
    ) :
        """
        Return a stellar evolution interpolator with the given configuration.

        All tracks that the interpolator should be based on must be
        pre-registered with the manager. Two ways are supported for
        identifying tracks: as a list of filenames or as a mass-metallicity
        grid combined with a suite. The first case always works, while the
        second requires that the set of identified tracks is unique, i.e. for
        none of the mass - metallicity combinations there are two or more
        tracks registered for the given suite.

        Args:
            - nodes:
                The number of nodes to use for the age interpolation of each
                quantity of each track. Should be a dictionary with keys
                VarChangingInterpolator.quantity_list. See the POET code
                StellarEvolution::Interpolator::create_from() documentation
                for a description of what this actually means.
            - smoothing:
                The amount of smoothing to use for the age interpolation of
                each quantity of each track. Should be a dictionary with keys
                VarChangingInterpolator.quantity_list. See the POET code
                StellarEvolution::Interpolator::create_from() documentation
                for a description of what this actually means.
            - vs_log_age:
                Use log(age) instead of age as the independent argument for
                the intperpolation? Should be a dictionary with keys
                VarChangingInterpolator.quantity_list.
            - log_quantity:
                Interpolate log(quantity) instead of quantity? Should be a
                dictionary with keys VarChangingInterpolator.quantity_list.
            - track_fnames:
                A list of files containing stellar evolution tracks the
                interpolator should be based on.
            - masses:
                A list of the stellar masses to include in the interpolation.
                Unique tracks with those masses and all selected
                metallicities (see next argument) must already be registered
                with the database for the given suite. If None, all track
                masses from the given suite are used.
            - metallicities:
                A list of the stellar metallicities to include in the
                interpolation. If None, all track metallicities from the
                given suite are used.
            - model_suite:
                The software suite used to generate the stellar evolution
                tracks. May be omitted if tracks are specified by filename,
                but must be supplied if using masses and metallicities.
            - new_interp_name:
                Name to assign to the a newly generated interolator. Ignored
                if an interpolator matching all other arguments already
                exists. If not specified, and no interpolator exists matching
                the remining arguments, a new interpolator is not generated.
            - num_threads:
                If a new interpolator is created this many simultaneous
                interpolation threads are used.

        Returns:
            An instance of VarChangingInterpolator (see stellar_evolution
            python module) configured per the arguments supplied or None if
            no existing interpolator is found and creating a new one is
            forbidden (see new_interp_name argument).
        """

        with db_session_scope() as db_session :
            if track_fnames is None :
                track_grid = self._track_grid_from_grid(masses,
                                                        metallicities,
                                                        model_suite,
                                                        db_session)
            else :
                track_grid = self._track_grid_from_files(track_fnames,
                                                         db_session,
                                                         model_suite)

            result = self._find_existing_interpolator(
                track_grid = track_grid,
                nodes = nodes,
                smoothing = smoothing,
                vs_log_age = vs_log_age,
                log_quantity = log_quantity,
                db_session = db_session
            )
            if result is not None :
                return result
            if new_interp_name is None :
                return None
            else : 
                return self._create_new_interpolator(
                    track_grid = track_grid,
                    nodes = nodes,
                    smoothing = smoothing,
                    vs_log_age = vs_log_age,
                    log_quantity = log_quantity,
                    db_session = db_session,
                    name = new_interp_name,
                    num_threads = num_threads
                )

    def get_interpolator_by_name(self, name) :
        """Return the interpolator with the given name."""

        with db_session_scope() as db_session :
            return ManagedInterpolator(
                db_interpolator = db_session.query(
                    SerializedInterpolator
                ).filter_by(
                    name = name
                ).one(),
                serialization_path = self._serialization_path,
                db_session = db_session
            )

    def list_interpolator_names(self) :
        """Return a list of all intorpolator names."""

        with db_session_scope() as db_session :
            return [record[0] for record in 
                    db_session.query(SerializedInterpolator.name).all()]

    def list_suites(self) :
        """Return a list of all software suites with available tracks."""

        with db_session_scope() as db_session :
            return [record[0] for record in 
                    db_session.query(ModelSuite.name).all()]

    def get_suite_tracks(self, model_suite = 'MESA') :
        """Return all tracks from a given suite."""

        with db_session_scope() as db_session :
            return db_session.query(Track).filter(
                and_(Track.model_suite_id == ModelSuite.id,
                     ModelSuite.name == model_suite)
            )

    def register_track(self,
                       track_fname,
                       mass,
                       metallicity,
                       model_suite = 'MESA') :
        """Register a track for use in creating interpolators."""

        with db_session_scope() as db_session :
            self._add_track(track_fname,
                            self._get_decimal(mass),
                            self._get_decimal(metallicity),
                            model_suite,
                            db_session)

    def register_track_collection(self,
                                  track_fnames,
                                  fname_rex = library_track_fname_rex,
                                  model_suite = 'MESA') :
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

        for fname in track_fnames :
            parsed_fname = fname_rex.match(
                os.path.basename(fname)
            ).groupdict()
            mass = self._get_decimal(parsed_fname['MASS'])
            if 'Z' in parsed_fname :
                metallicity = self._get_decimal(
                    library.feh_from_z(float(parsed_fname['Z']))
                )
            else :
                metallicity = self._get_decimal(parsed_fname['FeH'])

            with db_session_scope() as db_session :
                self._add_track(fname,
                                mass,
                                metallicity,
                                model_suite,
                                db_session)
