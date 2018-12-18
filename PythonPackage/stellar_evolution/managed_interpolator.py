"""Define a stellar evolution interpolator class managed through a database."""

import os
import hashlib
import sys

import numpy

sys.path.append('..')

#Need to add POET package to module search path before importing
#pylint: disable=wrong-import-position
from stellar_evolution.basic_utils import db_session_scope
from stellar_evolution.change_variables import VarChangingInterpolator

from .manager_data_model import \
    VarchangeGrid,\
    VarchangeDependentValue,\
    VarchangeAgeNode,\
    VarchangeMassNode,\
    VarchangeFeHNode,\
    VarchangeDependentVariable
#pylint: enable=wrong-import-position

def checksum_filename(fname):
    """Return a str checksum of the file with the given name."""

    assert os.path.exists(fname)
    with open(fname, 'rb') as opened_file:
        return hashlib.sha1(opened_file.read()).hexdigest()

def verify_checksum(filename, checksum, what):
    """
    Check if the given file has the expected checksum.

    Args:
        filename:    The name of the file whose checksum to verify.

        checksum:    The expected value of the checksum.

        what:    What is being verified (only used if error message if checksums
            do not match).

    Returns:
        None
    """

    if checksum != checksum_filename(filename):
        raise IOError(
            '%s with filename %s registered, with a different checksum!'
            %
            (what.title(), repr(filename))
        )

#This may be something to try to pick apart at a later time.
#pylint: disable=too-many-instance-attributes
class ManagedInterpolator(VarChangingInterpolator):
    """Add properties describing the configuration of an interpolator."""

    def _new_var_change_grid(self,
                             *,
                             grid_name,
                             feh,
                             masses,
                             ages,
                             db_session):
        """
        Create a new grid with the given nodes and register it with the DB.

        Args:
            grid_name:    The name to assign to the new grid in the database.

            feh:    The [Fe/H] values at which to tabulate the dependent
                variables.

            masses:    The stellar masses at which to tabulate the dependent
                variables.

            ages:    The ages (in Gyrs) at which to tabulate the dependent
                variables.

            db_session:    A database session to submit queries to.

        Returns:
            None
        """

        self._varchange_grid_name = grid_name
        self._define_var_change_grid(feh, masses, ages)

        self._grid_db_id = db_session.query(VarchangeGrid).count() + 1
        db_grid = VarchangeGrid(
            id=self._grid_db_id,
            name=grid_name,
            interpolator_id=self._db_id
        )
        db_grid.feh_nodes = [
            VarchangeFeHNode(index=index, value=value)
            for index, value in enumerate(feh)
        ]
        db_grid.mass_nodes = [
            VarchangeMassNode(index=index, value=value)
            for index, value in enumerate(masses)
        ]
        db_grid.age_nodes = [
            VarchangeAgeNode(index=index, value=value)
            for index, value in enumerate(ages)
        ]

        db_session.add(db_grid)
        db_session.add_all(db_grid.feh_nodes)
        db_session.add_all(db_grid.mass_nodes)
        db_session.add_all(db_grid.age_nodes)

    @staticmethod
    def _variable_db_id(variable, db_session, must_exist=True):
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
        ).filter_by(name=variable).all()

        if must_exist or variable_db_id:
            return variable_db_id[0][0]

        variable_db_id = (
            db_session.query(VarchangeDependentVariable).count() + 1
        )
        db_session.add(VarchangeDependentVariable(id=variable_db_id,
                                                  name=variable))
        return variable_db_id

    def _read_variable_from_db(self, variable, db_session):
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
            numpy.empty(
                (
                    #Fales positive
                    #pylint: disable=no-member
                    self.grid.masses.size,
                    self.grid.ages.size,
                    self.grid.feh.size
                    #pylint: enable=no-member
                ),
                dtype=(bool if variable == 'weights' else float)
            )
        )
        grid_var = getattr(self.grid, variable)
        for feh_i, mass_i, age_i, value in db_session.query(
                VarchangeDependentValue.feh_node_index,
                VarchangeDependentValue.mass_node_index,
                VarchangeDependentValue.age_node_index,
                VarchangeDependentValue.value
        ).filter_by(
            variable_id=variable_db_id,
            grid_id=self._grid_db_id
        ):
            grid_var[mass_i, age_i, feh_i] = value

    def _add_variable_to_db(self, variable, db_session):
        """Add pre-calculated node values of a variable to DB."""

        variable_db_id = self._variable_db_id(variable, db_session, False)
        grid_var = getattr(self.grid, variable)
        db_session.add_all(
            (
                VarchangeDependentValue(
                    variable_id=variable_db_id,
                    grid_id=self._grid_db_id,
                    feh_node_index=feh_i,
                    mass_node_index=mass_i,
                    age_node_index=age_i,
                    value=grid_var[mass_i, age_i, feh_i]
                )
                #Fales positive
                #pylint: disable=no-member
                for feh_i in range(self.grid.feh.size)
                for mass_i in range(self.grid.masses.size)
                for age_i in range(self.grid.ages.size)
                #pylint: enable=no-member
            )
        )

    def _add_grid_variable(self, variable):
        """
        Prepares to use another dependent variable to change from.

        Args: see VarChangingInterpolator._add_grid_variable.

        Returns: None
        """

        with db_session_scope() as db_session:
            variable_db_id = self._variable_db_id(variable,
                                                  db_session,
                                                  False)
            if (
                    db_session.query(
                        VarchangeDependentValue
                    ).filter_by(
                        variable_id=variable_db_id,
                        grid_id=self._grid_db_id
                    ).count() > 0
            ):
                self._read_variable_from_db(variable, db_session)
                #Attribute defined_weights defined by parent class.
                #pylint: disable=access-member-before-definition
                #pylint: disable=attribute-defined-outside-init
                if not self.defined_weights:
                    self._read_variable_from_db('weights', db_session)
                    self.defined_weights = True
                #pylint: enable=access-member-before-definition
            else:
                new_weights = not self.defined_weights
                super()._add_grid_variable(variable)
                assert self.defined_weights
                self._add_variable_to_db(variable, db_session)
                if new_weights:
                    self._add_variable_to_db('weights', db_session)

    def _set_var_change_grid(self, grid_name, db_session):
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
            name=grid_name,
            interpolator_id=self._db_id
        ).all()
        if not grid_db_id:
            return False

        self._grid_db_id = grid_db_id[0][0]
        self._define_var_change_grid(
            feh=numpy.array(
                db_session.query(
                    VarchangeFeHNode.value
                ).filter_by(
                    grid_id=self._grid_db_id,
                ).order_by(
                    VarchangeFeHNode.index
                ).all()
            ).flatten(),
            masses=numpy.array(
                db_session.query(
                    VarchangeMassNode.value
                ).filter_by(
                    grid_id=self._grid_db_id,
                ).order_by(
                    VarchangeMassNode.index
                ).all()
            ).flatten(),
            ages=numpy.array(
                db_session.query(
                    VarchangeAgeNode.value
                ).filter_by(
                    grid_id=self._grid_db_id,
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
                 **kwargs):
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

            - track_feh:
                List of stellar [Fe/H] on whose tracks the interpolation is
                based.

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

        if db_interpolator.checksum is not None:
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
        assert len(suite) == 1
        self.suite = suite.pop()
        for param in db_interpolator.parameters:
            quantity = self.quantity_names[param.quantity_id]
            self.smoothing[quantity] = param.smoothing or float('nan')
            self.nodes[quantity] = param.nodes
            self.vs_log_age[quantity] = param.vs_log_age
            self.log_quantity[quantity] = param.log_quantity
        self.track_masses = sorted(
            {track.mass  for track in db_interpolator.tracks}
        )
        self.track_feh = sorted(
            {track.feh for track in db_interpolator.tracks}
        )
        if kwargs:
            super().__init__(grid_feh=numpy.array([]),
                             grid_masses=numpy.array([]),
                             grid_ages=numpy.array([]),
                             **kwargs)
        else:
            super().__init__(grid_feh=numpy.array([]),
                             grid_masses=numpy.array([]),
                             grid_ages=numpy.array([]),
                             interpolator_fname=interpolator_fname)
        if not self._set_var_change_grid('default', db_session):
            self._new_var_change_grid(
                grid_name='default',
                feh=numpy.linspace(
                    float(self.track_feh[0]) * 0.99,
                    float(self.track_feh[-1]) * 0.99,
                    3 * len(self.track_feh)
                ),
                masses=numpy.linspace(
                    float(self.track_masses[0]),
                    float(self.track_masses[-1]),
                    10 * len(self.track_masses)
                ),
                ages=numpy.linspace(1e-2, 13.71, 412),
                db_session=db_session
            )

    def __str__(self):
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
                '], [Fe/H]: ['
                +
                ', '.join([str(feh) for feh in self.track_feh])
                +
                ']')
#pylint: enable=too-many-instance-attributes
