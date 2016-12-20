"""Define class for managing many stellar evolution interpolations."""

from .library_interface import MESAInterpolator
from .manager_data_model import\
    DataModelBase,\
    SerializedInterpolator,\
    InterpolationParameters,\
    Quantity,\
    Track,\
    ModelSuite
from . import Session
from . import __path__
from .basic_utils import db_session_scope, Structure, tempdir_scope
from sqlalchemy import exists, and_
import hashlib
import re
import os
import os.path
import shutil
import numpy
from uuid import uuid4 as get_uuid
from decimal import Decimal
import ctypes
from sqlalchemy import create_engine
from glob import glob

class StellarEvolutionManager :
    """
    Class for managing a collection of stellar evolution inteprolations.
    """

    _solarZ = 0.015

    def _get_decimal(self, value) :
        """Prepare a value for db storage as a limited precision decimal."""

        return Decimal(value).quantize(Decimal('0.0001'))

    default_track_dir = os.path.abspath(os.path.join(__path__[0],
                                                     '../../MESA_tracks'))

    def _define_evolution_quantities(self, db_session) :
        """
        Define the set of quantities tracked by MESAInterpolator instances.

        Args:
            - db_session: The currently active database session.

        Returns: None
        """

        db_quantities = [
            Quantity(id = id, name = name)
            for name, id in MESAInterpolator.quantity_ids.items()
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

    def _get_track_grid(self, track_fnames, fname_rex) :
        """
        Organize checksums for all given files in a mass - [Fe/H] grid.

        Args:
            - track_fnames:
                See get_interpolator track_fnames argument.
            - fname_rex:
                See get_interpolator fname_rex argument.

        Returns:
            A dictionary with keys the masses of tracks and values
            further dictionaries with keys the metallicities of tracks
            and values (filename, checksum) of each track. 
        """

        track_grid = dict()
        for fname in track_fnames :
            parsed_fname = fname_rex.match(os.path.basename(fname))
            mass_key = self._get_decimal(parsed_fname.group('MASS'))
            metallicity_key = self._get_decimal(
                numpy.log10(float(parsed_fname.group('METALLICITY'))
                            /
                            self._solarZ) 
            )
            if mass_key not in track_grid : track_grid[mass_key] = dict()
            assert(metallicity_key not in track_grid[mass_key])
            with open(fname, 'rb') as track_content : 
                track_grid[mass_key][metallicity_key] = (
                    fname,
                    self._checksum(track_content.read()).hexdigest()
                )
        return track_grid

    def _find_existing_interpolator(self,
                                    nodes,
                                    smoothing,
                                    track_grid,
                                    db_session) :
        """
        Return the specified interpolation if already exists, otherwise None.

        Args:
            - nodes:
                see get_interpolator.
            - smoothing:
                see get_interpolator.
            - track_grid:
                See result of _get_track_grid()
            - db_session:
                The currently active database session.

        Returns:
            An instance of MESAInterpolator created from a pre-serialized
            interpolation matching the given arguments if one is found in the
            interpolation archive. If no pre-serialized interpolation exists,
            returns None.
        """

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
                        self._get_decimal(smoothing[quantity.name])
                        ,
                        InterpolationParameters.nodes
                        ==
                        nodes[quantity.name]
                    )
                )
                for quantity in self._quantities
            ]
            +
            [
                SerializedInterpolator.tracks.any(
                    and_(Track.mass == self._get_decimal(mass),
                         Track.metallicity == self._get_decimal(metallicity),
                         Track.checksum == checksum)
                )
                for mass, mass_row in track_grid.items()
                for metallicity, (track_fname, checksum) in mass_row.items()
            ]
        )
        return db_session.query(
            SerializedInterpolator
        ).filter(*match_config).one_or_none()

    def _create_new_interpolator(self,
                                 track_fnames,
                                 nodes,
                                 smoothing,
                                 track_grid,
                                 model_suite,
                                 db_session) :
        """
        Generate the specified interpolation and add it to the archive.

        Args:
            - track_fnames:
                see get_interpolator()
            - nodes:
                see get_interpolator().
            - smoothing:
                see get_interpolator().
            - track_grid:
                See result of _get_track_grid()
            - model_suite:
                The stellar evolution suite used to generate the tracks.

        Returns: 
            An instance of MESAInterpolator created from scratch based on the
            given arguments.
        """

        interp_fname = os.path.join(self._serialization_path, 
                                    str(get_uuid()))

        db_interpolator = SerializedInterpolator(filename = interp_fname)

        with tempdir_scope() as track_dir :
            for mass, mass_row in track_grid.items() :
                for metallicity, (track_fname, checksum) in mass_row.items():
                    track = db_session.query(Track).filter_by(
                        mass = mass,
                        metallicity = metallicity,
                        checksum = checksum
                    ).one_or_none()
                    if track is None :
                        track = Track(id = self._new_track_id,
                                      mass = mass,
                                      metallicity = metallicity,
                                      checksum = checksum,
                                      suite = model_suite)

                        print(70*'-')
                        print('TRACK: ' + repr(track))
                        print(70*'-')

                        db_session.add(track)
                        self._new_track_id += 1
                    db_interpolator.tracks.append(track)
                    shutil.copy(track_fname, track_dir)
            interp_smoothing = numpy.empty(
                len(MESAInterpolator.quantity_list),
                dtype = ctypes.c_double
            )
            interp_nodes = numpy.empty(
                len(MESAInterpolator.quantity_list),
                dtype = ctypes.c_int
            )
            for q_name, q_index in MESAInterpolator.quantity_ids.items() :
                interp_smoothing[q_index] = smoothing[q_name]
                interp_nodes[q_index] = nodes[q_name]
            actual_interpolator = MESAInterpolator(
                mesa_dir = track_dir,
                smoothing = interp_smoothing,
                nodes = interp_nodes
            )
            actual_interpolator.save(interp_fname)

        db_interpolator.parameters = [
            InterpolationParameters(quantity_id = q.id,
                                    nodes = nodes[q.name],
                                    smoothing = smoothing[q.name],
                                    interpolator = db_interpolator)
            for q in self._quantities
        ]

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
        self._checksum = hashlib.sha1
        self._serialization_path = serialization_path
        with db_session_scope() as db_session :
            self._initialize_database(db_engine, db_session)
        with db_session_scope() as db_session :
            self._get_db_config(db_session)

    def get_interpolator(self,
                         track_fnames = glob(os.path.join(default_track_dir,
                                                          '*.csv')),
                         nodes = MESAInterpolator.default_nodes,
                         smoothing = MESAInterpolator.default_smoothing,
                         fname_rex = re.compile(
                             'M(?P<MASS>[0-9.E+-]+)'
                             '_'
                             'Z(?P<METALLICITY>[0-9.E+-]+)'
                             '.csv'
                         ),
                         model_suite = 'MESA',
                         masses = None,
                         metallicities = None) :
        """
        Return a stellar evolution interpolator with the given configuration.

        Args:
            - track_fnames:
                A list of filenames to base the interpolation on. The base
                name should follow the format specified by the fname_rex
                argument to allow extracting the mass and metallicity it
                corresponds to. Must not be given simultaneously with masses
                and metallicities arguments.
            - nodes:
                The number of nodes to use for the age interpolation of each
                quantity of each track. Should be a dictionary with keys
                MESAInterpolator.quantity_list. See the POET code
                StellarEvolution::Interpolator::create_from() documentation
                for a description of what this actually means.
            - smoothing:
                The amount of smoothing to use for the age interpolation of
                each quantity of each track. Should be a dictionary with keys
                MESAInterpolator.quantity_list. See the POET code
                StellarEvolution::Interpolator::create_from() documentation
                for a description of what this actually means.
            - fname_rex:
                A regular expression defining groups named 'MASS' and
                'METALLICITY' used to parse the filename for the stellar mass
                and metallicity each track applies to.
            - model_suite:
                The software suite used to generate the stellar evolution
                tracks.
            - masses:
                A list of the stellar masses to include in the interpolation.
                Unique tracks with those masses and all selected
                metallicities (see next argument) must already be registered
                with the database for the given suite. Must not be supplied
                simultaneously with track_fnames argument.
            - metallicities:
                A list of the stellar metallicities to include in the
                interpolation. Must not be supplied simultaneously with
                track_fnames argument.

        Returns:
            An instance of MESAInterpolator (see stellar_evolution python
            module) configured per the arguments supplied.
        """

        if type(fname_rex) is str : fname_rex = re.compile(fname_rex)
        track_grid = self._get_track_grid(track_fnames, fname_rex)
        with db_session_scope() as db_session :
            return (
                self._find_existing_interpolator(nodes,
                                                 smoothing,
                                                 track_grid,
                                                 db_session)
                or
                self._create_new_interpolator(
                    track_fnames,
                    nodes,
                    smoothing,
                    track_grid,
                    db_session.query(ModelSuite).filter_by(
                        name = model_suite
                    ).one(),
                    db_session
                )
            )


