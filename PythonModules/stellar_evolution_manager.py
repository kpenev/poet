#!/usr/bin/python3 -u

"""Define class for managing many stellar evolution interpolations."""

from stellar_evolution import MESAInterpolator
from stellar_evolution_manager_data_model import\
    DataModelBase,\
    SerializedInterpolators,\
    Nodes,\
    Smoothing,\
    Quantities
from SQLAlchemyUtil import Session, db_session_scope
from sqlalchemy import create_engine
import hashlib
import re
import os.path

class StellarEvolutionManager :
    """
    Class for managing a collection of stellar evolution inteprolations.
    """

    _solarZ = 0.015

    def _define_evolution_quantities(self, db_session) :
        """
        Define the set of quantities tracked by MESAInterpolator instances.

        Args:
            - db_session: The currently active database session.

        Returns: None
        """

        db_session.add_all([
            Quantities(id = id, name = name)
            for name, id in MESAInterpolator.quantity_ids.items()
        ])

    def _get_track_grid(self, track_fnames, fname_rex) :
        """
        Organize checksums for all given files in a mass - [Fe/H] grid.

        Args:
            - track_fnames:
                See get_interpolator mesa_track_fnames argument.
            - fname_rex:
                See get_interpolator fname_rex argument.

        Returns:
            A dictionary with keys the masses of tracks and values
            further dictionaries with keys the metallicities of tracks
            and values the checksum of each track. 
        """

        track_grid = dict()
        for fname in track_fnames :
            parsed_fname = fname_rex.match(os.path.basename(fname))
            mass_key = round(float(parsed_fname.group('MASS')), 3)
            metallicity_key = round(
                numpy.log10(float(parsed_fname.group('METALLICITY'))
                            /
                            self._solarZ), 
                3
            )
            if mass_key not in trac_grid : track_grid[mass_key] = dict()
            assert(metallicity_key not in track_grid[mass_key])
            with open(fname, 'rb') as track_content : 
                track_grid[metallicity_key] = self._checksum(
                    track_content
                ).hexdigest()
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

        match_nodes = [
            db_session.query(Nodes).
        ]
        interpolator = db_session.query(SerializedInterpolators).filter(
            

    def _create_new_interpolator(self,
                                 mesa_track_fnames,
                                 nodes,
                                 smoothing,
                                 track_grid) :
        """
        Generate the specified interpolation and add it to the archive.

        Args:
            - mesa_track_fnames:
                see get_interpolator()
            - nodes:
                see get_interpolator().
            - smoothing:
                see get_interpolator().
            - track_grid:
                See result of _get_track_grid()

        Returns: 
            An instance of MESAInterpolator created from scratch based on the
            given arguments.
        """

    def __init__(self, db_engine) :
        """
        Attach the manager to the given archive (sqlite file).

        Args:
            - interpolation_archive:
                The name of an sqlite3 database holding the information about
                the current set of serialized stellar evolution
                interpolations. Created if it does not exist.

        Returns: None.
        """

        self._checksum = hashlib.sha1

        for table in DataModelBase.metadata.sorted_tables :
            if not db_engine.has_table(table.name) :
                table.create(db_engine)
                if table.name == 'quantities' :
                    with db_session_scope() as db_session :
                        self._define_evolution_quantities(db_session)

    def get_interpolator(self,
                         mesa_track_fnames,
                         nodes,
                         smoothing, 
                         fname_rex = re.compile(
                             'M(?P<MASS>[0-9.E+-]+)'
                             '_'
                             'Z(?P<METALLICITY>[0-9.E+-]+)'
                             '.csv'
                         )) :
        """
        Return a stellar evolution interpolator with the given configuration.

        Args:
            - mesa_tracks_fnames:
                A list of filenames to base the interpolation on. The base
                name should follow the format specified by the fname_rex
                argument to allow extracting the mass and metallicity it
                corresponds to.
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

        Returns:
            An instance of MESAInterpolator (see stellar_evolution python
            module) configured per the arguments supplied.
        """

        if type(fname_rex) is str : fname_rex = re.compile(fname_rex)
        track_grid = self._get_track_grid(mesa_track_fnames, fname_rex)
        return (self._find_existing_interpolator(nodes,
                                                 smoothing,
                                                 track_grid)
                or
                self._create_new_interpolator(mesa_track_fnames,
                                              nodes,
                                              smoothing,
                                              track_grid))

if __name__ == '__main__' :
    db_engine = create_engine('sqlite:///test.sqlite', echo = True)
    Session.configure(bind = db_engine)
    manager = StellarEvolutionManager(db_engine)
