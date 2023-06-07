#!/usr/bin/env python3

"""An interface to the poet single tidal term body libray."""

from ctypes import cdll, c_double, c_int
from ctypes.util import find_library

from orbital_evolution.evolve_interface import\
    library as orbital_evolution_library
from orbital_evolution.dissipating_body import\
    c_dissipating_body_p,\
    DissipatingBody

def initialize_library(library_fname=None):
    """Prepare the planet library for use and return it."""

    if library_fname is None:
        library_fname = find_library('singleTermNonEvolvingBody')
    if library_fname is None:
        raise OSError('Unable to find POET\'s single term objecct library.')
    result = cdll.LoadLibrary(library_fname)

    result.create_single_term_non_evolving_body.argtypes = [c_double, c_double]
    result.create_single_term_non_evolving_body.restype = c_dissipating_body_p

    result.destroy_single_term_non_evolving_body.argtypes = [
        result.create_single_term_non_evolving_body.restype
    ]
    result.destroy_single_term_non_evolving_body.restype = None

    result.set_single_term_non_evolving_body_dissipation.argtypes = [
        result.create_single_term_non_evolving_body.restype,#the body
        c_int,      #obit frequency multiplier
        c_int,      #spin_frequency_multiplier
        c_double,   #the phase lag of the only dissipative term
    ]
    result.get_single_term_non_evolving_body_inertia.argtyes = [
        result.create_single_term_non_evolving_body.restype
    ]
    result.get_single_term_non_evolving_body_inertia.restype = c_double

    return result

library = initialize_library()

class SingleTermBody(DissipatingBody):
    """A class for object with single dissipative tidal term."""

    lib_configure_body = orbital_evolution_library.configure_single_term_body

    def __init__(self, mass, radius):
        """
        Create an object with a single dissipative tidal term.

        Args:
            - mass:
                The mass of the planet in solar masses.

            - radius:
                The radius of the planet in solar radii.

        Returns: None
        """

        super().__init__()
        self.mass = mass
        self.radius = radius
        self.c_body = library.create_single_term_non_evolving_body(mass, radius)

    def delete(self):
        """
        Destroy the library planet created at construction.
        """

        library.destroy_single_term_non_evolving_body(self.c_body)

    #The parent method simply saves the parameters, so it need not name them.
    #pylint: disable=arguments-differ
    def set_dissipation(self,
                        orbit_frequency_multiplier,
                        spin_frequency_multiplier,
                        phase_lag):
        """Set the only dissipative term of the object."""

        library.set_single_term_non_evolving_body_dissipation(
            self.c_body,
            orbit_frequency_multiplier,
            spin_frequency_multiplier,
            phase_lag
        )

        super().set_dissipation(
            zone_index=0,
            orbit_frequency_multiplier=orbit_frequency_multiplier,
            spin_frequency_multiplier=spin_frequency_multiplier,
            phase_lag=phase_lag
        )
    #pylint: enable=arguments-differ

    def inertia(self):
        """Return the moment of inertia of the object."""

        return library.get_single_term_non_evolving_body_inertia(self.c_body)


if __name__ == '__main__':
    #False positive.
    #pylint: disable=no-member
    body = SingleTermBody(1.0, 1.0)
    body.set_dissipation(2, 2, 1e-6)
    body.delete()
    #pylint: enable=no-member
