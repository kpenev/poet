#!/usr/bin/env python3

"""An interface to the POET planet library."""

from orbital_evolution.evolve_interface import\
    c_dissipating_body_p,\
    DissipatingBody
from orbital_evolution.evolve_interface import\
    library as orbital_evolution_library
from ctypes import cdll, c_void_p, c_double
from ctypes.util import find_library
from astropy import units, constants

def initialize_library() :
    """Prepare the planet library for use."""

    library_fname = find_library('planet')
    if(library_fname is None) :
        raise OSError('Unable to find POET\'s planet library.') 
    library = cdll.LoadLibrary(library_fname)

    library.create_planet.argtypes = [c_double, c_double]
    library.create_planet.restype = c_dissipating_body_p

    library.destroy_planet.argtypes = [library.create_planet.restype]
    library.destroy_planet.restype = None

    return library

library = initialize_library()

class LockedPlanet(DissipatingBody) :
    """A class for tidally locked and thus non-dissipative planets."""

    lib_configure_body = orbital_evolution_library.configure_planet

    def __init__(self, mass, radius) :
        """
        Create a planet with the given mass and radius.

        Args:
            - mass:
                The mass of the planet in solar masses.
            - radius:
                The radius of the planet in solar radii. 

        Returns: None
        """

        self.mass = mass
        self.radius = radius
        self.c_body = library.create_planet(
            (mass * constants.M_sun / constants.M_jup).to('').value,
            (radius * constants.R_sun / constants.R_jup).to('').value
        )

    def delete(self) :
        """
        Destroy the library planet created at construction.
        """

        library.destroy_planet(self.c_body)

if __name__ == '__main__' :
    planet = LockedPlanet((constants.M_jup / constants.M_sun).to(''),
                          (constants.R_jup / constants.R_sun).to(''))
