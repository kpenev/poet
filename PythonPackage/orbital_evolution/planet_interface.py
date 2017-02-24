#!/usr/bin/python3 -u

"""An interface to the POET planet library."""

from ctypes import cdll, c_void_p, c_double
from ctypes.util import find_library

class c_planet_p(c_void_p) : pass

def initialize_library() :
    """Prepare the planet library for use."""

    library = cdll.LoadLibrary(find_library('planet'))

    library.create_planet.argtypes = [c_double, c_double]
    library.create_planet.restype = c_planet_p

    library.destroy_planet.argtypes = [library.create_planet.restype]
    library.destroy_planet.restype = None

    return library

library = initialize_library()

class LockedPlanet :
    """A class for tidally locked and thus non-dissipative planets."""

    def __init__(self, mass, radius) :
        """
        Create a planet with the given mass and radius.

        Args:
            - mass:
                The mass of the planet in Jovian masses.
            - radius:
                The radius of the planet in Jovian radii. 

        Returns: None
        """

        self.mass = mass
        self.radius = radius
        self.planet = library.create_planet(mass, radius)

    def delete(self) :
        """
        Destroy the library planet created at construction.
        """

        library.destroy_planet(self.planet)

if __name__ == '__main__' :
    planet = LockedPlanet(1.0, 1.0)
