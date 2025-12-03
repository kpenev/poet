#!/usr/bin/env python3

"""An interface to the POET planet library."""

from ctypes import cdll, c_double, c_uint
from ctypes.util import find_library

import numpy
from astropy import constants

from orbital_evolution.c_interface_util import ndpointer_or_null
from orbital_evolution.dissipating_body import\
    c_dissipating_body_p,\
    DissipatingBody
from orbital_evolution.evolve_interface import\
    library as orbital_evolution_library

def initialize_library(library_fname=None):
    """Prepare the planet library for use and return it."""

    if library_fname is None:
        library_fname = find_library('planet')
    if library_fname is None:
        raise OSError('Unable to find POET\'s planet library.')
    result = cdll.LoadLibrary(library_fname)

    result.create_planet.argtypes = [c_double, c_double]
    result.create_planet.restype = c_dissipating_body_p

    result.destroy_planet.argtypes = [result.create_planet.restype]
    result.destroy_planet.restype = None

    result.set_planet_dissipation.argtypes = [
        result.create_planet.restype,                   #planet
        c_uint,                                         #num_tidal_freq_breaks
        c_uint,                                         #num_spin_freq_breaks
        c_uint,                                         #num_age_breaks
        ndpointer_or_null(dtype=c_double,               #tidal_freq_breaks
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,               #spin_freq_breaks
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_double,       #tidal_freq_powers
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_double,       #spin_freq_powers
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,       #age_breaks
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_double,       #reference_phase_lags
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_double,                                       #inertial_mode_enhancmnt
        c_double                                        #inertial_mode_sharpness
    ]
    result.get_planet_inertia.argtypes = [result.create_planet.restype]
    result.get_planet_inertia.restype = c_double

    return result

library = initialize_library()

class LockedPlanet(DissipatingBody):
    """A class for tidally locked and thus non-dissipative planets."""

    lib_configure_body = orbital_evolution_library.configure_planet

    def __init__(self, mass, radius):
        """
        Create a planet with the given mass and radius.

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
        self.c_body = library.create_planet(mass, radius)

    def delete(self):
        """
        Destroy the library planet created at construction.
        """

        library.destroy_planet(self.c_body)

    #The parent method simply saves the parameters, so it need not name them.
    #pylint: disable=arguments-differ
    def set_dissipation(self,
                        *,
                        tidal_frequency_breaks,
                        spin_frequency_breaks,
                        tidal_frequency_powers,
                        spin_frequency_powers,
                        age_breaks,
                        reference_phase_lags,
                        inertial_mode_enhancement=1.0,
                        inertial_mode_sharpness=10.0):
        """
        Set the dissipaation of the only zone of the planet.

        See EvolvingStar.set_dissipation() for description of the arguments.

        Returns:
            None
        """

        library.set_planet_dissipation(
            self.c_body,
            tidal_frequency_powers.size - 1,
            spin_frequency_powers.size - 1,
            0 if age_breaks is None else age_breaks.size,
            tidal_frequency_breaks,
            spin_frequency_breaks,
            tidal_frequency_powers,
            spin_frequency_powers,
            age_breaks,
            reference_phase_lags,
            inertial_mode_enhancement,
            inertial_mode_sharpness
        )
        super().set_dissipation(
            zone_index=0,
            tidal_frequency_breaks=tidal_frequency_breaks,
            spin_frequency_breaks=spin_frequency_breaks,
            tidal_frequency_powers=tidal_frequency_powers,
            spin_frequency_powers=spin_frequency_powers,
            reference_phase_lag=reference_phase_lags,
            inertial_mode_enhancement=inertial_mode_enhancement,
            inertial_mode_sharpness=inertial_mode_sharpness
        )


    def inertia(self):
        """Return the moment of inertia of the planet."""

        return library.get_planet_inertia(self.c_body)

    #pylint: enable=arguments-differ

if __name__ == '__main__':
    #False positive.
    #pylint: disable=no-member
    planet = LockedPlanet((constants.M_jup / constants.M_sun).to(''),
                          (constants.R_jup / constants.R_sun).to(''))
    #pylint: enable=no-member
