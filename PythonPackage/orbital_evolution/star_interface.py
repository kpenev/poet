#!/usr/bin/python3 -u

"""An interface to the POET star library."""

import sys
sys.path.append('..')

from stellar_evolution.library_interface import\
    library as stellar_evolution_library
from stellar_evolution.manager import StellarEvolutionManager
from ctypes import cdll, c_int, c_double, c_void_p, c_uint, c_bool, byref
import numpy

def initialize_library() :
    """Prepare the planet library for use."""

    library = cdll.LoadLibrary('libstar.so')

    library.create_star.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        stellar_evolution_library.create_interpolator.restype
    ]
    library.create_star.restype = c_void_p

    library.destroy_star.argtypes = [library.create_star.restype]
    library.destroy_star.restype = None

    library.set_dissipation.argtypes = [
        library.create_star.restype,
        c_uint,
        c_uint,
        c_uint,
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        c_double
    ]
    library.set_dissipation.restype = None

    library.detect_stellar_wind_saturation.argtypes = [
        library.create_star.restype
    ]
    library.detect_stellar_wind_saturation.restype = None

    library.modified_phase_lag.restype = c_double

    return library

library = initialize_library()

class EvolvingStar :
    """A class for stars following interpolated stellar evolution tracks."""

    deriv_list = ['NO', 'AGE', 'SPIN_FREQUENCY', 'ORBITAL_FREQUENCY']
    deriv_ids = {d: c_int.in_dll(library, d + '_DERIV').value
                 for d in deriv_list}

    def __init__(self,
                 mass,
                 metallicity,
                 wind_strength,
                 wind_saturation_frequency,
                 diff_rot_coupling_timescale,
                 interpolator) :
        """
        Create a star with the given properties and evolution.

        Args:
            - mass:
                The mass of the star in solar masses.
            - metallicity:
                The metallicity ([Fe/H]) of the star.
            - wind_strength:
                The efficiency of the wind carrying away angular momentum.
            - wind_saturation_frequency:
                The frequency at which the wind loss saturates in rad/day.
            - diff_rot_coupling_timescale:
                The timescale for differential rotation coupling.
            - interpolator:
                An instance of stellar_evolution.MESAInterpolator to base the
                stellar evolution on.

        Returns: None
        """

        self.mass = mass
        self.metallicity = metallicity
        self.wind_strength = wind_strength
        self.wind_saturation_frequency = wind_saturation_frequency
        self.diff_rot_coupling_timescale = diff_rot_coupling_timescale
        self.star = library.create_star(mass,
                                        metallicity,
                                        wind_strength,
                                        wind_saturation_frequency,
                                        diff_rot_coupling_timescale,
                                        interpolator.interpolator)

    def delete(self) :
        """Destroy the library star object created at construction."""

        library.destroy_star(self.star)

    def set_dissipation(self,
                        zone_index, 
                        tidal_frequency_breaks,
                        spin_frequency_breaks,
                        tidal_frequency_powers,
                        spin_frequency_powers,
                        reference_phase_lag) :
        """
        Set the dissipative properties of one of the zones of a star.

        Args:
            - zone_index:
                Which zone to set the dissiaption for (0 - envelope, 1 -
                core).
            - tidal_frequency_breaks:
                The locations of the breaks in tidal frequency in rad/day.
                Entries should be sorted.
            - spin_frequency_breaks:
                The locations of the breaks in spin frequency in rad/day.
                Entries should be sorted.
            - tidal_frequency_powers:
                The powerlaw indices for the tidal frequency dependence.
                Should be indexed in the same order as
                tidal_frequency_breaks, but must contain an additional
                starting entry for the powerlaw index before the first break.
            - spin_frequency_powers:
                The powerlaw indices for the spin frequency dependence.
                Should be indexed in the same order as spin_frequency_breaks,
                but must contain an additional starting entry for the
                powerlaw index before the first break.
            - reference_phase_lag:
                The phase lag at the first tidal and first spin frequency
                break. The rest are calculated by imposing continuity.

        Returns: None
        """

        library.set_dissipation(self.star,
                                zone_index,
                                tidal_frequency_breaks.size,
                                spin_frequency_breaks.size,
                                tidal_frequency_breaks,
                                spin_frequency_breaks,
                                tidal_frequency_powers,
                                spin_frequency_powers,
                                reference_phase_lag)

    def detect_stellar_wind_saturation(self) :
        """Tell a fully configured star to set its wind saturation state."""

        library.detect_stellar_wind_saturation(self.star)

    def modified_phase_lag(self,
                           zone_index,
                           orbital_frequency_multiplier,
                           spin_frequency_multiplier,
                           forcing_frequency,
                           deriv) :
        """
        Return the phase lag times the love number.

        The spin of the star must be set by calling the configure method.

        Args:
            - zone_index:
                The index of the zone whose lag to return.
            - orbital_frequency_multiplier:
                The multiplier of the orbital frequency in the expression for
                the forcing frequency.
            - spin_frequency_multiplier:
                The multiplier of the spin frequency in the expression for
                the forcing frequency.
            - forcing_frequency:
                The forcing frequency for which to return the phase lag.
            - deriv:
                One of the derivative IDs in self.deriv_ids identifying what
                derivative of the phase lag to return.

        Returns: The phase lag times the love number for the given
        parameters. If the forcing frequency is exactly 0.0, two values are
        returned, the first for infinitesimal positive and the second for
        infinitesimal negative forcing frequencies.
        """

        above_lock_value = c_double()
        below_lock_value = library.modified_phase_lag(
            self.star,
            c_uint(zone_index),
            c_int(orbital_frequency_multiplier),
            c_int(spin_frequency_multiplier),
            c_double(forcing_frequency),
            c_int(deriv),
            byref(above_lock_value)
        )
        if(forcing_frequency == 0) :
            return below_lock_value, above_lock_value.value
        else :
            return below_lock_value

if __name__ == '__main__' :
    serialized_dir = '../../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')
    star = EvolvingStar(1.0, 0.0, 0.15, 2.5, 5.0, interpolator)
    star.set_dissipation(0,
                         numpy.array([]),
                         numpy.array([]),
                         numpy.array([0.0]),
                         numpy.array([0.0]),
                         0.1)
    star.set_dissipation(1,
                         numpy.array([1.0]),
                         numpy.array([]),
                         numpy.array([0.0, 1.0]),
                         numpy.array([0.0]),
                         0.1)
    print('%25s %25s %25s' % ('w', 'Env(kDt)', 'Core(kDt)'))
    for w in numpy.linspace(0.1, 2.0, 20) :
        print('%25s %25s %25s'
              %
              (w, 
               repr(star.modified_phase_lag(0,
                                            1,
                                            1,
                                            w,
                                            star.deriv_ids['NO'])),
               repr(star.modified_phase_lag(1,
                                            1,
                                            1,
                                            w,
                                            star.deriv_ids['NO']))))


