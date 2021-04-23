#!/usr/bin/env python3

"""An interface to the POET star library."""

import sys
from ctypes import cdll, c_int, c_double, c_uint, byref
from ctypes.util import find_library

import numpy

sys.path.append('..')

#Need to add POEP packageto module search path before importing.
#pylint: disable=wrong-import-position
from stellar_evolution.library_interface import\
    library as stellar_evolution_library
from orbital_evolution.dissipating_body import\
    DissipatingBody
from orbital_evolution.evolve_interface import\
    library as orbital_evolution_library,\
    c_evolving_star_p
from orbital_evolution.c_interface_util import ndpointer_or_null
#pylint: enable=wrong-import-position

def initialize_library():
    """Prepare the planet library for use."""

    library_fname = find_library('star')
    if library_fname is None:
        raise OSError('Unable to find POET\'s star library.')
    result = cdll.LoadLibrary(library_fname)

    result.create_star.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        stellar_evolution_library.create_interpolator.restype
    ]
    result.create_star.restype = c_evolving_star_p

    result.destroy_star.argtypes = [result.create_star.restype]
    result.destroy_star.restype = None

    result.set_star_dissipation.argtypes = [
        result.create_star.restype,
        c_uint,
        c_uint,
        c_uint,
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_double,
        c_double,
        c_double
    ]
    result.set_star_dissipation.restype = None

    result.detect_stellar_wind_saturation.argtypes = [
        result.create_star.restype
    ]
    result.detect_stellar_wind_saturation.restype = None

    result.select_interpolation_region.argtypes = [
        result.create_star.restype,
        c_double
    ]
    result.select_interpolation_region.restype = None

    result.modified_phase_lag.restype = c_double

    result.core_formation_age.argtypes = [result.create_star.restype]
    result.core_formation_age.restype = c_double

    result.lifetime.argtypes = [result.create_star.restype]
    result.lifetime.restype = c_double

    result.luminosity.argtypes = [result.create_star.restype, c_double]
    result.luminosity.restype = c_double

    result.luminosity_array.argtypes = [
        result.create_star.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]
    result.luminosity_array.restype = None

    result.core_inertia.argtypes = [result.create_star.restype, c_double]
    result.core_inertia.restype = c_double

    result.core_inertia_deriv.argtypes = [result.create_star.restype,
                                          c_double,
                                          c_int]
    result.core_inertia_deriv.restype = c_double

    result.core_inertia_array.argtypes = [
        result.create_star.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]
    result.core_inertia_array.restype = None

    result.core_inertia_deriv_array.argtypes = [
        result.create_star.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_int,
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]
    result.core_inertia_array.restype = None

    result.envelope_inertia.argtypes = [result.create_star.restype,
                                        c_double]
    result.envelope_inertia.restype = c_double

    result.envelope_inertia_deriv.argtypes = [result.create_star.restype,
                                              c_double,
                                              c_int]
    result.envelope_inertia_deriv.restype = c_double

    result.envelope_inertia_array.argtypes = [
        result.create_star.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]
    result.envelope_inertia_array.restype = None

    result.envelope_inertia_deriv_array.argtypes = [
        result.create_star.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_int,
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]
    result.envelope_inertia_deriv_array.restype = None

    result.star_radius.argtypes = [result.create_star.restype, c_double]
    result.star_radius.restype = c_double

    result.star_radius_array.argtypes = [
        result.create_star.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]
    result.star_radius_array.restype = None

    return result

library = initialize_library()

class EvolvingStar(DissipatingBody):
    """A class for stars following interpolated stellar evolution tracks."""

    deriv_list = ['NO', 'AGE', 'SPIN_FREQUENCY', 'ORBITAL_FREQUENCY']
    deriv_ids = {d: c_int.in_dll(library, d + '_DERIV').value
                 for d in deriv_list}

    lib_configure_body = orbital_evolution_library.configure_star

    def _evaluate_stellar_property(self, property_name, age, deriv_order=None):
        """Evaluate a library function at a single age or array of ages."""

        if isinstance(age, numpy.ndarray):
            result = numpy.empty(dtype=c_double,
                                 shape=(age.size,),
                                 order='C')
            c_function = getattr(library, property_name + '_array')
            if deriv_order is None:
                c_function(self.c_body, age, age.size, result)
            else:
                c_function(self.c_body, age, deriv_order, age.size, result)
            return result

        if deriv_order is None:
            return getattr(library, property_name)(self.c_body, age)

        return getattr(library, property_name)(self.c_body,
                                               age,
                                               deriv_order)

    def __init__(self,
                 *,
                 mass,
                 metallicity,
                 wind_strength,
                 wind_saturation_frequency,
                 diff_rot_coupling_timescale,
                 interpolator):
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

        super().__init__()
        if (
                mass < interpolator.mass_range()[0]
                or
                mass > interpolator.mass_range()[1]
        ):
            raise ValueError(
                (
                    'Stellar mass: %s is outside the range supported by the '
                    'stelalr evolution interpolator: %s - %s'
                )
                %
                (
                    repr(mass),
                    repr(interpolator.mass_range()[0]),
                    repr(interpolator.mass_range()[1])
                )
            )
        if (
                metallicity < interpolator.feh_range()[0]
                or
                metallicity > interpolator.feh_range()[1]
        ):
            raise ValueError(
                (
                    'Stellar metallicity: %s is outside the range supported by '
                    'the stelalr evolution interpolator: %s - %s'
                )
                %
                (
                    repr(metallicity),
                    repr(interpolator.feh_range()[0]),
                    repr(interpolator.feh_range()[1])
                )
            )

        self.mass = mass
        self.metallicity = metallicity
        self.wind_strength = wind_strength
        self.wind_saturation_frequency = wind_saturation_frequency
        self.diff_rot_coupling_timescale = diff_rot_coupling_timescale
        self.interpolator = interpolator
        self.c_body = library.create_star(mass,
                                          metallicity,
                                          wind_strength,
                                          wind_saturation_frequency,
                                          diff_rot_coupling_timescale,
                                          interpolator.interpolator)

    def delete(self):
        """Destroy the library star object created at construction."""

        library.destroy_star(self.c_body)

    #The parent method simply saves the parameters, so it need not name them.
    #pylint: disable=arguments-differ
    def set_dissipation(self,
                        *,
                        zone_index,
                        tidal_frequency_breaks,
                        spin_frequency_breaks,
                        tidal_frequency_powers,
                        spin_frequency_powers,
                        reference_phase_lag,
                        inertial_mode_enhancement=1.0,
                        inertial_mode_sharpness=10.0):
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

            - inertial_mode_enhancement:
                A factor by which the dissipation is enhanced in the inertial
                mode range. Must be >= 1 (1 for no enhancement).

            - inertial_mode_sharpness:
                Parameter controlling how sharp the transition between enhanced
                and non-enhanced dissipation is.

        Returns: None
        """

        library.set_star_dissipation(self.c_body,
                                     zone_index,
                                     tidal_frequency_powers.size - 1,
                                     spin_frequency_powers.size - 1,
                                     tidal_frequency_breaks,
                                     spin_frequency_breaks,
                                     tidal_frequency_powers,
                                     spin_frequency_powers,
                                     reference_phase_lag,
                                     inertial_mode_enhancement,
                                     inertial_mode_sharpness)
        super().set_dissipation(
            zone_index=zone_index,
            tidal_frequency_breaks=tidal_frequency_breaks,
            tidal_frequency_powers=tidal_frequency_powers,
            spin_frequency_breaks=spin_frequency_breaks,
            spin_frequency_powers=spin_frequency_powers,
            reference_phase_lag=reference_phase_lag,
            inertial_mode_enhancement=inertial_mode_enhancement,
            inertial_mode_sharpness=inertial_mode_sharpness
        )
    #pylint: enable=arguments-differ

    def detect_stellar_wind_saturation(self):
        """Tell a fully configured star to set its wind saturation state."""

        library.detect_stellar_wind_saturation(self.c_body)

    def select_interpolation_region(self, age):
        """
        Prepare for interpolating stellar quantities around the given age.

        Args:
            - age:
                The age around which interpolation will be needed.

        Returns: None
        """

        library.select_interpolation_region(self.c_body, age)

    def modified_phase_lag(self,
                           *,
                           zone_index,
                           orbital_frequency_multiplier,
                           spin_frequency_multiplier,
                           forcing_frequency,
                           deriv):
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
            self.c_body,
            c_uint(zone_index),
            c_int(orbital_frequency_multiplier),
            c_int(spin_frequency_multiplier),
            c_double(forcing_frequency),
            c_int(deriv),
            byref(above_lock_value)
        )
        if forcing_frequency == 0:
            return below_lock_value, above_lock_value.value
        return below_lock_value

    def core_formation_age(self):
        """Return the age at which the core of the star forms in Gyrs."""

        return library.core_formation_age(self.c_body)

    def lifetime(self):
        """Return the maximum age at which the star can be queried."""

        return library.lifetime(self.c_body)

    def luminosity(self, age):
        """Return the luminosity of the star at the given age."""

        return self._evaluate_stellar_property('luminosity', age)

    def core_inertia(self, age, deriv_order=0):
        """
        Return the moment of inertia of the stellar core at the given age.
        """

        return self._evaluate_stellar_property('core_inertia_deriv',
                                               age,
                                               deriv_order)

    def envelope_inertia(self, age, deriv_order=0):
        """
        Return the moment of inertia of the stellar env. at the given age.
        """

        return self._evaluate_stellar_property('envelope_inertia_deriv',
                                               age,
                                               deriv_order)

    def radius(self, age):
        """Return the luminosity of the star at the given age."""

        return self._evaluate_stellar_property('star_radius', age)

if __name__ == '__main__':
    #Only needed for example under __main__
    #pylint: disable=ungrouped-imports
    from stellar_evolution.manager import StellarEvolutionManager
    #pylint: enable=ungrouped-imports

    def example():
        """An example of how to use this module."""

        serialized_dir = '../../stellar_evolution_interpolators'
        manager = StellarEvolutionManager(serialized_dir)
        interpolator = manager.get_interpolator_by_name('default')
        star1 = EvolvingStar(mass=1.0,
                             metallicity=0.0,
                             wind_strength=0.15,
                             wind_saturation_frequency=2.5,
                             diff_rot_coupling_timescale=5.0,
                             interpolator=interpolator)
        star2 = EvolvingStar(mass=1.0,
                             metallicity=0.0,
                             wind_strength=0.15,
                             wind_saturation_frequency=2.5,
                             diff_rot_coupling_timescale=5.0,
                             interpolator=interpolator)
        star1.set_dissipation(
            zone_index=0,
            tidal_frequency_breaks=numpy.array([]),
            spin_frequency_breaks=numpy.array([]),
            tidal_frequency_powers=numpy.array([0.0]),
            spin_frequency_powers=numpy.array([0.0]),
            reference_phase_lag=0.1
        )
        star1.set_dissipation(
            zone_index=1,
            tidal_frequency_breaks=numpy.array([1.0]),
            spin_frequency_breaks=numpy.array([]),
            tidal_frequency_powers=numpy.array([0.0, 1.0]),
            spin_frequency_powers=numpy.array([0.0]),
            reference_phase_lag=0.1
        )
        star2.set_dissipation(
            zone_index=0,
            tidal_frequency_breaks=numpy.array([0.6, 1.2]),
            spin_frequency_breaks=numpy.array([]),
            tidal_frequency_powers=numpy.array([0.0, 1.0, 0.0]),
            spin_frequency_powers=numpy.array([0.0]),
            reference_phase_lag=0.1
        )
        star2.set_dissipation(
            zone_index=1,
            tidal_frequency_breaks=numpy.array([0.5, 1.0, 1.5]),
            spin_frequency_breaks=numpy.array([]),
            tidal_frequency_powers=numpy.array([1.0, 2.0, 3.0, 4.0]),
            spin_frequency_powers=numpy.array([0.0]),
            reference_phase_lag=0.1
        )
        print('%25s %25s %25s %25s %25s'
              %
              ('w', 'Env1(kDt)', 'Core1(kDt)', 'Env2(kDt)', 'Core2(kDt)'))
        for wtide in numpy.linspace(0.1, 2.0, 20):
            print(
                '%25s %25s %25s %25s %25s'
                %
                (
                    wtide,
                    repr(star1.modified_phase_lag(
                        zone_index=0,
                        orbital_frequency_multiplier=1,
                        spin_frequency_multiplier=1,
                        forcing_frequency=wtide,
                        deriv=star1.deriv_ids['NO']
                    )),
                    repr(star1.modified_phase_lag(
                        zone_index=1,
                        orbital_frequency_multiplier=1,
                        spin_frequency_multiplier=1,
                        forcing_frequency=wtide,
                        deriv=star1.deriv_ids['NO']
                    )),
                    repr(star2.modified_phase_lag(
                        zone_index=0,
                        orbital_frequency_multiplier=1,
                        spin_frequency_multiplier=1,
                        forcing_frequency=wtide,
                        deriv=star2.deriv_ids['NO']
                    )),
                    repr(star2.modified_phase_lag(
                        zone_index=1,
                        orbital_frequency_multiplier=1,
                        spin_frequency_multiplier=1,
                        forcing_frequency=wtide,
                        deriv=star2.deriv_ids['NO']
                    ))
                )
            )

    example()
