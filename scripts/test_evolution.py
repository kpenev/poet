#!/usr/bin/env python3

import sys
sys.path.append('../PythonPackage')

from matplotlib import pyplot
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.evolve_interface import\
    DissipatingBody,\
    Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from stellar_evolution.manager import StellarEvolutionManager
import numpy

def create_planet(stellar_mass) :
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(mass = 1.0, radius = 1.0)
    planet.configure(age = 5e-3,
                     companion_mass = stellar_mass,
                     semimajor = 10.0,
                     eccentricity = 0.0,
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     locked_surface = False,
                     zero_outer_inclination = True,
                     zero_outer_periapsis = True)
    return planet

def create_star() :
    """Create the star to use in the evolution."""

    serialized_dir = '../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')
    star = EvolvingStar(mass = 1.0,
                        metallicity = 0.0,
                        wind_strength = 0.15,
                        wind_saturation_frequency = 2.5,
                        diff_rot_coupling_timescale = 5.0,
                        interpolator = interpolator)
    star.select_interpolation_region(star.core_formation_age())
    star.set_dissipation(zone_index = 0,
                         tidal_frequency_breaks = None,
                         spin_frequency_breaks = None,
                         tidal_frequency_powers = numpy.array([0.0]),
                         spin_frequency_powers = numpy.array([0.0]),
                         reference_phase_lag = phase_lag(6))
    star.set_dissipation(zone_index = 1,
                         tidal_frequency_breaks = None,
                         spin_frequency_breaks = None,
                         tidal_frequency_powers = numpy.array([0.0]),
                         spin_frequency_powers = numpy.array([0.0]),
                         reference_phase_lag = 0.0)
    return star

def create_system(star, planet) :
    """Create the system which to evolve from the given star and planet."""

    binary = Binary(primary = star,
                    secondary = planet,
                    initial_semimajor = 10.0,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = 2.0 * numpy.pi / 7.0,
                    disk_dissipation_age = 5e-3,
                    secondary_formation_age = 0.0)
    binary.configure(age = star.core_formation_age(),
                     semimajor = float('nan'),
                     eccentricity = float('nan'),
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     evolution_mode = 'LOCKED_SURFACE_SPIN')
    star.detect_stellar_wind_saturation()
    return binary

if __name__ == '__main__' :
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )

    star = create_star()
    planet = create_planet(star.mass)
    binary = create_system(star, planet)

    binary.evolve(10.0, 0.01, 1e-4, numpy.array([1.0, 2.0, 3.0]))
    print('====== FINAL STATE ======')
    print(binary.final_state().format())
    print('=========================')
    evolution = binary.get_evolution(['age', 'semimajor', 'envelope_angmom'])
    pyplot.plot(evolution.age,
                binary.orbital_period(evolution.semimajor),
                'xr')
    pyplot.plot(evolution.age,
                2.0 * numpy.pi *
                binary.primary.envelope_inertia(evolution.age)
                /
                evolution.envelope_angmom,
                'xg')
    pyplot.show()
