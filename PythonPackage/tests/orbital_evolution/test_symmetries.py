#!/usr/bin/env python3

"""Test some evolution symmetries (i.e. exchangi primary<->secondary) etc."""

import matplotlib
matplotlib.use('TkAgg')

import os.path
import sys
sys.path.insert(0, os.path.abspath('../../'))

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
import numpy
from astropy import units, constants

def create_planet(mass = (constants.M_jup / constants.M_sun).to('')) :
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass = mass,
        radius = (constants.R_jup / constants.R_sun).to('')
    )
    return planet

def create_star(interpolator, convective_phase_lag, wind = True) :
    """Create the star to use in the evolution."""

    star = EvolvingStar(mass = 1.0,
                        metallicity = 0.0,
                        wind_strength = 0.17 if wind else 0.0,
                        wind_saturation_frequency = 2.45,
                        diff_rot_coupling_timescale = 5.0e-3,
                        interpolator = interpolator)
    star.select_interpolation_region(star.core_formation_age())
    star.set_dissipation(zone_index = 0,
                         tidal_frequency_breaks = None,
                         spin_frequency_breaks = None,
                         tidal_frequency_powers = numpy.array([0.0]),
                         spin_frequency_powers = numpy.array([0.0]),
                         reference_phase_lag = convective_phase_lag)
    star.set_dissipation(zone_index = 1,
                         tidal_frequency_breaks = None,
                         spin_frequency_breaks = None,
                         tidal_frequency_powers = numpy.array([0.0]),
                         spin_frequency_powers = numpy.array([0.0]),
                         reference_phase_lag = 0.0)
    return star

def create_binary_system(primary,
                         secondary,
                         disk_lock_frequency,
                         porb_initial = 2.5,
                         disk_dissipation_age = 4e-3) :
    """Create a binary system to evolve from the given objects."""

    binary = Binary(primary = primary,
                    secondary = secondary,
                    initial_orbital_period = porb_initial,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = disk_lock_frequency,
                    disk_dissipation_age = disk_dissipation_age,
                    secondary_formation_age = disk_dissipation_age)
    binary.configure(age = primary.core_formation_age(),
                     semimajor = float('nan'),
                     eccentricity = float('nan'),
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     evolution_mode = 'LOCKED_SURFACE_SPIN')
    if isinstance(secondary, LockedPlanet) :
        secondary_config = dict(spin_angmom = numpy.array([0.0]),
                                inclination = None,
                                periapsis = None)
    else :
        secondary.select_interpolation_region(disk_dissipation_age)
        secondary_config = dict(spin_angmom = numpy.array([1.0e-10, 1.0e-10]),
                                inclination = numpy.array([0.0]),
                                periapsis = numpy.array([0.0]))
    secondary.configure(age = disk_dissipation_age,
                        companion_mass = primary.mass,
                        semimajor = binary.semimajor(porb_initial),
                        eccentricity = 0.0,
                        locked_surface = False,
                        zero_outer_inclination = True,
                        zero_outer_periapsis = True,
                        **secondary_config)
    if isinstance(secondary, EvolvingStar) :
        secondary.detect_stellar_wind_saturation()
    primary.detect_stellar_wind_saturation()
    return binary

def plot_evolution(binary,
                   color = dict(orb = 'r', core = 'b', env = 'g')) :
    """Calculate and plot the evolution of a properly constructed binary."""

    wsun = 2.0 * numpy.pi / 25.34

    binary.evolve(10.0, 0.1, 1e-6, None)
    evolution = binary.get_evolution()

    worb = (2.0 * numpy.pi / binary.orbital_period(evolution.semimajor) 
            /
            wsun)
    wenv = (evolution.envelope_angmom
            /
            binary.primary.envelope_inertia(evolution.age)) / wsun
    wcore = (evolution.core_angmom
             /
             binary.primary.core_inertia(evolution.age)) / wsun

    pyplot.loglog(evolution.age, worb, '-' + color['orb'])
    pyplot.loglog(evolution.age, wenv, '.' + color['env'])
    pyplot.loglog(evolution.age, wcore, ':' + color['core'])

def test_primary_only_dissipative(interpolator) :
    """Compare planet-star evolution to 2 stars - one non dissipative."""

    print('Evolving star-planet system.')
    star = create_star(interpolator = interpolator,
                       convective_phase_lag = phase_lag(5.5))
    planet = create_planet(1.0)
    binary = create_binary_system(star, planet, 2.0 * numpy.pi / 3.0)
    plot_evolution(binary)
    planet.delete()
    star.delete()
    binary.delete()

    print('Evolving binary star system.')
    primary = create_star(interpolator = interpolator,
                          convective_phase_lag = phase_lag(5.5))
    secondary = create_star(interpolator = interpolator,
                            convective_phase_lag = 0.0)
    binary = create_binary_system(primary, secondary, 2.0 * numpy.pi / 3.0)
    plot_evolution(binary)
    planet.delete()
    star.delete()
    binary.delete()
    pyplot.show()

if __name__ == '__main__' :
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '../../../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_primary_only_dissipative(interpolator)
