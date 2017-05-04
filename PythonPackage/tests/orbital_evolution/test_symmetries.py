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
                        wind_saturation_frequency = 2.78,
                        diff_rot_coupling_timescale = 5.0e-3,
                        interpolator = interpolator)
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
                         initial_semimajor,
                         disk_dissipation_age,
                         secondary_angmom = None) :
    """Create a binary system to evolve from the given objects."""

    if isinstance(secondary, LockedPlanet) :
        secondary_config = dict(spin_angmom = numpy.array([0.0]),
                                inclination = None,
                                periapsis = None)
    else :
        secondary.select_interpolation_region(disk_dissipation_age)
        secondary_config = dict(spin_angmom = secondary_angmom,
                                inclination = numpy.array([0.0]),
                                periapsis = numpy.array([0.0]))
    secondary.configure(age = disk_dissipation_age,
                        companion_mass = primary.mass,
                        semimajor = initial_semimajor,
                        eccentricity = 0.0,
                        locked_surface = False,
                        zero_outer_inclination = True,
                        zero_outer_periapsis = True,
                        **secondary_config)
    if isinstance(secondary, EvolvingStar) :
        secondary.detect_stellar_wind_saturation()

    primary.select_interpolation_region(primary.core_formation_age())
    binary = Binary(primary = primary,
                    secondary = secondary,
                    initial_semimajor = initial_semimajor,
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

    primary.detect_stellar_wind_saturation()
    return binary

def output_evolution(evolution) :
    """Write the given evolution to stdout organized in columns."""

    quantities = list(
        filter(lambda q: q[0]!='_' and q != 'format', dir(evolution))
    )
    print(' '.join(['%25s' % q for q in quantities]))
    for i in range(len(evolution.age)) :
        print(' '.join(['%25s' % repr(getattr(evolution, q)[i])
                        for q in quantities]))


def plot_evolution(binary,
                   wsat,
                   style = dict(orb = '-r', core = '-b', env = '-g')) :
    """Calculate and plot the evolution of a properly constructed binary."""

    wsun = 2.0 * numpy.pi / 25.34

    binary.evolve(10.0, 1e-3, 1e-6, None)
    evolution = binary.get_evolution()
    
    worb = (2.0 * numpy.pi / binary.orbital_period(evolution.semimajor) 
            /
            wsun)
    wenv = (getattr(evolution, 'envelope_angmom',
                    getattr(evolution, 'primary_envelope_angmom', None))
            /
            binary.primary.envelope_inertia(evolution.age)) / wsun
    wcore = (getattr(evolution, 'core_angmom',
                     getattr(evolution, 'primary_core_angmom', None))
             /
             binary.primary.core_inertia(evolution.age)) / wsun

    pyplot.loglog(evolution.age, worb, style['orb'])
    pyplot.loglog(evolution.age, wenv, style['env'])
    pyplot.loglog(evolution.age, wcore, style['core'])
    pyplot.axhline(y = wsat / wsun, color = 'black')

    if isinstance(binary.secondary, EvolvingStar) :
        wenv = (evolution.secondary_envelope_angmom
                /
                binary.secondary.envelope_inertia(evolution.age)) / wsun
        wcore = (evolution.secondary_core_angmom
                 /
                 binary.secondary.core_inertia(evolution.age)) / wsun

        pyplot.loglog(evolution.age, wenv, style['sec_env'])
        pyplot.loglog(evolution.age, wcore, style['sec_core'])
    return evolution

def test_no_dissipation(interpolator, wind) :
    """Compare planet-non dissipative star to 2 non-dissipative stars."""

    tdisk = 5e-3

    star = create_star(interpolator = interpolator,
                       convective_phase_lag = 0.0,
                       wind = wind)
    planet = create_planet(1.0)
    binary = create_binary_system(star,
                                  planet,
                                  2.0 * numpy.pi / 3.0,
                                  10.0,
                                  tdisk)
    evolution = plot_evolution(binary, wsat = 2.78)
    planet.delete()
    star.delete()
    binary.delete()

    tdisk_index = evolution.age.searchsorted(tdisk)

    primary = create_star(interpolator = interpolator,
                          convective_phase_lag = 0.0,
                          wind = wind)
    secondary = create_star(interpolator = interpolator,
                            convective_phase_lag = 0.0,
                            wind = wind)
    binary = create_binary_system(
        primary,
        secondary,
        2.0 * numpy.pi / 3.0,
        10.0,
        tdisk,
        secondary_angmom = numpy.array([
            evolution.envelope_angmom[tdisk_index],
            evolution.core_angmom[tdisk_index]
        ])
    )
    evolution = plot_evolution(binary,
                               wsat = 2.78,
                               style = dict(orb = 'xr',
                                            core = 'xb',
                                            env = 'xg',
                                            sec_env = '+c',
                                            sec_core = '+m'))

    output_evolution(evolution)

    primary.delete()
    secondary.delete()
    binary.delete()
    pyplot.show()

def test_single_dissipator(interpolator, wind) :
    """Compare planet-dissipative star to 2 stars, one non-dissipative."""

    tdisk = 5e-3

    star = create_star(interpolator = interpolator,
                       convective_phase_lag = phase_lag(6.0),
                       wind = wind)
    planet = create_planet(1.0)
    binary = create_binary_system(star,
                                  planet,
                                  2.0 * numpy.pi / 3.0,
                                  10.0,
                                  tdisk)
    planet_evolution = plot_evolution(binary, wsat = 2.78)

    tdisk_index = planet_evolution.age.searchsorted(tdisk)
    initial_angmom = numpy.array([
        planet_evolution.envelope_angmom[tdisk_index],
        planet_evolution.core_angmom[tdisk_index]
    ])

    planet.delete()
    star.delete()
    binary.delete()

    print(80 * '=')
    print('Done with planet-star evolution')
    print(80 * '=')


    primary = create_star(interpolator = interpolator,
                          convective_phase_lag = phase_lag(6.0),
                          wind = wind)
    secondary = create_star(interpolator = interpolator,
                            convective_phase_lag = 0.0,
                            wind = wind)
    binary = create_binary_system(
        primary,
        secondary,
        2.0 * numpy.pi / 3.0,
        10.0,
        tdisk,
        secondary_angmom = initial_angmom
    )
    evolution = plot_evolution(binary,
                               wsat = 2.78,
                               style = dict(orb = 'xr',
                                            core = 'xb',
                                            env = 'xg',
                                            sec_env = ':c',
                                            sec_core = ':m'))

    output_evolution(evolution)

    primary.delete()
    secondary.delete()
    binary.delete()

    print(80 * '=')
    print('Done with dissipative primary star-star evolution')
    print(80 * '=')

    primary = create_star(interpolator = interpolator,
                          convective_phase_lag = 0.0,
                          wind = wind)
    secondary = create_star(interpolator = interpolator,
                            convective_phase_lag = phase_lag(6.0),
                            wind = wind)
    binary = create_binary_system(
        primary,
        secondary,
        2.0 * numpy.pi / 3.0,
        10.0,
        tdisk,
        secondary_angmom = initial_angmom
    )
    evolution = plot_evolution(binary,
                               wsat = 2.78,
                               style = dict(orb = '+r',
                                            core = 'oc',
                                            env = 'om',
                                            sec_env = '+g',
                                            sec_core = '+b'))

    output_evolution(evolution)

    primary.delete()
    secondary.delete()
    binary.delete()

    print(80 * '=')
    print('Done with dissipative secondary star-star evolution')
    print(80 * '=')

    pyplot.show()

if __name__ == '__main__' :
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '../../../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_single_dissipator(interpolator, True)
