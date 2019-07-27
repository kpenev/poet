#!/usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')

import sys
sys.path.append('../PythonPackage')
sys.path.append('../scripts')

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from orbital_evolution.initial_condition_solver import InitialConditionSolver
from basic_utils import Structure
import numpy
from astropy import units, constants

wsun = 2.0 * numpy.pi / 25.34

def create_planet(mass = (constants.M_jup / constants.M_sun).to(''),
                  phase_lag=0.0) :
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass = mass,
        radius = (constants.R_jup / constants.R_sun).to('')
    )
    if phase_lag:
        print('Setting planet dissipation')
        planet.set_dissipation(tidal_frequency_breaks = None,
                               spin_frequency_breaks = None,
                               tidal_frequency_powers = numpy.array([0.0]),
                               spin_frequency_powers = numpy.array([0.0]),
                               reference_phase_lag = phase_lag)
    return planet

def create_star(interpolator, convective_phase_lag) :
    """Create the star to use in the evolution."""

    star = EvolvingStar(mass = 1.0,
                        metallicity = 0.0,
                        wind_strength = 0.17,
                        wind_saturation_frequency = 2.45,
                        diff_rot_coupling_timescale = 5.0e-3,
                        interpolator = interpolator)
    if convective_phase_lag:
        print('Setting star dissipation')
        star.set_dissipation(zone_index = 0,
                             tidal_frequency_breaks = None,
                             spin_frequency_breaks = None,
                             tidal_frequency_powers = numpy.array([0.0]),
                             spin_frequency_powers = numpy.array([0.0]),
                             reference_phase_lag = convective_phase_lag)
    star.select_interpolation_region(star.core_formation_age())
    return star

def create_system(star,
                  planet,
                  disk_lock_frequency,
                  initial_eccentricity = 0.0) :
    """Create the system which to evolve from the given star and planet."""

    porb_initial = 3.5
    disk_dissipation_age = 4e-3
    binary = Binary(primary = star,
                    secondary = planet,
                    initial_orbital_period = porb_initial,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = disk_lock_frequency,
                    disk_dissipation_age = disk_dissipation_age,
                    secondary_formation_age = disk_dissipation_age)
    binary.configure(age = star.core_formation_age(),
                     semimajor = float('nan'),
                     eccentricity = float('nan'),
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     evolution_mode = 'LOCKED_SURFACE_SPIN')
    planet.configure(age = disk_dissipation_age,
                     companion_mass = star.mass,
                     semimajor = binary.semimajor(porb_initial),
                     eccentricity = initial_eccentricity,
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     locked_surface = False,
                     zero_outer_inclination = True,
                     zero_outer_periapsis = True)
    star.detect_stellar_wind_saturation()
    return binary

def test_evolution(interpolator,
                   convective_phase_lag = phase_lag(5.5),
                   planet_phase_lag = 0.0) :
    """Run a single orbital evolution calculation and plot the results."""

    for pdisk, color, wsat_enabled in [(1.4, 'r', '1')] :#,
#                                       (3.0, 'g', '2'),
#                                       (7.0, 'b', '3')] :
        star = create_star(interpolator = interpolator,
                           convective_phase_lag=convective_phase_lag)
        planet = create_planet(phase_lag=planet_phase_lag)
        binary = create_system(star, planet, 2.0 * numpy.pi / pdisk)

        binary.evolve(10.0, 0.001, 1e-6, None)
        print('====== FINAL STATE ======')
        print(binary.final_state().format())
        print('=========================')
        evolution_quantities = ['age',
                                'semimajor',
                                'eccentricity',
                                'envelope_angmom',
                                'core_angmom',
                                'wind_saturation',
                                'planet_angmom']
        evolution = binary.get_evolution(evolution_quantities)
        worb = (2.0 * numpy.pi / binary.orbital_period(evolution.semimajor)
                /
                wsun)
        wenv = (evolution.envelope_angmom
                /
                binary.primary.envelope_inertia(evolution.age)) / wsun
        wcore = (evolution.core_angmom
                 /
                 binary.primary.core_inertia(evolution.age)) / wsun

        planet_inertia = 0.3 * planet.mass * planet.radius**2

        print('Lplanet = ' + repr(evolution.planet_angmom))

        wplanet = (evolution.planet_angmom / planet_inertia) / wsun

        print('Wplanet = ' + repr(wplanet))

        numpy.savetxt(
            'Pdisk=%f.evol' % pdisk,
            numpy.dstack([evolution.age,
                          worb,
                          wenv,
                          wcore,
                          wplanet,
                          evolution.wind_saturation])[0],
            fmt = '%25s',
            header = ' '.join(
                ['%25s' % q
                 for q in ['t', 'worb', 'wenv', 'wcore', 'wplanet', 'wind_sat']]
            )
        )
#        pyplot.loglog(
#            evolution.age,
#            (
#                2.0 * numpy.pi / binary.orbital_period(evolution.semimajor)
#                /
#                wsun
#            ),
#            '-r'
#        )
        pyplot.loglog(evolution.age,  worb, '-' + color)
        pyplot.loglog(evolution.age[evolution.wind_saturation],
                      wenv[evolution.wind_saturation],
                      'o' + color,
                      markerfacecolor=color)
        pyplot.loglog(
            evolution.age[numpy.logical_not(evolution.wind_saturation)],
            wenv[numpy.logical_not(evolution.wind_saturation)],
            'o' + color,
            markerfacecolor='none'
        )
        pyplot.loglog(evolution.age,
                      wcore,
                      'x' + color)
        pyplot.loglog(evolution.age,
                      wplanet,
                      '.' + color)

#        wind_sat = numpy.zeros(evolution.age.shape)
#        wind_sat[evolution.wind_saturation] = wsat_enabled
#        pyplot.loglog(evolution.age, wind_sat, '.')

        pyplot.show()

        pyplot.plot(evolution.age, evolution.eccentricity)
        pyplot.show()

        star.delete()
        planet.delete()
        binary.delete()
    pyplot.loglog([4.6], [1.0], 'o')
    pyplot.loglog([4.6], [1.0], 'x')
#    pyplot.axhline(2.45 / wsun)
#    pyplot.ylim((0.1, 100))

def test_ic_solver(interpolator) :
    """Find initial condition to reproduce some current state and plot."""

    find_ic = InitialConditionSolver(disk_dissipation_age = 5e-3,
                                     evolution_max_time_step = 1e-2)
    target = Structure(age = 5.0,
                       Porb = 3.0,
                       Psurf = 10.0,
                       planet_formation_age = 5e-3)
    star = create_star(interpolator)
    planet = create_planet()
    initial_porb, initial_psurf = find_ic(target = target,
                                          star = star,
                                          planet = planet)
    print('IC: Porb0 = %s, P*0 = %s' % (repr(initial_porb),
                                        repr(initial_psurf)))

if __name__ == '__main__' :
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_evolution(interpolator,
                   planet_phase_lag=phase_lag(7.0),
                   convective_phase_lag=0.0*phase_lag(6.0))
#    test_ic_solver()
