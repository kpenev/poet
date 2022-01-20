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

from evolution_utils import creaty_star, create_planet, create_system

wnorm = 2.0 * numpy.pi / 25.34

def plot_evolution(binary, color):
    """Create and display plots of the calculated evolution."""

    print('====== FINAL STATE ======')
    print(binary.final_state().format())
    print('=========================')
    evolution_quantities = ['age',
                            'semimajor',
                            'eccentricity',
                            'wind_saturation']
    if isinstance(binary.secondary, EvolvingStar):
        evolution_quantities.extend([
            'primary_envelope_angmom',
            'primary_core_angmom',
            'secondary_envelope_angmom',
            'secondary_core_angmom'
        ])
    else:
        evolution_quantities.extend([
            'envelope_angmom',
            'core_angmom',
            'planet_angmom'
        ])

    evolution = binary.get_evolution(evolution_quantities)
    worb = (2.0 * numpy.pi / binary.orbital_period(evolution.semimajor)
            /
            wnorm)
    #False positive for astropy units & constants
    #pylint: disable=no-member
    orbital_energy = -(
        (
            (
                constants.G
                *
                binary.primary.mass * units.M_sun
                *
                binary.secondary.mass * units.M_sun
            )
            /
            (2.0 * evolution.semimajor * units.au)
        )
        /
        (
            units.M_sun * units.R_sun**2 / units.day**2
        )
    ).to('')
    #pylint: enable=no-member
    wenv = (
        getattr(
            evolution,
            'primary_envelope_angmom',
            getattr(evolution, 'envelope_angmom', None)
        )
        /
        binary.primary.envelope_inertia(evolution.age)
    ) / wnorm
    wcore = (
        getattr(
            evolution,
            'primary_core_angmom',
            getattr(evolution, 'core_angmom', None)
        )
        /
        binary.primary.core_inertia(evolution.age)
    ) / wnorm

    primary_spin_energy = (
        binary.primary.envelope_inertia(evolution.age)
        *
        (wenv * wnorm)**2
        +
        binary.primary.core_inertia(evolution.age)
        *
        (wcore * wnorm)**2
    ) / 2.0

    if isinstance(binary.secondary, EvolvingStar):
        wsecondary = (
            evolution.secondary_envelope_angmom
            /
            binary.secondary.envelope_inertia(evolution.age)
        ) / wnorm
        secondary_spin_energy = (
            evolution.secondary_envelope_angmom**2
            /
            binary.secondary.envelope_inertia(evolution.age)
            +
            evolution.secondary_core_angmom**2
            /
            binary.secondary.core_inertia(evolution.age)
        ) / 2.0
    else:
        secondary_spin_energy = numpy.zeros(evolution.age.shape)
        wsecondary = (
            evolution.planet_angmom
            /
            0.3 * binary.secondary.mass * binary.secondary.radius**2
        ) / wnorm

#    numpy.savetxt(
#        'Pdisk=%f.evol' % pdisk,
#        numpy.dstack([evolution.age,
#                      worb,
#                      wenv,
#                      wcore,
#                      wsecondary,
#                      evolution.wind_saturation])[0],
#        fmt = '%25s',
#        header = ' '.join(
#            ['%25s' % q
#             for q in ['t', 'worb', 'wenv', 'wcore', 'wsecondary', 'wind_sat']]
#        )
#    )
#        pyplot.loglog(
#            evolution.age,
#            (
#                2.0 * numpy.pi / binary.orbital_period(evolution.semimajor)
#                /
#                wnorm
#            ),
#            '-r'
#        )
    pyplot.loglog(evolution.age, worb, '-' + color)
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
                  wsecondary,
                  '.' + color)

#        wind_sat = numpy.zeros(evolution.age.shape)
#        wind_sat[evolution.wind_saturation] = wsat_enabled
#        pyplot.loglog(evolution.age, wind_sat, '.')

    pyplot.show()
    pyplot.cla()

    pyplot.semilogx(evolution.age, evolution.eccentricity)
    pyplot.show()
    pyplot.cla()

    pyplot.semilogx(evolution.age, orbital_energy, 'o' + color)
    pyplot.semilogx(evolution.age, primary_spin_energy, '--' + color)
    pyplot.semilogx(evolution.age, secondary_spin_energy, ':' + color)
    pyplot.semilogx(
        evolution.age,
        orbital_energy + primary_spin_energy + secondary_spin_energy,
        '-' + color
    )
    pyplot.show()
#    pyplot.axhline(2.45 / wnorm)
#    pyplot.ylim((0.1, 100))



def test_evolution(interpolator,
                   convective_phase_lag = phase_lag(5.5),
                   planet_phase_lag = 0.0,
                   create_c_code='',
                   eccentricity_expansion_fname=None):
    """Run a single orbital evolution calculation and plot the results."""

    for pdisk, color, wsat_enabled in [
#            (1.4, 'r', '1'),
#            (3.0, 'g', '2'),
            (7.0, 'b', '3')
    ]:
        star = create_star(interpolator=interpolator,
                           convective_phase_lag=convective_phase_lag)
        planet = create_planet(phase_lag=planet_phase_lag)
        binary = create_system(star, planet, 2.0 * numpy.pi / pdisk)

        binary.evolve(10.0,
                      0.001,
                      1e-6,
                      None,
                      create_c_code=create_c_code,
                      eccentricity_expansion_fname=eccentricity_expansion_fname)
        plot_evolution(binary, color)

        star.delete()
        planet.delete()
        binary.delete()


def test_binary_star_evolution(*,
                               interpolator,
                               convective_phase_lag=phase_lag(5.5),
                               create_c_code='',
                               eccentricity_expansion_fname=None):
    """Calculate and plot the evolution of a system of two stars."""

    pdisk = 7.0
    porb = 11.9
    star1_mass = 1.14
    star2_mass = 0.9
    feh = 0.21
    age = 6.3
    disk_lifetime = 5e-3
    eccentricity = 0.45

    primary = create_star(interpolator,
                          convective_phase_lag,
                          mass=star1_mass,
                          metallicity=feh)
    secondary = create_star(interpolator,
                            convective_phase_lag,
                            mass=star2_mass,
                            metallicity=feh,
                            interp_age=disk_lifetime)
    binary = create_system(primary,
                           secondary,
                           2.0 * numpy.pi / pdisk,
                           eccentricity,
                           porb_initial=porb,
                           disk_dissipation_age=disk_lifetime)
    binary.evolve(age,
                  0.001,
                  1e-6,
                  None,
                  create_c_code=create_c_code,
                  eccentricity_expansion_fname=eccentricity_expansion_fname)

    plot_evolution(binary, 'k')

    primary.delete()
    secondary.delete()
    binary.delete()

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

if __name__ == '__main__':
    eccentricity_expansion_fname = b"eccentricity_expansion_coef.txt"
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_expansion_fname
    )
    serialized_dir = '../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_binary_star_evolution(
        interpolator=interpolator,
#                               planet_phase_lag=phase_lag(5.0),
        convective_phase_lag=phase_lag(5.0),
        create_c_code='../poet_src/debug/test_evol.cpp',
        eccentricity_expansion_fname=eccentricity_expansion_fname
    )
#    test_ic_solver()
