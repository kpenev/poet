#!/usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')

from os import path
import sys
poet_root = path.dirname(
    path.dirname(
        path.abspath(__file__)
    )
)
sys.path.append(path.join(poet_root, 'PythonPackage'))
sys.path.append(path.join(poet_root, 'scripts'))
import numpy

from astropy import units, constants

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

if __name__ == '__main__':
    print('__file__ = ' + repr(__file__))
    interpolator = StellarEvolutionManager(
        path.join(poet_root, 'stellar_evolution_interpolators')
    ).get_interpolator_by_name(
        'default'
    )
    orbital_evolution_library.prepare_eccentricity_expansion(
        path.join(
            poet_root, 'eccentricity_expansion_coef_O400.sqlite'
        ).encode(
            'ascii'
        ),
        1e-6,
        True,
        False
    )

    star = EvolvingStar(mass = 1.0,
                        metallicity = 0.0,
                        wind_strength = 0.17,
                        wind_saturation_frequency = 2.45,
                        diff_rot_coupling_timescale = 5.0e-3,
                        interpolator = interpolator)
    star.select_interpolation_region(star.core_formation_age())
    star.set_dissipation(zone_index = 0,
                         tidal_frequency_breaks = None,
                         spin_frequency_breaks = None,
                         tidal_frequency_powers = numpy.array([0.0]),
                         spin_frequency_powers = numpy.array([0.0]),
                         reference_phase_lag = phase_lag(6.0))
    planet = LockedPlanet(
        mass = constants.M_jup.to_value('M_sun'),
        radius = constants.R_jup.to_value('R_sun')
    )

    disk_dissipation_age = 4e-3
    porb_initial = 3.0
    binary = Binary(primary = star,
                    secondary = planet,
                    initial_orbital_period = porb_initial,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.5,
                    disk_lock_frequency = 2.0 * numpy.pi / 7.0,
                    disk_dissipation_age = disk_dissipation_age,
                    secondary_formation_age = 4e-3)

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
                     eccentricity = 0.0,
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     locked_surface = False,
                     zero_outer_inclination = True,
                     zero_outer_periapsis = True)

    star.detect_stellar_wind_saturation()
    binary.evolve(10.0, 0.001, 1e-6, None)
    evolution = binary.get_evolution()
    print(repr(evolution))
    pyplot.plot(evolution.age, evolution.envelope_inclination)
    pyplot.plot(evolution.age, evolution.core_inclination)
    pyplot.show()
