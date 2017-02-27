#!/usr/bin/env python3

import sys
sys.path.append('../PythonPackage')

from orbital_evolution.evolve_interface import *
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from stellar_evolution.manager import StellarEvolutionManager

if __name__ == '__main__' :

    planet = LockedPlanet(1.0, 1.0)

    serialized_dir = '../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')
    star = EvolvingStar(1.0, 0.0, 0.15, 2.5, 5.0, interpolator)
    star.set_dissipation(0,
                         numpy.array([]),
                         numpy.array([]),
                         numpy.array([0.0]),
                         numpy.array([0.0]),
                         phase_lag(5))

    planet.configure(5e-3,
                     star.mass,
                     4.0,
                     0.0,
                     numpy.array([]),
                     numpy.array([]), 

    binary = Binary(primary = star,
                    secondary = planet,
                    initial_semimajor = 4.0,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = 2.0 * numpy.pi / 7.0,
                    disk_dissipation_age = 5e-3,
                    secondary_formation_age = 0.0)
    binary.evolve(8.0, 0.01, 1e-6, numpy.array([1.0, 2.0, 3.0]))
