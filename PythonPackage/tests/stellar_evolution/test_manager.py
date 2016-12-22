#!/usr/bin/python3 -u

import os.path
import sys
sys.path.insert(0, os.path.abspath('../../'))

from stellar_evolution.manager import StellarEvolutionManager

if __name__ == '__main__' :
    manager = StellarEvolutionManager(
        '../../../stellar_evolution_interpolators'
    )

    interp = manager.get_interpolator(
        masses = [1.0],
        metallicities = [0.0],
        nodes = dict(RADIUS = 100,
                     ICONV = 100,
                     LUM = 100,
                     IRAD = 100,
                     MRAD = 100,
                     RRAD = 100),
        smoothing = dict(RADIUS = 0.1,
                         ICONV = 0.2,
                         LUM = 0.3,
                         IRAD = 0.4,
                         MRAD = 0.5,
                         RRAD = 0.6),
        new_interp_name = 'test_interp_nonan'
    )
    print(interp)

    interp = manager.get_interpolator(
        nodes = dict(RADIUS = 100,
                     ICONV = 100,
                     LUM = 100,
                     IRAD = 100,
                     MRAD = 100,
                     RRAD = 100),
        smoothing = dict(RADIUS = 0.1,
                         ICONV = 0.2,
                         LUM = float('nan'),
                         IRAD = 0.4,
                         MRAD = 0.5,
                         RRAD = 0.6),
        new_interp_name = 'test_interp_lumnan'
    )
    print(interp)
    exit(0)

    print(manager.get_by_name('test_interp_nonan'))

    print('Getting default interpolator.')
    print(manager.get_interpolator(new_interp_name = 'default'))
    print('Re-getting default interpolator')
    print(manager.get_interpolator())
    print('Done')
