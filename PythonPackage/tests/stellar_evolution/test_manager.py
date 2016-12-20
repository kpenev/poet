#!/usr/bin/python3 -u

import os.path
import sys
sys.path.insert(0, os.path.abspath('../../'))

from stellar_evolution.manager import StellarEvolutionManager

if __name__ == '__main__' :
    print('Creating a manager')
    manager = StellarEvolutionManager(
        '../../../stellar_evolution_interpolators'
    )
    print('Getting default interpolator.')
    default_interp = manager.get_interpolator()
    print('Re-getting default interpolator')
    manager.get_interpolator()
    print('Done')
    exit(0)

    manager.get_interpolator(
        [os.path.join(StellarEvolutionManager.default_track_dir,
                      'M1.0_Z0.015.csv')],
        dict(RADIUS = 100,
             ICONV = 100,
             LUM = 100,
             IRAD = 100,
             MRAD = 100,
             RRAD = 100),
        dict(RADIUS = 0.1,
             ICONV = 0.2,
             LUM = 0.3,
             IRAD = 0.4,
             MRAD = 0.5,
             RRAD = 0.6),
        db_session
    )
