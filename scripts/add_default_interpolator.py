#!/usr/bin/python3

import sys
sys.path.append('../PythonPackage')
from stellar_evolution.manager import StellarEvolutionManager
import os.path
from glob import glob
import re

if __name__ == '__main__' :

    default_track_dir = '../MESA_tracks'
    serialized_dir = '../stellar_evolution_interpolators'

    manager = StellarEvolutionManager(serialized_dir)

    if manager.get_suite_tracks() is None :
        manager.register_track_collection(
            glob(os.path.join(default_track_dir, '*.csv')),
            fname_rex = re.compile(
                'M(?P<MASS>[0-9.E+-]+)_Z(?P<Z>[0-9.E+-]+).csv'
            ),
            model_suite = 'MESA'
        )

#    manager.get_interpolator(new_interp_name = 'default')
