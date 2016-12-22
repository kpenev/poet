#!/usr/bin/python3

import sys
sys.path.append('../PythonPackage')
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.library_interface import MESAInterpolator
import os.path
from glob import glob
import re

if __name__ == '__main__' :

    default_track_dir = '../MESA_tracks'
    serialized_dir = '../stellar_evolution_interpolators'

    manager = StellarEvolutionManager(serialized_dir)

    if not list(manager.get_suite_tracks()) :
        manager.register_track_collection(
            track_fnames = glob(os.path.join(default_track_dir, '*.csv')),
        )

    print(manager.get_interpolator(new_interp_name = 'default'))
