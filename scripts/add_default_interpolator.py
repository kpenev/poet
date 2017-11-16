#!/usr/bin/env python3

import sys
sys.path.append('../PythonPackage')
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.library_interface import MESAInterpolator
import os.path
from glob import glob
import re
from argparse import ArgumentParser

def parse_command_line() :
    """Return a structure with the command line arguments as members."""

    parser = ArgumentParser(
        description = 'Create an interpolator repository and add a single '
        'interpolator to it.'
    )
    parser.add_argument(
        '--repository-path', '--repository', '-r',
        type = str,
        default = '../stellar_evolution_interpolators',
        help = 'The path to the newly created repository. '
        'Default: %(default)s'
    )
    parser.add_argument(
        '--track-path', '--tracks', '-t',
        type = str,
        default = '../MESA_tracks',
        help = 'The directory containing the tracks to interpolate. '
        'Default: %(default)s'
    )
    parser.add_argument(
        '--trivial',
        action = 'store_true',
        default = False,
        help = 'If passed, instead of adding the default interpolator, adds '
        'an interpolator which does no smoothing for any quantity and thus '
        'is very fast to generate. Useful for debugging purposes.'
    )
    parser.add_argument(
        '--num-threads', '--threads',
        type = int,
        default = 4,
        help = 'The number of simultaneous threads to use when generating '
        'the interepolator. Default: %(default)d'
    )
    return parser.parse_args()

if __name__ == '__main__' :
    cmdline_args = parse_command_line()
    assert(os.path.isdir(cmdline_args.track_path))

    manager = StellarEvolutionManager(cmdline_args.repository_path)

    if not list(manager.get_suite_tracks()) :
        manager.register_track_collection(
            track_fnames = glob(
                os.path.join(cmdline_args.track_path, '*.csv')
            )
        )

    creation_args = dict(num_threads = cmdline_args.num_threads)
    if cmdline_args.trivial :
        nan = float('nan')
        creation_args['new_interp_name'] = 'trivial' 
        creation_args['nodes'] = {
            q: 0 for q in MESAInterpolator.quantity_list
        }
        creation_args['smoothing'] = {
            q: float('nan') for q in MESAInterpolator.quantity_list
        }
    else :
        creation_args['new_interp_name'] = 'default'

    print(manager.get_interpolator(**creation_args))
