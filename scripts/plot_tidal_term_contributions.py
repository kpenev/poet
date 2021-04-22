#!/usr/bin/env python3

"""Visualize contributions of individual tidal terms to overall torque/power."""

import os.path

from matplotlib import pyplot
from configargparse import ArgumentParser, DefaultsFormatter

from orbital_evolution.command_line_util import\
    add_binary_config,\
    set_up_library,\
    get_binary

def parse_config():
    """Return the configuration for what/how to plot."""

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=DefaultsFormatter,
        args_for_setting_config_path=['-c', '--config-file'],
        ignore_unknown_config_file_keys=False
    )
    parser.add_argument(
        '--eccentricity_expansion_fname',
        default='eccentricity_expansion_coef.txt',
        help='The filename storing the eccentricity expansion coefficienst.'
    )
    parser.add_argument(
        '--stellar-evolution',
        nargs=2,
        metavar=('INTERP_DIR', 'INTERP_NAME'),
        default=(
            os.path.join(
                os.path.dirname(os.path.dirname(__file__)),
                'stellar_evolution_interpolators'
            ),
            'default'
        ),
        help='The directory containing serialized stellar evolutions and the '
        'name of the interpolator to use.'
    )

    add_binary_config(parser,
                      skip={'Tdisk', 'Wdisk', 'secondary_formation_age'})

    return parser.parse_args()

def main(config):
    """Avoid polluting global namespace."""

    interpolator = set_up_library(config)
    binary = get_binary(config, interpolator)

if __name__ == '__main__':
    main(parse_config())
