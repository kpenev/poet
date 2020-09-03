#!/usr/bin/env python3

"""Run evolutions with prescribed parameters and IC for debugging."""

import matplotlib
matplotlib.use('TkAgg')

#pylint: disable=wrong-import-position
#pylint: disable=wrong-import-order
import sys
sys.path.append('../PythonPackage')
sys.path.append('../scripts')

from math import pi

from matplotlib import pyplot
from configargparse import ArgumentParser, DefaultsFormatter
import asteval

from orbital_evolution.command_line_util import\
    add_binary_config,\
    add_evolution_config,\
    run_evolution

#pylint: enable=wrong-import-position
#pylint: enable=wrong-import-order

wnorm = 2.0 * pi / 25.34

def parse_configuration():
    """Return the configuration for what evolution to run and how."""

    parser = ArgumentParser(
        description=__doc__,
        default_config_files=['system.params'],
        formatter_class=DefaultsFormatter,
        ignore_unknown_config_file_keys=False
    )
    parser.add_argument(
        '--config', '-c',
        is_config_file=True,
        help='Config file to use instead of default.'
    )

    add_binary_config(parser)
    add_evolution_config(parser)

    parser.add_argument(
        '--plot',
        nargs=2,
        action='append',
        metavar=('X_EXPR', 'Y_EXPR'),
        default=[],
        help='Add another plot to create. Each quantity can be a mathematical '
        'expression involving evolution quantities.'
    )

    return parser.parse_args()

def main(cmdline_args):
    """Avoid polluting the global namespace."""

    evolution = run_evolution(cmdline_args)
    print(evolution.format())
    evaluator = asteval.Interpreter()
    evaluator.symtable.update(vars(evolution))
    for plot in cmdline_args.plot:
        plot_x, plot_y = (evaluator(expression) for expression in plot)
        pyplot.plot(plot_x, plot_y)
        pyplot.show()
        pyplot.cla()

if __name__ == '__main__':
    main(parse_configuration())
