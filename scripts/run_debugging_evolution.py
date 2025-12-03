#/#!/usr/bin/env python3

"""Run evolutions with prescribed parameters and IC for debugging."""

import logging
import matplotlib
#matplotlib.use('TkAgg')

#pylint: disable=wrong-import-position
#pylint: disable=wrong-import-order
import sys

import matplotlib
from matplotlib import pyplot
from configargparse import ArgumentParser, DefaultsFormatter
import asteval
import numpy

sys.path.append('../PythonPackage')
sys.path.append('../scripts')

#pylint: disable=wrong-import-position
#pylint: disable=wrong-import-order
from orbital_evolution.command_line_util import\
    add_binary_config,\
    add_evolution_config,\
    run_evolution
#from stellar_evolution.library_interface import MESAInterpolator

from multiprocessing_util import setup_process

#pylint: enable=wrong-import-position
#pylint: enable=wrong-import-order

wnorm = 2.0 * numpy.pi / 25.34

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
    parser.add_argument(
        '--create-config',
        default=None,
        help='Filename to create a config file where all options are set per '
        'what is currently parsed.'
    )

    add_binary_config(parser)
    add_evolution_config(parser)

    parser.add_argument(
        '--plot',
        nargs='+',
        action='append',
        metavar=('FILENAME X_EXPR Y_EXPR [YEXPR ...]'),
        default=[],
        help='Add another plot to create. Each quantity can be a mathematical '
        'expression involving evolution quantities.'
    )

    parser.add_argument(
        '--plot-with-tangents',
        nargs=5,
        action='append',
        metavar=('FILENAME', 'X_EXPR', 'Y_EXPR', 'DYDX_EXPR', 'NUM_TANGENTS'),
        default=[],
        help='Add another plot that will also show tangent lines calculated '
        'assuming `DYDX_EXPR` evaluates to the slope at a given point. Tangent '
        'lines are drawn at the tabulated evolution points closest to '
        '`NUM_TANGENTS` evenly spaced values of `X_EXPR` covering the full '
        'range.'
    )
    parser.add_argument(
        '--logging-level',
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        default='info',
        help='The verbosity level of logging messages to use.'
    )

    result = parser.parse_args()
    logging.basicConfig(level=getattr(logging, result.logging_level.upper()))
    if result.create_config:
        print('Creating config file: ' + repr(result.create_config))
        parser.write_config_file(result,
                                 [result.create_config],
                                 exit_after=True)
    return result


def plot_tangents(plot_x, plot_y, plot_dydx, num_tangents):
    """Add tangent lines to the current plot."""

    def get_plot_data():
        """Return a list of the indices at which to plot tangent lines."""

        tangent_x_grid = numpy.linspace(numpy.min(plot_x),
                                        numpy.max(plot_x),
                                        num_tangents)
        plot_indices = numpy.array([0, plot_x.size - 1])
        for tangent_x in tangent_x_grid:
            abs_x_difference = numpy.abs(plot_x - tangent_x)
            plot_indices = numpy.concatenate(
                (
                    plot_indices,
                    numpy.argwhere(
                        numpy.logical_and(
                            abs_x_difference[:-2] > abs_x_difference[1:-1],
                            abs_x_difference[2:] > abs_x_difference[1:-1]
                        )
                    ).flatten() + 1
                )
            )
        plot_indices = numpy.unique(plot_indices)
        tangent_x = plot_x[plot_indices]
        plot_order = numpy.argsort(tangent_x)
        tangent_x = tangent_x[plot_order]
        plot_indices = plot_indices[plot_order]
        tangent_y = plot_y[plot_indices]
        tangent_slope = plot_dydx[plot_indices]

        return tangent_x, tangent_y, tangent_slope

    def add_tangent(previous_x, this_x, next_x, this_y, this_slope):
        """Add a single tangent line to the plot."""

        x_0 = 0.5 * (this_x + previous_x)
        x_1 = 0.5 * (this_x + next_x)
        y_0 = this_y - this_slope * (this_x - x_0)
        y_1 = this_y + this_slope * (x_1 - this_x)

        pyplot.plot([x_0, x_1], [y_0, y_1], '-')


    plot_limits = pyplot.xlim(), pyplot.ylim()
    tangent_x, tangent_y, tangent_slope = get_plot_data()

    previous_x = 2.0 * tangent_x[0] - tangent_x[1]
    this_x = tangent_x[0]
    for next_x, this_y, this_slope in zip(tangent_x[1:],
                                          tangent_y,
                                          tangent_slope):
        add_tangent(previous_x, this_x, next_x, this_y, this_slope)
        previous_x = this_x
        this_x = next_x

    add_tangent(tangent_x[-2],
                tangent_x[-1],
                2 * tangent_x[-1] - tangent_x[-2],
                tangent_y[-1],
                tangent_slope[-1])
    pyplot.xlim(plot_limits[0])
    pyplot.ylim(plot_limits[1])

def main(cmdline_args):
    """Avoid polluting the global namespace."""
    systemname='9881258'
    setup_process(
                    fname_datetime_format='%Y%m%d%H%M%S',
                    system=systemname,
                    std_out_err_fname='josh_output_10/{task}/{system}_{now}_{pid:d}.outerr',
                    logging_fname='josh_output_10/{task}/{system}_{now}_{pid:d}.log',
                    logging_verbosity='debug',
                    logging_message_format='%(levelname)s %(asctime)s %(name)s: %(message)s | %(pathname)s.%(funcName)s:%(lineno)d'
                  )

    evolution = run_evolution(cmdline_args, print_progress=True)
    #print(repr(evolution))
    for item in evolution.__dict__.items():
        print(item[0])
        print(item[1][-1])
    evaluator = asteval.Interpreter()
    evaluator.symtable.update(vars(evolution))
    for plot in cmdline_args.plot:
        plot_data = [evaluator(expression) for expression in plot[1:]]
        for plot_y, label in zip(plot_data[1:], plot[2:]):
            pyplot.plot(plot_data[0], plot_y, label=label)
        pyplot.legend()
        pyplot.savefig(plot[0])
        pyplot.clf()


    for plot in cmdline_args.plot_with_tangents:
        plot_tangents(*[evaluator(expression) for expression in plot[1:]])
        pyplot.savefig(plot[0])
        pyplot.clf()

if __name__ == '__main__':
    main(parse_configuration())
