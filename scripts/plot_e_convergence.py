#!/usr/bin/env python3

"""
Plot eccentricity dependence of torque & power for different expansion orders.
"""

#matplotlib needs to be imported first.
#pylint: disable=wrong-import-order
from matplotlib import pyplot

from argparse import ArgumentParser
import re
#pylint: enable=wrong-import-order

import numpy

def parse_command_line():
    """Parse the command line."""

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        'e_dependence_fname',
        type=str,
        help='The name of the file containing the eccentricity dependence of '
        'the tidal torque and power for various expansion orders.'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='eccentricity_expansion_convergence.eps',
        help='The filename under which to save the generated plot. Default: '
        '%(default)s.'
    )
    return parser.parse_args()

def parse_column_name(column_name):
    """Return the quantity and eccentricity expansion order per column_name."""

    colname_rex = re.compile(
        r'^d(?P<quantity>[A-Za-z]+)_dt\[O\(e\^(?P<eorder>[0-9]+)\)\]$'
    )
    match = colname_rex.match(column_name)
    assert match
    return match.group('quantity'), int(match.group('eorder'))

def get_subplot_grid(num_quantities):
    """Return the number of rows/cols of sub-plots per number of quantities."""

    if num_quantities <= 2:
        return num_quantities, 1
    if num_quantities <= 6:
        return (num_quantities + 1) // 2, 2
    if num_quantities <= 12:
        return (num_quantities + 2) // 3, 3
    raise ValueError('Too many quantities to plot!')

def plot_quantity(tabulated_eccentricity, tabulated_values):
    """
    Plot a single quantity in the currently selected subplot.

    Args:
        tabulated_eccentricity (numpy.array):    The eccentricity values where
            the quantity values are known.

        tabulated_values (tuple(int, numpy.array)):    The quantity values
            calculated with different expansion orders in eccentricity. The
            first entry is the expansion order and the second should be numpy
            array matching in shape tabulated_eccentricity.

    Returns:
        None
    """

    tabulated_values.sort()
    reference_order, reference_values = tabulated_values[-1]
    for e_order, values in tabulated_values[:-1]:
        pyplot.semilogy(tabulated_eccentricity,
                        values / reference_values,
                        label=(r'$\mathcal{O}\left(e^{%d}\right)$' % e_order),
                        linewidth=1)

    pyplot.ylim(0.1, 10.0)
    pyplot.xlabel('eccentricity')
    pyplot.ylabel('expansion ratio')

def create_plots(cmdline_args):
    """Create the plots specified by the given command line arguments."""

    data = numpy.genfromtxt(cmdline_args.e_dependence_fname,
                            names=True,
                            deletechars='')
    tabulated_quantities = dict()
    for colname in data.dtype.names:
        if colname != 'e':
            quantity, e_order = parse_column_name(colname)
            if quantity not in tabulated_quantities:
                tabulated_quantities[quantity] = []
            tabulated_quantities[quantity].append((e_order, data[colname]))

    subplot_grid = get_subplot_grid(len(tabulated_quantities))

    for subplot_index, quantity in enumerate(tabulated_quantities):
        pyplot.subplot(*subplot_grid, subplot_index + 1)
        plot_quantity(data['e'], tabulated_quantities[quantity])
        pyplot.title(quantity)

    pyplot.figlegend()

    pyplot.savefig(cmdline_args.output)

if __name__ == '__main__':
    create_plots(parse_command_line())
