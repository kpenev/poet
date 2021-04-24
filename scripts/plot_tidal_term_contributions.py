#!/usr/bin/env python3

"""Visualize contributions of individual tidal terms to overall torque/power."""

import os.path

from matplotlib import pyplot
from configargparse import ArgumentParser, DefaultsFormatter
import numpy

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
        default='eccentricity_expansion_coef_O200.txt',
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
    parser.add_argument(
        '--eccentricity-expansion-order', '--e-order',
        type=int,
        default=10,
        help='The maximum power of eccentricity to include in the tidal '
        'expansion.'
    )

    add_binary_config(parser,
                      skip={'secondary_formation_age'})

    parser.add_argument(
        '--y-quantity',
        choices=['torque_power', 'e_rate'],
        default='torque_power',
        help='What quantity to plot vs the tidal term. For `torque_power`, '
        'four plots are generated showing the tidal power and the 3 components '
        'of the torque.'
    )

    return parser.parse_args()

def get_torque_power_rates(config):
    """
    Return an array of tidal torque and power for all tidal terms.

    Args:
    """

    config.final_age = config.disk_dissipation_age
    interpolator = set_up_library(config)
    binary = get_binary(config, interpolator)
    binary.change_e_order(config.eccentricity_expansion_order)

    orbital_multipliers = numpy.arange(-config.eccentricity_expansion_order - 2,
                                       config.eccentricity_expansion_order + 3)
    spin_multipliers = numpy.arange(-2, 3)
    tidal_rates = numpy.empty(
        (spin_multipliers.size, orbital_multipliers.size),
        dtype=[
            (quantity, float)
            for quantity in ['power', 'torque_x', 'torque_y', 'torque_z']
        ]
    )
    for spin_index, spin_mul in enumerate(spin_multipliers):
        for orbital_index, orbital_mul in enumerate(orbital_multipliers):
            tidal_rates[
                spin_index,
                orbital_index
            ] = binary.envelope_tidal_torque_power(
                'primary',
                semimajor=binary.semimajor(config.initial_orbital_period),
                eccentricity=config.initial_eccentricity,
                age=config.disk_dissipation_age,
                spin=config.disk_lock_frequency,
                obliquity=config.initial_obliquity,
                tidal_term=(spin_mul, orbital_mul)
            )

    return tidal_rates, spin_multipliers, orbital_multipliers

def plot_spin_index(orbital_multipliers,
                    plot_y,
                    plot_func,
                    spin_multiplier,
                    add_label=True):
    """Add the contributinos per orbital multiplier for a single spin mult."""

    markers = ['o*', 's+', '^v']

    selected = plot_y > 0
    color = plot_func(
        orbital_multipliers[selected],
        plot_y[selected],
        markers[abs(spin_multiplier)][0],
        markerfacecolor='none',
        label=('m=%+d' % spin_multiplier if add_label else None)
    )[0].get_color()

    selected = numpy.logical_not(selected)
    plot_func(orbital_multipliers[selected],
              -plot_y[selected],
              markers[abs(spin_multiplier)][1],
              color=color,
              markerfacecolor='none')


def plot_torque_power(config):
    """Create plots showing the tidal torque and power vs tidal term."""

    _, all_axes = pyplot.subplots(2, 2)
    all_axes = all_axes.flatten()
    all_axes[0].set_title('Power')
    for axis, component in zip(all_axes[1:], 'XYZ'):
        axis.set_title(component + ' Torque')

    tidal_rates, spin_multipliers, orbital_multipliers = get_torque_power_rates(
        config
    )

    for spin_index, spin_mul in enumerate(spin_multipliers):
        for axis, quantity in zip(all_axes, tidal_rates.dtype.names):
            plot_y = tidal_rates[spin_index][quantity]
            plot_spin_index(
                orbital_multipliers,
                plot_y,
                axis.semilogy,
                spin_mul,
                quantity == 'power'
            )

def plot_e_rate(config):
    """Create a plot of the circularization timescale vs tidal term."""

    (
        tidal_rates,
        spin_multipliers,
        orbital_multipliers
    ) = get_torque_power_rates(
        config
    )

    e_rates = (
        tidal_rates['power']
        -
        (
            2.0 * tidal_rates['torque_z']
            /
            numpy.sqrt(1.0 - config.initial_eccentricity**2)
        )
    )
    for spin_index, spin_mul in enumerate(spin_multipliers):
        plot_y = e_rates[spin_index]
        plot_spin_index(
            orbital_multipliers,
            plot_y,
            pyplot.semilogy,
            spin_mul
        )

def main(config):
    """Avoid polluting global namespace."""

    globals()['plot_' + config.y_quantity](config)

    pyplot.figlegend()
    pyplot.show()

if __name__ == '__main__':
    main(parse_config())
