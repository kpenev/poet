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
    parser.add_argument(
        '--eccentricity-expansion-order', '--e-order',
        type=int,
        default=10,
        help='The maximum power of eccentricity to include in the tidal '
        'expansion.'
    )

    add_binary_config(parser,
                      skip={'secondary_formation_age'})

    return parser.parse_args()

def main(config):
    """Avoid polluting global namespace."""

    config.final_age = config.disk_dissipation_age
    interpolator = set_up_library(config)
    binary = get_binary(config, interpolator)
    binary.change_e_order(config.eccentricity_expansion_order)

    _, all_axes = pyplot.subplots(2, 2)
    all_axes = all_axes.flatten()
    all_axes[0].set_title('Power')
    for axis, component in zip(all_axes[1:], 'XYZ'):
        axis.set_title(component + ' Torque')
    plot_x = numpy.arange(-config.eccentricity_expansion_order - 2,
                          config.eccentricity_expansion_order + 3)
    for spin_multiplier in range(-2, 3):
        tidal_rates = numpy.empty(shape=(4, plot_x.size), dtype=float)
        for rate_index, orbital_multiplier in enumerate(plot_x):
            tidal_rates[:, rate_index] = binary.envelope_tidal_torque_power(
                'primary',
                semimajor=binary.semimajor(config.initial_orbital_period),
                eccentricity=config.initial_eccentricity,
                age=config.disk_dissipation_age,
                spin=config.disk_lock_frequency,
                obliquity=config.initial_obliquity,
                tidal_term=(spin_multiplier, orbital_multiplier)
            )
        label='m=%d' % spin_multiplier
        for axis, plot_y in zip(all_axes, tidal_rates):
            selected = plot_y > 0
            color = axis.semilogy(
                plot_x[selected],
                plot_y[selected],
                'o',
                markerfacecolor='none',
                label=label
            )[0].get_color()
            label=None
            selected = numpy.logical_not(selected)
            axis.semilogy(plot_x[selected],
                          -plot_y[selected],
                          '*',
                          color=color,
                          markerfacecolor='none')

    total_rates = binary.envelope_tidal_torque_power(
        'primary',
        semimajor=binary.semimajor(config.initial_orbital_period),
        eccentricity=config.initial_eccentricity,
        age=config.disk_dissipation_age,
        spin=config.disk_lock_frequency,
        obliquity=config.initial_obliquity,
    )
    for axis, rate in zip(all_axes, total_rates):
        axis.axhline(y=rate, color='black')

    pyplot.figlegend()
    pyplot.show()

if __name__ == '__main__':
    main(parse_config())
