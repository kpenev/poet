#!/usr/bin/env python3

"""Generate plots of the modified phase lag."""

from os.path import dirname, abspath, join as join_paths
import logging

from matplotlib import pyplot
import numpy
from configargparse import ArgumentParser, DefaultsFormatter

from basic_utils import calc_semimajor
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.command_line_util import\
    add_star_config,\
    get_phase_lag_config

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
    parser.add_argument(
        '--stellar-evolution',
        nargs=2,
        metavar=('INTERP_DIR', 'INTERP_NAME'),
        default=(
            join_paths(
                dirname(dirname(abspath(__file__))),
                'stellar_evolution_interpolators'
            ),
            'default'
        ),
        help='The directory containing serialized stellar evolutions and the '
        'name of the interpolator to use.'
    )
    parser.add_argument(
        '--eccentricity-expansion-fname',
        default=join_paths(
            dirname(dirname(__file__)),
            'eccentricity_expansion_coef_O400.sqlite'
        ),
        help='The filename storing the eccentricity expansion coefficients.'
    )


    add_star_config(parser)
    parser.add_argument(
        '--metallicity', '--feh',
        type=float,
        default=0.0,
        help='The metallicity of the system to evolve.'
    )

    parser.add_argument(
        '--age',
        nargs='+',
        type=float,
        default=[1.0],
        help='The age or range of ages in Gyrs to evaluate the phase lag at. '
        'Only 1,2 or 3 arguments must be specified. If 1, the age of the star '
        'is fixed at the given value and plots are generated vs some other '
        'quantity. If 2, age is the x axis of the plot and phase lag is plotted'
        ' over the given age range. If 3 their meaning is min_age, max_age, '
        'nages, a series of plots is generated, each one with the age held '
        'fixed at one of the values given by linspace(min_age, max_age, nages).'
        ' The filename should then contain a {age} format substitution which '
        'will be replaced by the age the plot was evaluated at.'
    )
    parser.add_argument(
        '--spin-frequency',
        nargs='+',
        type=float,
        default=[1.0],
        help='Same as `--age` but sets the spin frequency.'
    )
    parser.add_argument(
        '--orbital-frequency',
        nargs='+',
        type=float,
        default=[-5.0, 5.0],
        help='Same as `--age` but sets the orbital frequency.'
    )
    parser.add_argument(
        '--spin-multiplier',
        type=int,
        nargs=2,
        default=[2, 2],
        help='The range of spin frequency multipliers to generate plots for '
        '(inclusive).'
    )
    parser.add_argument(
        '--orbital-multiplier',
        type=int,
        nargs=2,
        default=[2, 2],
        help='The range of orbital frequency multipliers to generate plots for '
        '(inclusive).'
    )
    parser.add_argument(
        '--plot-resolution',
        type=int,
        default=1000,
        help='The number of x values to evaluate the phase lag at for each '
        'plot.'
    )
    parser.add_argument(
        '--plot-fname',
        default='phase_lag_vs_worb_t{age!s}_wspin{spin_frequency!s}_'
        'm{spin_multiplier:d}_s{orbital_multiplier:d}.pdf',
        help='The filenames to save the plots as. Should contain format '
        'substitutions for all variables that will be evaluated at multiple '
        'values.'
    )
    result = parser.parse_args()
    plot_x = None
    for var_name in ['age', 'spin_frequency', 'orbital_frequency']:
        to_validate = getattr(result, var_name)
        if (
            len(to_validate) not in [1, 2, 3]
            or
            (
                len(to_validate) == 3
                and
                int(to_validate[2]) != to_validate[2]
            )
        ):
            parser.print_usage()
            print(
                'Invalid {!r} option {!r}. Should be 1 or 2 floats or 2 floats '
                'followed by an integer!'.format(
                    '--' + var_name.replace('_', '-'),
                    to_validate
                )
            )
        if len(to_validate) == 2:
            if plot_x is None:
                plot_x = var_name
            else:
                parser.print_help()
                print(
                    (
                        'Only one of `age`, `spin frequency`, or '
                        '`orbital frequency` can be the plot x axis. '
                        'Both {} and {} selected!'
                    ).format(
                        plot_x,
                        var_name
                    )
                )

    if plot_x is None:
        parser.print_help()
        print('No x axis selected!')

    return result


def get_plot_params(config):
    """
    Return x to plot vs and [(age, spin, orbital, m, s), ...] for each plot.
    """

    plot_params = [()]
    param_names = ['age',
                   'spin_frequency',
                   'orbital_frequency',
                   'spin_multiplier',
                   'orbital_multiplier']
    for var_name in param_names[:3]:
        var_config = getattr(config, var_name)
        if len(var_config) == 1:
            var_values = var_config
        elif len(var_config) == 2:
            var_values = [None]
            plot_x = numpy.linspace(*var_config, config.plot_resolution)
        else:
            var_values = numpy.linspace(var_config[0],
                                        var_config[1],
                                        int(var_config[2]))
        plot_params = [
            current + (new,)
            for current in plot_params
            for new in var_values
        ]
    for multiplier_config in [config.spin_multiplier,
                              config.orbital_multiplier]:
        plot_params = [
            current + (new, )
            for current in plot_params
            for new in range(multiplier_config[0], multiplier_config[1] + 1)
        ]
    plot_params = [
        dict(zip(param_names, param_values)) for param_values in plot_params
    ]
    return plot_x, plot_params


def get_plot_y(star, param_values, plot_x):
    """Return the y coordinates of the points for a single plot."""

    x_var_name = None
    for var_name in param_values.keys():
        if param_values[var_name] is None:
            x_var_name = var_name
            break
    assert x_var_name is not None

    assume_m2 = 1.0
    plot_y = numpy.empty(plot_x.shape)
    for y_dest, x_value in enumerate(plot_x):
        param_values[x_var_name] = x_value
        star.configure(
            age=param_values['age'],
            semimajor=calc_semimajor(
                star.mass,
                assume_m2,
                2.0 * numpy.pi / param_values['orbital_frequency']
            ),
            spin_angmom=numpy.array([
                (
                    star.envelope_inertia(param_values['age'])
                    *
                    param_values['spin_frequency']
                ),
                (
                    star.core_inertia(param_values['age'])
                    *
                    param_values['spin_frequency']
                )
            ]),
            companion_mass=assume_m2,
            eccentricity=0.0,
            inclination=numpy.array([0.0]),
            periapsis=numpy.array([0.0]),
            locked_surface=False,
            zero_outer_inclination=True,
            zero_outer_periapsis=True
        )
        plot_y[y_dest] = star.modified_phase_lag(
            zone_index=0,
            orbital_frequency_multiplier=param_values['orbital_multiplier'],
            spin_frequency_multiplier=param_values['spin_multiplier'],
            forcing_frequency=(
                (
                    param_values['spin_multiplier']
                    *
                    param_values['spin_frequency']
                )
                -
                (
                    param_values['orbital_multiplier']
                    *
                    param_values['orbital_frequency']
                )
            ),
            deriv=star.deriv_ids['NO']
        )

    param_values[x_var_name] = None
    return plot_y


def main(config):
    """Avoid polluting global namespace."""

    orbital_evolution_library.prepare_eccentricity_expansion(
        config.eccentricity_expansion_fname.encode('ascii'),
        0.1,
        True,
        True
    )

    star_interp = StellarEvolutionManager(
        config.stellar_evolution[0]
    ).get_interpolator_by_name(
        config.stellar_evolution[1]
    )
    star = EvolvingStar(
        metallicity=config.metallicity,
        interpolator=star_interp,
        **{
            var: getattr(config, 'primary_' + var)
            for var in ['mass',
                        'wind_strength',
                        'wind_saturation_frequency',
                        'diff_rot_coupling_timescale']
        }
    )
    star.set_dissipation(zone_index=0, **get_phase_lag_config(config))

    plot_x, plot_params = get_plot_params(config)
    for param_values in plot_params:
        plot_fname = config.plot_fname.format_map(param_values)
        plot_y = get_plot_y(star, param_values, plot_x)
        negative_y = (plot_y < 0)
        pyplot.plot(plot_x[numpy.logical_not(negative_y)],
                    plot_y[numpy.logical_not(negative_y)],
                    '.g',
                    markersize=0.2)
        pyplot.plot(plot_x[negative_y],
                    -plot_y[negative_y],
                    '.r',
                    markersize=0.2)
        pyplot.savefig(plot_fname)
        logging.info('Generated %s', plot_fname)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main(parse_configuration())
