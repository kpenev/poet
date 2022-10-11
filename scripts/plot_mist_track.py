#!/usr/bin/env python3

"""Utilities for creating custom MIST tracks on the fly."""

from os import path
import logging
from multiprocessing import cpu_count

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from configargparse import ArgumentParser, DefaultsFormatter
from mesa_reader import MesaLogDir, MesaData
import numpy

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution import MISTTrackMaker, TemporaryMESAWorkDirectory
from stellar_evolution.library_interface import MESAInterpolator

#TODO: 1.5Msun needs 2000 nodes not 1000

def parse_command_line():
    """Return the configuration for what track to plot and how."""

    parser = ArgumentParser(
        description=__doc__,
        default_config_files=['plot_mist_track.cfg'],
        formatter_class=DefaultsFormatter,
        ignore_unknown_config_file_keys=False
    )
    parser.add_argument(
        '--stellar-evolution',
        nargs=2,
        metavar=('INTERP_DIR', 'INTERP_NAME'),
        default=(
            path.join(
                path.dirname(path.dirname(path.abspath(__file__))),
                'stellar_evolution_interpolators'
            ),
            'default'
        ),
        help='The directory containing serialized stellar evolutions and the '
        'name of the interpolator to use.'
    )
    parser.add_argument(
        '--mass',
        type=float,
        default=1.0,
        help='The mass of the star to plot the track for.'
    )
    parser.add_argument(
        '--feh',
        type=float,
        default=0.0,
        help='The [Fe/H] of the star to plot the track for.'
    )
    parser.add_argument(
        '--max-age',
        type=float,
        default=10.0,
        help='The maximum age to include in the plot in Gyr (also max age for '
        'the generated MIST track).'
    )
    parser.add_argument(
        '--mist-interp-min-age',
        type=float,
        default=1e-4,
        help='Ages earlier than this (in Gyrs) are not used when building the '
        'MIST interpolator to plot.'
    )
    parser.add_argument(
        '--profile-interval',
        type=int,
        default=None,
        help='How frequently to generate profiles. Plots are shown for each '
        'profile.'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='The filename to save the plots as (multi-page PDF).'
    )
    parser.add_argument(
        '--ncpu',
        type=int,
        default=cpu_count(),
        help='The maximum number of CPUs to use when running MESA.'
    )
    parser.add_argument(
        '--interp-config',
        nargs=5,
        metavar=('QUANTITY', 'NODES', 'SMOOTHING', 'LOGAGE', 'LOGQUANTITY'),
        action='append',
        default=[],
        help='Change the default configuration for interpolating one of the '
        'quantities. By default all quantities are interpolated with 1000 '
        'nodes, with smoothing 3 and using log(quantity) vs log(age). LOGAGE '
        'and LOGQUANTITY should be either ``T`` or ``F``.'
    )
    parser.add_argument(
        '--plot-logy',
        action='store_true',
        default=False,
        help='If passed, y axis of the plots is in log scale.'
    )
    result = parser.parse_args()
    for cfg in result.interp_config:
        assert cfg[3].upper() in 'TF'
        assert cfg[4].upper() in 'TF'

    return result


def plot_history(config, history, compare_interpolators, pdf):
    """Create plots involing the MESA history columns."""

    pyplot.plot(history.star_age,
                10.0**history.log_R,
                label=r'$R_\star$')
    pyplot.plot(history.star_age, history.rcore, label='$R_{core}$')

    compare_interpolated = [None for interpolator in compare_interpolators]
    interp_ages = [slice(False) for interpolator in compare_interpolators]
    for interp_i, interpolator in enumerate(compare_interpolators):
        if (
                hasattr(interpolator, 'in_range')
                and
                not interpolator.in_range(config.mass, config.feh)
        ):
            continue

        compare_interpolated[interp_i] = {
            q: interpolator(q, config.mass, config.feh)
            for q in interpolator.quantity_list
        }
        interp_ages[interp_i] = numpy.copy(
            history.star_age[
                numpy.logical_and(
                    (
                        history.star_age
                        >
                        compare_interpolated[interp_i]['RADIUS'].min_age
                    ),
                    (
                        history.star_age
                        <
                        compare_interpolated[interp_i]['RADIUS'].max_age
                    )
                )
            ]
        )


    for interp_i, interpolated in enumerate(compare_interpolated):
        if interpolated  is None:
            continue

        pyplot.plot(
            interp_ages[interp_i],
            interpolated['RRAD'](interp_ages[interp_i]),
            label=r'Interp %d: $R_{core}$' % interp_i
        )
        pyplot.plot(
            interp_ages[interp_i],
            interpolated['RADIUS'](interp_ages[interp_i]),
            label=r'Interp %d: $R_\star$' % interp_i
        )

    pyplot.axhline(y=0)
    pyplot.legend()
    finalize_plot(pdf, config)

    pyplot.plot(
        history.star_age,
        history.mcore / history.star_mass,
        label=r'$M_{core}/M_\star$'
    )
    for interp_i, interpolated in enumerate(compare_interpolated):
        if interpolated  is None:
            continue

        pyplot.plot(
            interp_ages[interp_i],
            (
                interpolated['MRAD'](interp_ages[interp_i])
                /
                config.mass
            ),
            label=r'Interp %d: $M_{core}/M_\star$' % interp_i
        )

    pyplot.axhline(y=0)
    pyplot.axhline(y=1)
    pyplot.legend()
    finalize_plot(pdf, config)

    pyplot.plot(history.star_age,
                history.core_inertia,
                label='$I_{core}$')
    pyplot.plot(history.star_age,
                history.env_inertia,
                label='$I_{env}$')
    pyplot.plot(history.star_age,
                history.core_inertia + history.env_inertia,
                label='$I_{tot}$')
    for interp_i, interpolated in enumerate(compare_interpolated):
        if interpolated  is None:
            continue

        pyplot.plot(
            interp_ages[interp_i],
            interpolated['IRAD'](interp_ages[interp_i]),
            label='Interp %d: $I_{core}$' % interp_i
        )
        pyplot.plot(
            interp_ages[interp_i],
            interpolated['ICONV'](interp_ages[interp_i]),
            label='Interp %d: $I_{env}$' % interp_i
        )

        pyplot.plot(
            interp_ages[interp_i],
            (
                interpolated['ICONV'](interp_ages[interp_i])
                +
                interpolated['IRAD'](interp_ages[interp_i])
            ),
            label='Interp %d: $I_{tot}$' % interp_i
        )


    pyplot.legend()
    finalize_plot(pdf, config)


def plot_profiles(config, mesa_data, pdf):
    """Create plots involing the MESA profiles."""

    history = mesa_data.history_data
    for model_number in mesa_data.model_numbers[-2:]:
        #False positive
        #pylint: disable=no-member
        profile_age = history.star_age[model_number - 1]
        #pylint: enable=no-member

        profile = mesa_data.profile_data(model_number=model_number)
    #    print('Mcore=%.3f, Rcore=%.3f, Ienv=%.3e, Icore=%.3e'
    #          %
    #          process_profile(profile))
        pyplot.plot(
            profile.radius / 10.0**history.log_R[model_number - 1],
            profile.mixing_type,
            '-g'
        )
        pyplot.axvline(
            (
                history.rcore[model_number - 1]
                /
                10.0**history.log_R[model_number - 1]
            ),
            color='green'
        )
        pyplot.plot(profile.mass / history.star_mass[model_number - 1],
                    profile.mixing_type,
                    '-r')
        pyplot.axvline(
            (
                history.mcore[model_number - 1]
                /
                history.star_mass[model_number - 1]
            ),
            color='red'
        )
        pyplot.title('Mixing($M$ and $R$) at $t=%.3e$ Gyr'
                     %
                     (profile_age/1e9))

        finalize_plot(pdf, config)


def finalize_plot(pdf, config):
    """Save or show the plot and get ready for next one."""

    pyplot.xscale('log')
    if config.plot_logy:
        pyplot.yscale('log')
    if config.output is None:
        pyplot.show()
        pyplot.clf()
    else:
        pdf.savefig()
        pyplot.close()


def main(config):
    """Avoid polluting global namespace."""

    manager = StellarEvolutionManager(config.stellar_evolution[0])
    compare_interpolators = [
        manager.get_interpolator_by_name(config.stellar_evolution[1])
    ]

    if config.output is None:
        pdf = None
    else:
        pdf = PdfPages(config.output)

    create_mist_track = MISTTrackMaker()
    interpolation_config = {
        q: dict(
            nodes=1000,
            smoothing=3,
            log_quantity=True,
            vs_log_age=True
        ) for q in MESAInterpolator.quantity_list
    }
    for (
            quantity,
            nodes,
            smoothing,
            vs_log_age,
            log_quantity
    ) in config.interp_config:
        assert quantity.upper() in interpolation_config
        interpolation_config[quantity.upper()] = dict(
            nodes=int(nodes),
            smoothing=float(smoothing),
            vs_log_age=(vs_log_age.upper() == 'T'),
            log_quantity=(log_quantity.upper() == 'T')
        )

    with TemporaryMESAWorkDirectory() as mesa_workdir:
        compare_interpolators.append(
            create_mist_track.create_interpolator(
                mass=config.mass,
                feh=config.feh,
                max_age=config.max_age,
                mesa_workdir=mesa_workdir,
                profile_interval=config.profile_interval,
                min_age=config.mist_interp_min_age,
                interpolation_config=interpolation_config,
                num_threads=config.ncpu
            )
        )
        if config.profile_interval:
            mesa_data = MesaLogDir(
                path.join(mesa_workdir, 'LOGS')
            )

            history = mesa_data.history_data
        else:
            history = MesaData(
                path.join(mesa_workdir, 'LOGS/history.data')
            )
        #False positive
        #pylint: disable=no-member
        history.star_age /= 1e9
        #pylint: enable=no-member

        plot_history(config, history, compare_interpolators, pdf)

        if config.profile_interval:
            plot_profiles(config, mesa_data, pdf)


    if config.output is not None:
        pdf.close()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main(parse_command_line())
