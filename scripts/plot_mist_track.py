#!/usr/bin/env python3

"""Utilities for creating custom MIST tracks on the fly."""

from os import path, makedirs
import logging

from matplotlib import pyplot
from configargparse import ArgumentParser, DefaultsFormatter
from mesa_reader import MesaLogDir, MesaData
import numpy

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution import MISTTrackMaker, TemporaryMESAWorkDirectory
from stellar_evolution.MIST_realtime.build_model import _init_mesa_workdir

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
        '--profile-interval',
        type=int,
        default=None,
        help='How frequently to generate profiles. Plots are shown for each '
        'profile.'
    )
    return parser.parse_args()


def main(config):
    """Avoid polluting global namespace."""

    manager = StellarEvolutionManager(config.stellar_evolution[0])
    interpolator = manager.get_interpolator_by_name(
        config.stellar_evolution[1]
    )

    interpolated = {
        q: interpolator(q, 1.0, 0.0)
        for q in interpolator.quantity_list
    }

    create_mist_track = MISTTrackMaker()
    with TemporaryMESAWorkDirectory() as mesa_workdir:
        create_mist_track.run_mesa(mass=config.mass,
                                   feh=config.feh,
                                   max_age=config.max_age,
                                   mesa_workdir=mesa_workdir,
                                   profile_interval=config.profile_interval)
        if config.profile_interval:
            mesa_data = MesaLogDir(
                path.join(mesa_workdir, 'LOGS')
            )

            history = mesa_data.history_data
        else:
            history = MesaData(
                path.join(mesa_workdir, 'LOGS/history.data')
            )
        history.star_age /= 1e9
        interp_ages = numpy.copy(
            history.star_age[
                numpy.logical_and(
                    history.star_age > interpolated['RADIUS'].min_age,
                    history.star_age < interpolated['RADIUS'].max_age
                )
            ]
        )

        pyplot.semilogx(history.star_age,
                        10.0**history.log_R,
                        label='$R_\star$')
        pyplot.semilogx(history.star_age, history.rcore, label='$R_{core}$')
        pyplot.semilogx(
            interp_ages,
            interpolated['RRAD'](interp_ages),
            label='Interp: $R_{core}$'
        )
        pyplot.semilogx(
            interp_ages,
            interpolated['RADIUS'](interp_ages),
            label='Interp: $R_\star$'
        )
        pyplot.axhline(y=0)
        pyplot.legend()
        pyplot.show()
        pyplot.clf()

        pyplot.semilogx(
            history.star_age,
            history.mcore / history.star_mass,
            label='$M_{core}/M_\star$'
        )
        pyplot.semilogx(
            history.star_age,
            (
                interpolated['MRAD'](numpy.copy(history.star_age))
                /
                history.star_mass
            ),
            label='Interp: $M_{core}/M_\star$'
        )
        pyplot.axhline(y=0)
        pyplot.axhline(y=1)
        pyplot.legend()
        pyplot.show()
        pyplot.clf()

        pyplot.semilogx(history.star_age,
                        history.core_inertia,
                        label='$I_{core}$')
        pyplot.semilogx(history.star_age,
                        history.env_inertia,
                        label='$I_{env}$')
        pyplot.semilogx(history.star_age,
                        history.core_inertia + history.env_inertia,
                        label='$I_{tot}$')
        pyplot.semilogx(
            interp_ages,
            interpolated['IRAD'](interp_ages),
            label='Interp: $I_{core}$'
        )
        pyplot.semilogx(
            interp_ages,
            interpolated['ICONV'](interp_ages),
            label='Interp: $I_{env}$'
        )
        pyplot.semilogx(
            interp_ages,
            (
                interpolated['ICONV'](interp_ages)
                +
                interpolated['IRAD'](interp_ages)
            ),
            label='Interp: $I_{tot}$'
        )


        pyplot.legend()
        pyplot.show()
        pyplot.clf()

        if not config.profile_interval:
            return

        for model_number in mesa_data.model_numbers[-2:]:
            profile_age = mesa_data.history_data.star_age[model_number - 1]
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
            pyplot.show()
            pyplot.clf()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main(parse_command_line())

