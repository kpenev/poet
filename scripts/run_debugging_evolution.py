#!/usr/bin/env python3

"""Run evolutions with prescribed parameters and IC for debugging."""

import matplotlib
matplotlib.use('TkAgg')

#pylint: disable=wrong-import-position
#pylint: disable=wrong-import-order
import sys
sys.path.append('../PythonPackage')
sys.path.append('../scripts')

from matplotlib import pyplot
import numpy
from configargparse import ArgumentParser, DefaultsFormatter
import asteval

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.star_interface import EvolvingStar
from evolution_utils import create_star, create_planet, create_system

#pylint: enable=wrong-import-position
#pylint: enable=wrong-import-order

wnorm = 2.0 * numpy.pi / 25.34

def parse_configuration():
    """Return the configuration for what evolution to run and how."""

    def add_star_config(parser, primary=True):
        """
        Add to the parser arguments to configure a star.

        By default add arguments to configure the primary, use primary=False to
        configure the secondary. Secondary arguments other than mass all have
        None default values to allow falling back to the primary's values.
        """

        component_name = ('primary' if primary else 'secondary')
        prefix = '--' + component_name
        suffix = ('1' if primary else '2')
        extra_help = (''
                      if primary else
                      ' Defaults to primary star value if not  specified.')
        parser.add_argument(
            prefix + '-mass', '--m' + suffix,
            type=float,
            default=(1.0 if primary else None),
            help=(
                'The mass of the ' + component_name + ' star.'
                +
                (
                    ''
                    if primary else
                    ' If not specified, a single star evolution is calculated.'
                )
            )
        )
        parser.add_argument(
            prefix + '-radius', '--R' + suffix,
            type=float,
            default=None,
            help=(
                'Make the '
                +
                component_name
                +
                ' component a planet of the given radius. If not specified, '
                'but the mass is specified, the '
                +
                component_name
                +
                ' will be a star.'
            )
        )

        parser.add_argument(
            prefix + '-wind-strength', '--Kw' + suffix,
            type=float,
            default=(0.17 if primary else None),
            help=('The wind strength constant to assume for the '
                  +
                  component_name
                  +
                  ' star.'
                  +
                  extra_help)
        )
        parser.add_argument(
            prefix + '-wind-saturation-frequency', '--Wsat' + suffix,
            type=float,
            default=(2.45 if primary else None),
            help=('The wind saturation frequency to assume for the '
                  +
                  component_name
                  +
                  ' star in rad/day.'
                  +
                  extra_help)
        )
        parser.add_argument(
            prefix + '-diff-rot-coupling-timescale', '--Tcoup' + suffix,
            type=float,
            default=(5e-3 if primary else None),
            help=(
                'The timescale (in Gyrs) on which the rotation of the core and'
                ' the envelope of the '
                +
                component_name
                +
                ' star decay toward synchonism.'
                +
                extra_help
            )
        )
        parser.add_argument(
            prefix + '-reference-dissipation',
            nargs=4,
            metavar=('PhaseLag',
                     'TidalFrequency',
                     'PowerlawBefore',
                     'PowerlawAfter'),
            type=float,
            default=None,
            help=(
                'Define the reference point for the dissipation in the '
                +
                component_name
                +
                '. This initializes the phase lag dependence on frequnecy to a '
                'broken powerlaw with a single break. Further breaks can be '
                'introduced using '
                +
                prefix
                +
                '-dissipation-break.'
            )
        )
        parser.add_argument(
            prefix + '-dissipation_break',
            nargs=2,
            metavar=('TidalFrequency',
                     'PowerlawIndex'),
            action='append',
            help='Add another powerlaw break to the phase lag dependence on '
            'frequency. All breaks specified are sorted by their distance from '
            'the reference frequency (closest to furthest) and the powerlaw '
            'index at each break is defined to apply for the frequency range '
            'away from the reference.'
        )

    parser = ArgumentParser(
        description=__doc__,
        default_config_files=['system.params'],
        formatter_class=DefaultsFormatter,
        ignore_unknown_config_file_keys=False
    )
    parser.add_argument(
        '--metallicity', '--feh',
        type=float,
        default=0.0,
        help='The metallicity of the system to evolve.'
    )
    parser.add_argument(
        '--disk-dissipation-age', '--Tdisk',
        type=float,
        default=5e-3,
        help='The at which the disk dissipates and the binary evolution stars.'
    )
    parser.add_argument(
        '--disk-lock-frequency', '--Wdisk',
        type=float,
        default=(2.0 * numpy.pi / 7.0),
        help='The fixed spin frequency of the surface of the primary until the '
        'disk dissipates.'
    )

    add_star_config(parser, True)
    add_star_config(parser, False)

    parser.add_argument(
        '--initial-eccentricity', '--e0', '-e',
        type=float,
        default=0.0,
        help='The initial eccentricity at the time the secondary appears.'
    )
    parser.add_argument(
        '--initial_obliquity', '--Lambda0', '--Lambda',
        type=float,
        default=0.0,
        help='The initial obliquity of the orbit relative to the surface spin '
        'of the primary at the time the secondary appears.'
    )
    parser.add_argument(
        '--secondary-formation-age',
        type=float,
        default=None,
        help='The age at which the binary is formed. If left unspecified, the '
        'binary forms at the disk dissipation age.'
    )
    parser.add_argument(
        '--final-age',
        type=float,
        default=10.0,
        help='The age at which to stop calculating the evolution in Gyr.'
    )
    parser.add_argument(
        '--max-time-step',
        type=float,
        default=1e-3,
        help='The largest timestep in Gyrs the evolution is allowed to take.'
    )
    parser.add_argument(
        '--precision',
        type=float,
        default=1e-6,
        help='The precision to require of the solution.'
    )
    parser.add_argument(
        '--create-c-code',
        default='',
        help='If specified a compile-able C code is saved to a file with the '
        'given name that will calculate the specified evolution.'
    )
    parser.add_argument(
        '--eccentricity_expansion_fname',
        default='eccentricity_expansion_coef.txt',
        help='The filename storing the eccentricity expansion coefficienst.'
    )
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

def get_phase_lag_config(cmdline_args, primary=True):
    """Return a phase lag configuration to pass directly to create_star."""

    component_name = ('primary' if primary else 'secondary')
    reference_dissipation = getattr(cmdline_args,
                                    component_name + '_reference_dissipation')
    if reference_dissipation is None:
        return None

    result = dict(
        reference_phase_lag=reference_dissipation[0],
        spin_frequency_breaks=None,
        spin_frequency_powers=numpy.array([0.0]),
    )
    dissipation_breaks = getattr(cmdline_args,
                                 component_name + '_dissipation_break')

    if (
            reference_dissipation[2] == reference_dissipation[3] == 0
            and
            not dissipation_breaks
    ):
        result['tidal_frequency_powers'] = numpy.array([0.0])
        result['tidal_frequency_breaks'] = None
    else:
        raise NotImplementedError('Not implemented yet')

    return result

def get_component(cmdline_args, interpolator, primary=True):
    """Return one of the components of the system."""

    component_name = ('primary' if primary else 'secondary')

    if getattr(cmdline_args, component_name + '_mass') is None:
        assert not primary
        return create_planet()

    radius = getattr(cmdline_args, component_name + '_radius')
    phase_lag_config = get_phase_lag_config(cmdline_args, primary)
    mass = getattr(cmdline_args, component_name + '_mass')
    if radius is not None:
        return create_planet(mass=mass,
                             radius=radius,
                             phase_lag=phase_lag_config)

    return create_star(
        interpolator=interpolator,
        convective_phase_lag=phase_lag_config,
        mass=mass,
        metallicity=cmdline_args.metallicity,
        wind_strength=getattr(
            cmdline_args,
            component_name + '_wind_strength'
        ),
        wind_saturation_frequency=getattr(
            cmdline_args,
            component_name + '_wind_saturation_frequency'
        ),
        diff_rot_coupling_timescale=getattr(
            cmdline_args,
            component_name + '_diff_rot_coupling_timescale'
        )
    )

def get_binary(cmdline_args, interpolator):
    """Return the fully constructed binary to evolve."""

    return create_system(
        primary=get_component(cmdline_args, interpolator, True),
        secondary=get_component(cmdline_args, interpolator, False),
        disk_lock_frequency=cmdline_args.disk_lock_frequency,
        initial_eccentricity=cmdline_args.initial_eccentricity,
        initial_inclination=cmdline_args.initial_obliquity,
        secondary_formation_age=(
            cmdline_args.secondary_formation_age
            if cmdline_args.secondary_mass else
            cmdline_args.final_age
        )
    )

def run_evolution(cmdline_args):
    """Run the evolution specified on the command line."""

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        cmdline_args.eccentricity_expansion_fname.encode('ascii')
    )
    serialized_dir = '../stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')
    binary = get_binary(cmdline_args, interpolator)
    binary.evolve(
        cmdline_args.final_age,
        cmdline_args.max_time_step,
        cmdline_args.precision,
        None,
        create_c_code=cmdline_args.create_c_code,
        eccentricity_expansion_fname=(
            cmdline_args.eccentricity_expansion_fname.encode('ascii')
        )
    )
    evolution = binary.get_evolution()
    for component_name in ['primary', 'secondary']:
        component = getattr(binary, component_name)
        if isinstance(component, EvolvingStar):
            for zone in ['core', 'envelope']:
                setattr(
                    evolution,
                    '_'.join([component_name, zone, 'inertia']),
                    #False positive
                    #pylint: disable=no-member
                    getattr(component, zone + '_inertia')(evolution.age)
                    #pylint: enable=no-member
                )
    return evolution

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
