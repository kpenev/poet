"""Add command line/config file options for defining an evolution to run."""

from math import pi

def add_star_config(parser, primary=True):
    """
    Add to the parser arguments to configure a star.

    By default add arguments to configure the primary, use primary=False to
    configure the secondary. Secondary arguments other than mass all have
    None default values to allow falling back to the primary's values.

    Args:
        parser:    The command line/cornfig file parser, or argument group, to
            add the options to.

        primary(bool):    Whether the star being configured is the primary or
            the secondary in the system.

    Returns:
        None
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

def add_binary_config(parser, skip=()):
    """
    Add command line/config file options to specify the binary to evolve.

    Args:
        parser:    The command line/cornfig file parser, or argument group, to
            add the options to.

        skip:    Collection of configuration options to exclude. Presumaly those
            will be denifen in some other way.

    Returns:
        None
    """

    if 'feh' not in skip and 'metallicity' not in skip:
        parser.add_argument(
            '--metallicity', '--feh',
            type=float,
            default=0.0,
            help='The metallicity of the system to evolve.'
        )

    if 'Tdisk' not in skip and 'disk_dissipation_age' not in skip:
        parser.add_argument(
            '--disk-dissipation-age', '--Tdisk',
            type=float,
            default=5e-3,
            help='The at which the disk dissipates and the binary evolution '
            'stars.'
        )

    if 'Wdisk' not in skip and 'disk_lock_frequency' not in skip:
        parser.add_argument(
            '--disk-lock-frequency', '--Wdisk',
            type=float,
            default=(2.0 * pi / 7.0),
            help='The fixed spin frequency of the surface of the primary until '
            'the disk dissipates.'
        )

    add_star_config(parser, True)
    add_star_config(parser, False)

    if (
            'e' not in skip
            and
            'e0' not in skip
            and
            'initial_eccentricity' not in skip
    ):
        parser.add_argument(
            '--initial-eccentricity', '--e0', '-e',
            type=float,
            default=0.0,
            help='The initial eccentricity at the time the secondary appears.'
        )

    if (
            'Lambda' not in skip
            and
            'Lambda0' not in skip
            and
            'initial_obliquity' not in skip
    ):
        parser.add_argument(
            '--initial_obliquity', '--Lambda0', '--Lambda',
            type=float,
            default=0.0,
            help='The initial obliquity of the orbit relative to the surface '
            'spin of the primary at the time the secondary appears.'
        )

    if 'secondary_formation_age' not in skip:
        parser.add_argument(
            '--secondary-formation-age',
            type=float,
            default=None,
            help='The age at which the binary is formed. If left unspecified, '
            'the binary forms at the disk dissipation age.'
        )

def add_evolution_config(parser):
    """
    Add command line/config file options to specify how to run the evolution.

    Args:
        parser:    The command line/cornfig file parser, or argument group, to
            add the options to.

    Returns:
        None
    """

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
