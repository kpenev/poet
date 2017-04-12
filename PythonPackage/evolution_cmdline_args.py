"""Provide functions for adding parameters defining the evolution model."""

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from math import pi
from os.path import join as join_paths
from os.path import dirname

def add_disk_args(parser) :
    """Add command line options for the initial disk spin locking."""

    parser.add_argument(
        '--disk-lifetime', '--tdisk',
        type = float,
        default = 0.03,
        help = 'The scaling and powerlaw index of the age at which the '
        'protoplanetary disk evaporates in Gyrs as a function of the disk '
        'locking frequency. Default: %(default)s.'
    )
    parser.add_argument(
        '--disk-lock-frequency', '--wdisk',
        type = float,
        default = 3.5 * wsun,
        help = 'The disk locking frequency in rad/day. Default: %s.'
    )

def add_spindown_args(parser) :
    """Add command line options for spindown model parameters to parser."""

    wsun = 2.0 * pi / 25.34 #rad/day

    parser.add_argument(
        '--core-env-coupling-timescale',
        type = float,
        default = 1e-3,
        help = 'The timescale in Gyrs on which the differential rotation in '
        'the star decays. Default: %(default)s.'
    )
    parser.add_argument(
        '--wind-strength', '--Kw',
        type = float,
        default = 0.13,
        help = 'The normalization constant of the angular momentum loss rate'
        ' due to magnetized stellar wind. Default: %(default)s.'
    )
    parser.add_argument(
        '--wind-saturation-frequency', '--wsat',
        type = float,
        default = 2.78,
        help = 'The frequency at which the angular momentum loss due to wind'
        ' saturates. Default: %(default)s.'
    )

def add_orbital_evolution_args(parser) :
    """Add to parser all arguments controlling the orbital evolution."""

    parser.add_argument(
        '--interpolator-dir',
        default = join_paths(dirname(__file__),
                             '..',
                             'stellar_evolution_interpolators'),
        help = 'The directory contaning managed seralized interpolators to '
        'use. Default: \'%(default)s\'.'
    )
    parser.add_argument(
        '--lgQ',
        type = float,
        default = 6.0,
        help = 'log10(Q\'*). Default: %(default)s.'
    )
    parser.add_argument(
        '--max-evolution-step', '--max-step',
        type = float,
        default = 0.1,
        help = 'The maximum time step allowed when calculating individual '
        'orbital evolutions in Gyr. Default: %(default)s.'
    )
    parser.add_argument(
        '--evolution-precision', '--precision',
        type = float,
        default = 1e-8,
        help = 'The precision to require of the orbital evolution (both '
        'relative and absolute). Default: %(default)s.'
    )
    parser.add_argument(
        '--interpolator', '-I',
        default = 'default',
        help = 'The name of the interpolator which to use for the stellar '
        'evolution. Default: \'%(default)s\'.'
    )
    parser.add_argument(
        '--eccentricity-coef',
        default = 'eccentricity_expansion_coef.txt',
        help = 'The file containing the tabulated eccentriciy expansion '
        'coefficients for the tidal potential. Default: \'%(default)s\'.'
    )
    parser.add_argument(
        '--planet-formation-age', '--formation-age',
        type = float,
        default = 5e-3,
        help = 'The age at which the planet is assumed to form in Gyr. '
        'Default: %(default)s.'
    )

def add_and_parse_evolution_args(parser, disk_args = True) :
    """
    Add arguments for the binary evolution, some setup & return interpolator.

    Args:
        - parser:
            An instance of argparse.ArgumentParser to fill with arguments
            defining the evolution. Must alraedy contain all other arguments.

        - disk_args:
            Should command line arguments be added for configuring the
            initial circumstellar disk holding the primary's spin locked?

    Returns:
        - cmdline_args:
            An object containing the parsed command line arguments as
            members.

        - interpolator:
            A stellar evolution interpolator instance.
    """

    add_spindown_args(parser)
    if disk_args : add_disk_args(parser)
    add_orbital_evolution_args(parser)
    cmdline_args = parser.parse_args()
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        cmdline_args.eccentricity_coef.encode('ascii')
    )
    print('Interp dir: ' + repr(cmdline_args.interpolator_dir))
    interpolator = StellarEvolutionManager(
        cmdline_args.interpolator_dir
    ).get_interpolator_by_name(cmdline_args.interpolator)
    return cmdline_args, interpolator
