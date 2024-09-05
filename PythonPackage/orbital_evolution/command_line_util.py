"""Add and use commandline/config file options for defining evolution to run."""

from math import pi
from os.path import dirname, join as join_paths
import logging

import numpy
from astropy import constants

from stellar_evolution.library_interface import MESAInterpolator
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

_logger = logging.getLogger(__name__)


def add_star_config(parser,
                    primary=True,
                    require_secondary=False,
                    dissipation=True):
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

        require_secondary(bool):    If true, removes the possibility of leaving
            secondary mass unspecified.

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
        default=(1.0 if (primary or require_secondary) else None),
        help=(
            'The mass of the ' + component_name + ' star.'
            +
            (
                ''
                if (primary or require_secondary) else
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
    if not primary:
        parser.add_argument(
            '--secondary-initial-angmom',
            type=float,
            nargs=2,
            default=None,
            help='The initial angular momentum of the secondary. Ignored if the'
            ' secondary is a planet. If the secondary is a star and this is not'
            ' specified the angular momentum is initialized like for the '
            'primary (disk formalism).'
        )
    if dissipation:
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
            prefix + '-dissipation-break',
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
        parser.add_argument(
            prefix + '-inertial-mode-enhancement',
            type=float,
            default=1.0,
            help='A factor by which the dissipation is larger in the inertial '
            'mode range relative to outside of it. This is applied on top of '
            'the spin and forcing frequency dependencies defined by the other '
            'arguments.'
        )
        parser.add_argument(
            prefix + '-inertial-mode-sharpness',
            type=float,
            default=10.0,
            help='A parameter controlling how suddenly the enhancement due to '
            'inertial modes gets turned on near the inertial mode range '
            'boundaries.'
        )
        parser.add_argument(
            prefix + '-age-dissipation-jump',
            nargs=2,
            type=float,
            action='append',
            metavar=('Age', 'JumpFactor'),
            help='Add a jump in the dissipation at a given age.'
        )


def add_binary_config(parser,
                      skip=(),
                      **star_config_kwargs):
    """
    Add command line/config file options to specify the binary to evolve.

    Args:
        parser:    The command line/cornfig file parser, or argument group, to
            add the options to.

        skip:    Collection of configuration options to exclude. Presumaly those
            will be denifen in some other way.

        require_secondary(bool):    If true, removes the possibility of leaving
            secondary mass unspecified.

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

    if 'dissipation' in skip:
        star_config_kwargs['dissipation'] = False

    add_star_config(parser, primary=True, **star_config_kwargs)
    add_star_config(parser, primary=False, **star_config_kwargs)

    if (
            'Porb' not in skip
            and
            'Porb0' not in skip
            and
            'initial_orbital_period' not in skip
    ):
        parser.add_argument(
            '--initial-orbital-period', '--Porb0', '--Porb',
            type=float,
            default=5.0,
            help='The initial orbital period (in days) in which the binary '
            'forms.'
        )

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
            '--initial-obliquity', '--Lambda0', '--Lambda',
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
        '--eccentricity-expansion-fname',
        default=join_paths(
            dirname(dirname(dirname(__file__))),
            'eccentricity_expansion_coef_O400.sqlite'
        ),
        help='The filename storing the eccentricity expansion coefficients.'
    )
    parser.add_argument(
        '--stellar-evolution',
        nargs=2,
        metavar=('INTERP_DIR', 'INTERP_NAME'),
        default=(
            join_paths(
                dirname(dirname(dirname(__file__))),
                'stellar_evolution_interpolators'
            ),
            'default'
        ),
        help='The directory containing serialized stellar evolutions and the '
        'name of the interpolator to use.'
    )
    parser.add_argument(
        '--interp-lower-limits',
        nargs='+',
        metavar='QUANTIY:LIMIT',
        default=[],
        help='Define lower limits to impose on interpolated stellar evolution '
        'quantities.'
    )
    parser.add_argument(
        '--max-evolution-runtime', '--timeout',
        type=float,
        default=0,
        help='The maximum number of seconds calculating the orbital evolution '
        'is allowed to take. Non-positive value results in no timeout. '
        'Partially cumputed evolutions that time out can still be querried.'
    )

def set_up_library(cmdline_args):
    """Define eccentricity expansion and return stellar evol interpolator."""

    print(f'Lower limits: {cmdline_args.interp_lower_limits!r}')
    for limit_str in cmdline_args.interp_lower_limits:
        quantity, limit = limit_str.split(':')
        MESAInterpolator.set_quantity_lower_limit(quantity, float(limit))

    orbital_evolution_library.prepare_eccentricity_expansion(
        cmdline_args.eccentricity_expansion_fname.encode('ascii'),
        1e-4,
        True,
        True
    )
    manager = StellarEvolutionManager(cmdline_args.stellar_evolution[0])
    return manager.get_interpolator_by_name(
        cmdline_args.stellar_evolution[1]
    )

#Pylint false positive for astropy constants
#pylint: disable=no-member
def create_planet(mass=(constants.M_jup / constants.M_sun).to(''),
                  radius=(constants.R_jup / constants.R_sun).to(''),
                  phase_lag=0.0):
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass=mass,
        radius=radius
    )
    if phase_lag:
        try:
            planet.set_dissipation(
                tidal_frequency_breaks=None,
                spin_frequency_breaks=None,
                age_breaks=None,
                tidal_frequency_powers=numpy.array([0.0]),
                spin_frequency_powers=numpy.array([0.0]),
                reference_phase_lags=numpy.array([float(phase_lag)])
            )
        except TypeError:
            planet.set_dissipation(**phase_lag)
    return planet
#pylint: enable=no-member

def create_star(interpolator,
                convective_phase_lag,
                *,
                mass=1.0,
                metallicity=0.0,
                wind_strength=0.17,
                wind_saturation_frequency=2.45,
                diff_rot_coupling_timescale=5.0e-3,
                interp_age=None):
    """Create the star to use in the evolution."""

    _logger.debug(
        'Creating %s Msun at t = %s star with dissipation %s',
        repr(mass),
        repr(interp_age),
        repr(convective_phase_lag)
    )
    star = EvolvingStar(mass=mass,
                        metallicity=metallicity,
                        wind_strength=wind_strength,
                        wind_saturation_frequency=wind_saturation_frequency,
                        diff_rot_coupling_timescale=diff_rot_coupling_timescale,
                        interpolator=interpolator)

    star.select_interpolation_region(star.core_formation_age()
                                     if interp_age is None else
                                     interp_age)

    if convective_phase_lag:
        try:
            star.set_dissipation(
                zone_index=0,
                tidal_frequency_breaks=None,
                spin_frequency_breaks=None,
                age_breaks=None,
                tidal_frequency_powers=numpy.array([0.0]),
                spin_frequency_powers=numpy.array([0.0]),
                reference_phase_lags=numpy.array([float(convective_phase_lag)])
            )
        except TypeError:
            star.set_dissipation(zone_index=0,
                                 **convective_phase_lag)
    return star

def create_system(primary,
                  secondary,
                  disk_lock_frequency,
                  *,
                  initial_eccentricity=0.0,
                  porb_initial=3.5,
                  disk_dissipation_age=4e-3,
                  initial_inclination=0.0,
                  secondary_formation_age=None,
                  secondary_initial_angmom=(0.1, 0.1)):
    """Combine the given primary and secondar in a system ready to evolve."""

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_orbital_period=porb_initial,
                    initial_eccentricity=initial_eccentricity,
                    initial_inclination=initial_inclination,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=(secondary_formation_age
                                             or
                                             disk_dissipation_age))
    binary.configure(
        age=(primary.core_formation_age() if isinstance(primary, EvolvingStar)
             else 0.5 * disk_dissipation_age),
        semimajor=float('nan'),
        eccentricity=float('nan'),
        spin_angmom=numpy.array([0.0]),
        inclination=None,
        periapsis=None,
        evolution_mode='LOCKED_SURFACE_SPIN'
    )

    if isinstance(secondary, EvolvingStar):
        initial_obliquity = numpy.array([0.0])
        initial_periapsis = numpy.array([0.0])
    else:
        initial_obliquity = None
        initial_periapsis = None
    secondary.configure(age=disk_dissipation_age,
                        companion_mass=primary.mass,
                        semimajor=binary.semimajor(porb_initial),
                        eccentricity=initial_eccentricity,
                        spin_angmom=(
                            numpy.array(secondary_initial_angmom)
                            if isinstance(secondary, EvolvingStar) else
                            numpy.array([0.0])
                        ),
                        inclination=initial_obliquity,
                        periapsis=initial_periapsis,
                        locked_surface=False,
                        zero_outer_inclination=True,
                        zero_outer_periapsis=True)

    if isinstance(primary, EvolvingStar):
        primary.detect_stellar_wind_saturation()
    if isinstance(secondary, EvolvingStar):
        secondary.detect_stellar_wind_saturation()

    return binary

def get_phase_lag_config(cmdline_args, primary=True):
    """Return a phase lag configuration to pass directly to create_star."""

    component_name = ('primary' if primary else 'secondary')
    age_breaks = getattr(cmdline_args, component_name + '_age_dissipation_jump')
    reference_dissipation = getattr(
        cmdline_args,
        component_name + '_reference_dissipation'
    )
    if reference_dissipation is None:
        return None

    reference_phase_lags = numpy.empty(
        (1 if age_breaks is None else len(age_breaks)),
        dtype=float
    )
    reference_phase_lags[0] = reference_dissipation[0]
    if age_breaks is not None:
        for i, (_, jump) in enumerate(age_breaks):
            reference_phase_lags[i + 1] = reference_phase_lags[i] * jump

    result = {
        'reference_phase_lags': reference_phase_lags,
        'spin_frequency_breaks': None,
        'spin_frequency_powers': numpy.array([0.0]),
        'age_breaks': (None if age_breaks is None
                       else numpy.array([age for age, _ in age_breaks]))
    }
    dissipation_breaks = (
        getattr(cmdline_args, component_name + '_dissipation_break')
        or
        []
    )

    if (
            reference_dissipation[2] == reference_dissipation[3] == 0
            and
            not dissipation_breaks
    ):
        result['tidal_frequency_powers'] = numpy.array([0.0])
        result['tidal_frequency_breaks'] = None
    else:
        result['tidal_frequency_breaks'] = numpy.empty(
            1 + len(dissipation_breaks),
            dtype=float
        )
        result['tidal_frequency_powers'] = numpy.empty(
            2 + len(dissipation_breaks),
            dtype=float
        )
        result['tidal_frequency_breaks'][0] = reference_dissipation[1]
        result['tidal_frequency_powers'][0] = reference_dissipation[2]
        result['tidal_frequency_powers'][1] = reference_dissipation[3]
        for index, (frequency, power) in enumerate(dissipation_breaks):
            result['tidal_frequency_breaks'][index + 1] = frequency
            result['tidal_frequency_powers'][index + 2] = power

    for param in ['inertial_mode_enhancement', 'inertial_mode_sharpness']:
        result[param] = getattr(cmdline_args, component_name + '_' + param)

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

    create_args = dict(
        interpolator=interpolator,
        convective_phase_lag=phase_lag_config,
        mass=mass,
        metallicity=cmdline_args.metallicity
    )
    for arg_name in ['wind_strength',
                     'wind_saturation_frequency',
                     'diff_rot_coupling_timescale']:
        create_args[arg_name] = getattr(
            cmdline_args,
            component_name + '_' + arg_name
        )
        if not primary and create_args[arg_name] is None:
            create_args[arg_name] = getattr(
                cmdline_args,
                'primary_' + arg_name
            )
    if not primary:
        create_args['interp_age'] = cmdline_args.disk_dissipation_age
    return create_star(**create_args)


def find_initial_secondary_angmom(cmdline_args,
                                  interpolator,
                                  **extra_evolve_args):
    """Run disk evolution for the secondary to find initial angmom."""

    _logger.debug('Finding initial secondary angular momentum using disks')
    binary = create_system(
        primary=get_component(cmdline_args, interpolator, False),
        secondary=get_component(cmdline_args, interpolator, True),
        disk_lock_frequency=cmdline_args.disk_lock_frequency,
        porb_initial=cmdline_args.initial_orbital_period,
        initial_eccentricity=cmdline_args.initial_eccentricity,
        initial_inclination=cmdline_args.initial_obliquity,
        disk_dissipation_age=cmdline_args.disk_dissipation_age,
        secondary_formation_age=cmdline_args.final_age,
        secondary_initial_angmom=(0.0, 0.0)
    )
    binary.evolve(
        cmdline_args.disk_dissipation_age,
        cmdline_args.max_time_step,
        cmdline_args.precision,
        create_c_code=False,
        eccentricity_expansion_fname=(
            cmdline_args.eccentricity_expansion_fname.encode('ascii')
        ),
        timeout=cmdline_args.max_evolution_runtime,
        **extra_evolve_args
    )
    final_state = binary.final_state()

    _logger.debug('Final state of secondary evolution: %s',
                  repr(final_state))

    binary.delete()
    return (
        final_state.primary_envelope_angmom,
        final_state.primary_core_angmom
    )


def get_binary(cmdline_args, interpolator):
    """Return the fully constructed binary to evolve."""

    return create_system(
        primary=get_component(cmdline_args, interpolator, True),
        secondary=get_component(cmdline_args, interpolator, False),
        disk_lock_frequency=cmdline_args.disk_lock_frequency,
        porb_initial=cmdline_args.initial_orbital_period,
        initial_eccentricity=cmdline_args.initial_eccentricity,
        initial_inclination=cmdline_args.initial_obliquity,
        disk_dissipation_age=cmdline_args.disk_dissipation_age,
        secondary_formation_age=(
            cmdline_args.secondary_formation_age
            if cmdline_args.secondary_mass else
            cmdline_args.final_age
        ),
        secondary_initial_angmom=cmdline_args.secondary_initial_angmom
    )

def run_evolution(cmdline_args,
                  interpolator=None,
                  required_ages_only=False,
                  final_state_only=False,
                  **extra_evolve_args):
    """Run the evolution specified on the command line."""

    if interpolator is None:
        interpolator = set_up_library(cmdline_args)

    if 'required_ages' not in extra_evolve_args:
        extra_evolve_args['required_ages'] = None

    if (
        cmdline_args.secondary_initial_angmom is None
        and
        cmdline_args.secondary_radius is None
    ):
        cmdline_args.secondary_initial_angmom = find_initial_secondary_angmom(
            cmdline_args,
            interpolator,
            **extra_evolve_args
        )
    binary = get_binary(cmdline_args, interpolator)
    _logger.debug('Evolving binary')
    binary.evolve(
        cmdline_args.final_age,
        cmdline_args.max_time_step,
        cmdline_args.precision,
        create_c_code=cmdline_args.create_c_code,
        eccentricity_expansion_fname=(
            cmdline_args.eccentricity_expansion_fname.encode('ascii')
        ),
        timeout=cmdline_args.max_evolution_runtime,
        **extra_evolve_args
    )
    if final_state_only:
        result = binary.final_state()
    else:
        result = binary.get_evolution(required_ages_only=required_ages_only)

    _logger.debug('Evolution result: %s', repr(result))

    for component_name in ['primary', 'secondary']:
        component = getattr(binary, component_name)
        if isinstance(component, EvolvingStar):
            for zone in ['core', 'envelope']:
                for deriv_order in range(1):
                    setattr(
                        result,
                        '_'.join(
                            [component_name, zone, 'inertia']
                            +
                            (['d%d' % deriv_order] if deriv_order else [])
                        ),
                        getattr(component, zone + '_inertia')(result.age,
                                                              deriv_order)
                    )
    result.orbital_period = binary.orbital_period(result.semimajor)
    binary.delete()
    return result
