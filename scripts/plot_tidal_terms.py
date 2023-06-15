#!/usr/bin/env python3

"""Create plots showing contribution of sepaarate tidal terms to evolution."""

from os import path

from matplotlib import pyplot
import numpy
from astropy import units as u
from configargparse import ArgumentParser, DefaultsFormatter

from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution import \
    SingleTermBody, \
    LockedPlanet, \
    phase_lag as calc_phase_lag

def parse_command_line():
    """Return the configuration for what evolution to run and how."""

    parser = ArgumentParser(
        description=__doc__,
        default_config_files=['plot_tidal_terms.cfg'],
        formatter_class=DefaultsFormatter,
        ignore_unknown_config_file_keys=False
    )
    parser.add_argument(
        '--config', '-c',
        is_config_file=True,
        help='Config file to use instead of default.'
    )
    parser.add_argument(
        '--eccentricity-expansion-fname',
        default=path.join(
            path.dirname(path.dirname(path.abspath(__file__))),
            'eccentricity_expansion_coef_O400.sqlite'
        ),
        help='The filename storing the eccentricity expansion coefficients.'
    )
    parser.add_argument(
        '--frequencies-vs-e',
        default=None,
        help='Specify a filename to save a plot of the leading tidal frequency,'
        ' the spin frequency, and the orbital period vs eccentricity. Empty '
        'string results in the plot being shown instead of saved. If not '
        'specified, the plot is not generated.'
    )
    parser.add_argument(
        '--top-terms-vs-e',
        default=None,
        help='Specify a filename to save a plot of the leading tidal term '
        '(orbit and spin multipliers) vs eccentricity. Empty string results in '
        'the plot being shown instead of saved. If not specified, the plot is '
        'not generated.'
    )
    parser.add_argument(
        '--spin-frequency', '--wspin',
        default='ps',
        help="The spin frequency to assume for the various plots. Use ``'ps'`` "
        '(default) to assume pseudo-synchronous rotation.'
    )
    parser.add_argument(
        '--semimajor-axis',
        default='const angmom',
        help='The semimajor axis to assume for frequencies plot. Use '
        "``'const_angmom'`` (default) to set it so angular momentum is the same"
        'for each point.'
    )
    parser.add_argument(
        '--single-term-rates-vs-spin',
        default=None,
        help='Enable plotting of the contribution of various tidal terms for '
        'circular orbit as a function of spin. Plot is saved with the specified'
        'filename, or shown if empty string.'
    )
    parser.add_argument(
        '--rate-spectra-eccentricities',
        type=float,
        nargs='+',
        default=[],
        help='Specify eccentricities for which to plot rate spectra. No figure '
        'is generated if empty.'
    )
    parser.add_argument(
        '--rate-spectra',
        default='',
        help='Filename to save the spectrum of tidal evolution waves. Plot is '
        'shown if empty string.'
    )
    parser.add_argument(
        '--min-lgQ',
        type=float,
        default=7.0,
        help='The log10(Q) value for tidal frequencies where dissiaption is '
        'saturated.'
    )
    parser.add_argument(
        '--Q-period-powerlaw',
        type=float,
        default=1.0,
        help='The powerlaw index of the dependence of Q on tidal period.'
    )
    parser.add_argument(
        '--Q-break-period',
        type=float,
        default=0.001,
        help='The break period to assume for Q vs Ptide in days.'
    )
    parser.add_argument(
        '--pseudo-synchronization',
        default=None,
        help='Enable plotting showing pseudo-synhcronization. Save the plots '
        'with the given filanem or show if empty string.'
    )
    return parser.parse_args()


#e is more readable than eccentricity in this case
#pylint: disable=invalid-name
def get_hut_pseudo_synchronous_spin(e):
    """
    Return the ratio of pseudo-synchronous to orbital angular velocity.

    Direct implementation of Hut (1981) Eq. 42.
    """

    return (
        1.0 + 7.5 * e**2 + 45.0 / 8.0 * e**4 + 5.0 / 16.0 * e**6
    ) / (
        (1.0 + 3.0 * e**2 + 3.0 / 8.0 * e**4) * (1.0 - e**2)**1.5
    )
#pylint: enable=invalid-name


def create_single_term_binary(*,
                              primary_mass=1.0 * u.M_sun,
                              secondary_mass=1.0 * u.M_sun,
                              primary_radius=1.0 * u.R_sun,
                              secondary_radius=1.0 * u.R_sun,
                              semimajor=10.0 * u.R_sun,
                              eccentricity=0.0,
                              expansion_order=100):
    """
    Create a binary with only a single dissipative tidal term in the primary.

    Args: The parameters of the system to create (all obliquities are zero).

    Returns:
        Binary:
            The binary containing the two objects configured in the given orbit.
    """

    primary = SingleTermBody(primary_mass.to_value(u.M_sun),
                             primary_radius.to_value(u.R_sun))
    secondary = SingleTermBody(secondary_mass.to_value(u.M_sun),
                               secondary_radius.to_value(u.R_sun))

    primary.set_dissipation(2, 2, 1.0)

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_semimajor=semimajor.to_value(u.R_sun),
                    initial_eccentricity=eccentricity,
                    disk_lock_frequency=1.0,
                    disk_dissipation_age=0.1)
    binary.set_expansion_order(expansion_order)

    worb = binary.orbital_frequency(semimajor.to_value(u.R_sun))

    binary_config = dict(
        age=1.0,
        semimajor=semimajor.to_value(u.R_sun),
        eccentricity=eccentricity,
        spin_angmom=(
            worb
            *
            get_hut_pseudo_synchronous_spin(eccentricity)
            *
            numpy.array([primary.inertia(), secondary.inertia()])
        ),
        inclination=numpy.array([0.0, 0.0]),
        periapsis=numpy.array([0.0]),
        evolution_mode='BINARY'
    )

    binary.configure(**binary_config)
    return binary, binary_config


def create_constant_time_lag_binary(*,
                                    time_lag,
                                    spin_orbit_ratio,
                                    primary_mass=1.0 * u.M_sun,
                                    secondary_mass=1.0 * u.M_sun,
                                    primary_radius=1.0 * u.R_sun,
                                    secondary_radius=1.0 * u.R_sun,
                                    semimajor=10.0 * u.R_sun,
                                    eccentricity=0.0):
    """
    Create a binary with only the primary having constant time lag dissipation.

    Args:
        time_lag:    The time lag to assume for the primary taides, including
            units (secondary is non-dissipative).

        The remaining arguments defined the parameters of the system to create
        (all obliquities are zero).

    Returns:
        Binary:
            The binary containing the two objects configured in the given orbit.

        dict:
            The configuration of the binary used to setup pseudo-synchronous
            rotation.
    """

    primary = LockedPlanet(primary_mass.to_value(u.M_sun),
                           primary_radius.to_value(u.R_sun))
    secondary = LockedPlanet(secondary_mass.to_value(u.M_sun),
                             secondary_radius.to_value(u.R_sun))

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_semimajor=semimajor.to_value(u.R_sun),
                    initial_eccentricity=eccentricity,
                    disk_lock_frequency=1.0,
                    disk_dissipation_age=0.1)

    primary.set_dissipation(
        tidal_frequency_breaks=numpy.array([2.0e2 * numpy.pi]),
        spin_frequency_breaks=None,
        tidal_frequency_powers=numpy.array([1.0, 0.0]),
        spin_frequency_powers=numpy.array([0.0]),
        reference_phase_lag=1e-2 * time_lag.to_value(u.day)
    )


    binary_config = dict(
        age=1.0,
        semimajor=semimajor.to_value(u.R_sun),
        eccentricity=eccentricity,
        spin_angmom=(
            spin_orbit_ratio
            *
            binary.orbital_frequency(semimajor.to_value(u.R_sun))
            *
            numpy.array([primary.inertia(), secondary.inertia()])
        ),
        inclination=numpy.array([0.0, 0.0]),
        periapsis=numpy.array([0.0]),
        evolution_mode='BINARY'
    )

    binary.configure(**binary_config)
    return binary, binary_config


def plot_pseudosynchronization(eccentricities=(0.1, 0.3, 0.5, 0.7),
                               scaled_spins=numpy.linspace(0.1, 10, 1000),
                               expansion_order=100,
                               log=False,
                               spin_scale='ps'):
    """Create a plot of torque vs spin rate for several eccentricities."""

    if log:
        pyplot.xscale('log')
        pyplot.yscale('log')
    for ecc in eccentricities:
        hut_pseudo_synchronous_factor = get_hut_pseudo_synchronous_spin(ecc)

        const_timelag_binary, binary_config = create_constant_time_lag_binary(
            spin_orbit_ratio=(1.0 if spin_scale == 'orbit'
                              else hut_pseudo_synchronous_factor),
            time_lag=1e-6 * u.day,
            eccentricity=ecc
        )
        const_timelag_binary.set_expansion_order(expansion_order)

        torque = numpy.empty(scaled_spins.shape)
        angmom_scale = binary_config['spin_angmom']
        print(
            'Const time lag rates at pseudo synchronization (e=%s): %s'
            %
            (
                ecc,
                repr(const_timelag_binary.calculate_rates(1.0)),
            )
        )

        for index, spin_factor in enumerate(scaled_spins):
            binary_config['spin_angmom'] = angmom_scale * spin_factor
            const_timelag_binary.configure(**binary_config)
            torque[index] = const_timelag_binary.calculate_rates(1.0)[-2]

        color = pyplot.plot(scaled_spins,
                            torque,
                            'o',
                            label='e=' + str(ecc))[0].get_color()
        if log:
            pyplot.plot(scaled_spins, -torque, 'x', color=color)
        if spin_scale == 'orbit':
            pyplot.axvline(x=hut_pseudo_synchronous_factor, color=color)

    if spin_scale == 'ps':
        pyplot.axvline(x=1, color='black')
    if not log:
        pyplot.axhline(y=0, color='black')
    pyplot.xlabel(r'$\Omega_{\star} / \Omega_{ps}$')
    pyplot.ylabel('Torque')

    pyplot.legend()


#No reasonable way found to simplify
#pylint: disable=too-many-locals
def get_rate_spectra(eccentricity,
                     *,
                     expansion_order=100,
                     max_phase_lag=1e-6,
                     tidal_powerlaw=0.0,
                     break_frequency=(2.0 * numpy.pi)):
    """Return rates and frequencies to plot."""

    semimajor = 10.0 * u.R_sun
    binary, binary_config = create_single_term_binary(
        eccentricity=eccentricity,
        semimajor=semimajor,
        expansion_order=expansion_order
    )

    worb = binary.orbital_frequency(semimajor.to_value(u.R_sun))
    wspin = binary_config['spin_angmom'][0] / binary.primary.inertia()
    num_terms = 10 * expansion_order + 5
    result = dict(
        frequencies=numpy.zeros(num_terms),
        phase_lags=numpy.zeros(num_terms),
        semimajor=numpy.zeros(num_terms),
        eccentricity=numpy.zeros(num_terms),
        spin=numpy.zeros(num_terms)
    )
    for worb_multiplier in range(-expansion_order, expansion_order + 1):
        for wspin_multiplier in range(-2, 3):
            if worb_multiplier == 0 and wspin_multiplier == 0:
                continue
            wtide = numpy.abs(worb_multiplier * worb - wspin_multiplier * wspin)
            phase_lag = max_phase_lag * min(
                1.0,
                (wtide / break_frequency)**tidal_powerlaw
            )
            binary.primary.set_dissipation(worb_multiplier,
                                           wspin_multiplier,
                                           phase_lag)
            rates = binary.calculate_rates(1.0)
            result_index = (5 * (worb_multiplier + expansion_order)
                            +
                            wspin_multiplier + 2)
            result['frequencies'][result_index] = wtide / worb
            result['phase_lags'][result_index] = phase_lag
            result['semimajor'][result_index] = rates[0]
            result['eccentricity'][result_index] = rates[1]
            result['spin'][result_index] = rates[-2] / binary.primary.inertia()

    return result
#pylint: enable=too-many-locals


def plot_rate_spectra(eccentricity, **rate_config):
    """Plot fourer components of the evolution rates."""

    rate_spectra = get_rate_spectra(eccentricity, **rate_config)
    sorted_frequencies = numpy.sort(rate_spectra['frequencies'])
    print('Min frequency difference: '
          +
          repr((sorted_frequencies[1:] - sorted_frequencies[:-1]).min()))
    color = None
    for plot_index, quantity in enumerate(['phase_lags',
                                           'semimajor',
                                           'eccentricity',
                                           'spin']):
        print('Max absolute %s rate: %s'
              %
              (repr(quantity), repr(numpy.abs(rate_spectra[quantity]).max())))
        pyplot.subplot(2, 2, plot_index + 1)
        pyplot.title(quantity)
        if plot_index == 0:
            plot_config = dict(label=('e = ' + str(eccentricity)))
        else:
            plot_config = dict(color=color)
        include = rate_spectra[quantity] > 1e-100
        if include.any():
            print('Positive range: (%s, %s)'
                  %
                  (
                      repr(rate_spectra[quantity][include].min()),
                      repr(rate_spectra[quantity][include].max())
                  ))
            color = pyplot.semilogy(rate_spectra['frequencies'][include],
                                    rate_spectra[quantity][include],
                                    'o',
                                    **plot_config)[0].get_color()
            plot_config = dict(color=color)
        print('Include %d positive values' % include.sum())


        include = rate_spectra[quantity] < -1e-100
        if include.any():
            print('Negative range: (%s, %s)'
                  %
                  (
                      repr(rate_spectra[quantity][include].min()),
                      repr(rate_spectra[quantity][include].max())
                  ))
            color = pyplot.semilogy(rate_spectra['frequencies'][include],
                                    -rate_spectra[quantity][include],
                                    'x',
                                    **plot_config)[0].get_color()
        print('Include %d negative values' % include.sum())


def plot_single_term_rates_vs_spin():
    """Plot evolution vs spin term by term."""

    semimajor = 10.0 * u.R_sun
    binary, binary_config = create_single_term_binary(
        eccentricity=0.0,
        semimajor=semimajor,
        expansion_order=10
    )

    worb = binary.orbital_frequency(semimajor.to_value(u.R_sun))
    print('worb = ' + repr(worb))
    scaled_spins = numpy.linspace(0.5, 1.5, 3)

    for worb_multiplier in [-2, 2]:
        for wspin_multiplier in [-2, 2]:
            binary.primary.set_dissipation(worb_multiplier,
                                           wspin_multiplier,
                                           1e-6)
            rates = dict(
                semimajor=numpy.empty(scaled_spins.shape),
                spin=numpy.empty(scaled_spins.shape)
            )
            for index, wspin_factor in enumerate(scaled_spins):
                binary_config['spin_angmom'] = (
                    wspin_factor
                    *
                    worb
                    *
                    numpy.array([binary.primary.inertia(),
                                 binary.secondary.inertia()])
                )
                binary.configure(**binary_config)
                spin_rates = binary.calculate_rates(1.0)
                rates['semimajor'][index] = spin_rates[0]
                rates['spin'][index] = spin_rates[-2] * binary.primary.inertia()
            pyplot.subplot(121)
            pyplot.plot(
                scaled_spins,
                rates['semimajor'],
                'o',
                label='%d worb - %d wspin' % (worb_multiplier, wspin_multiplier)
            )
            pyplot.subplot(122)
            pyplot.plot(scaled_spins, rates['spin'], 'o')
    pyplot.subplot(121)
    pyplot.axvline(x=1.0, color='black')
    pyplot.title('Semimajor rate')
    pyplot.legend()

    pyplot.subplot(122)
    pyplot.axvline(x=1.0, color='black')
    pyplot.title('Spin rate')
    pyplot.legend()


def get_leading_frequencies(binary, scaled_wspin, expansion_order):
    """Return the frequency of the most important tidal terms."""

    max_rates = dict(semimajor=-1.0, eccentricity=-1.0, spin=-1.0)
    leading_worb_multiplier = dict()
    leading_wspin_multiplier = dict()
    result = dict()
    for worb_multiplier in range(-expansion_order, expansion_order + 1):
        for wspin_multiplier in range(-2, 3):
            binary.primary.set_dissipation(worb_multiplier,
                                           wspin_multiplier,
                                           1e-6)
            abs_rates = numpy.absolute(binary.calculate_rates(1.0))
            for quantity, rate_ind in [('semimajor', 0),
                                       ('eccentricity', 1),
                                       ('spin', -2)]:
                if abs_rates[rate_ind] > max_rates[quantity]:
                    max_rates[quantity] = abs_rates[rate_ind]
                    result[quantity] = (float(worb_multiplier)
                                        -
                                        float(wspin_multiplier) * scaled_wspin)
                    leading_worb_multiplier[quantity] = worb_multiplier
                    leading_wspin_multiplier[quantity] = wspin_multiplier
    return result, leading_worb_multiplier, leading_wspin_multiplier


#This is straightforward to understand, no need to simplify
#pylint: disable=too-many-locals
def get_leading_frequencies_vs_e(eval_e,
                                 wspin='ps',
                                 semimajor='const angmom',
                                 expansion_order=100):
    """Return the the frequency of the most important term vs eccentricity."""

    result = dict(
        semimajor=numpy.empty(eval_e.shape),
        eccentricity=numpy.empty(eval_e.shape),
        spin=numpy.empty(eval_e.shape),
        orbital_freq=numpy.empty(eval_e.shape),
        spin_freq=numpy.empty(eval_e.shape),
        semimajor_worb_mult=numpy.empty(eval_e.shape, dtype=int),
        semimajor_wspin_mult=numpy.empty(eval_e.shape, dtype=int),
        eccentricity_worb_mult=numpy.empty(eval_e.shape, dtype=int),
        eccentricity_wspin_mult=numpy.empty(eval_e.shape, dtype=int),
        spin_worb_mult=numpy.empty(eval_e.shape, dtype=int),
        spin_wspin_mult=numpy.empty(eval_e.shape, dtype=int),
    )
    binary, binary_config = create_single_term_binary(
        eccentricity=eval_e[0],
        semimajor=(10.0 * u.R_sun if semimajor == 'const angmom'
                   else semimajor),
        expansion_order=expansion_order
    )

    if semimajor == 'const angmom':
        a_scale = 10.0 * u.R_sun * (1.0 + eval_e.max())
        worb_scale = binary.orbital_frequency(a_scale.to_value(u.R_sun))
    else:
        worb = binary.orbital_frequency(semimajor.to_value(u.R_sun))
        worb_scale = worb

    if wspin != 'ps':
        binary_config['spin_angmom'] = (
            wspin
            *
            numpy.array([binary.primary.inertia(), binary.secondary.inertia()])
        )


    for index, eccentricity in enumerate(eval_e):
        binary_config['eccentricity'] = eccentricity
        if semimajor == 'const angmom':
            binary_config['semimajor'] = (a_scale.to_value(u.R_sun)
                                          /
                                          (1.0 - eccentricity**2))
            worb = binary.orbital_frequency(binary_config['semimajor'])

        if wspin == 'ps':
            scaled_wspin = get_hut_pseudo_synchronous_spin(eccentricity)
            binary_config['spin_angmom'] = (
                worb
                *
                scaled_wspin
                *
                numpy.array([binary.primary.inertia(),
                             binary.secondary.inertia()])
            )

        binary.configure(**binary_config)
        (
            leading_wtide,
            leading_worb_multiplier,
            leading_wspin_multiplier
        ) = get_leading_frequencies(binary,
                                    scaled_wspin,
                                    expansion_order)
        for quantity, wtide in leading_wtide.items():
            result[quantity][index] = wtide * worb / worb_scale
            result[quantity + '_worb_mult'][index] = (
                leading_worb_multiplier[quantity]
            )
            result[quantity + '_wspin_mult'][index] = (
                leading_wspin_multiplier[quantity]
            )


        result['orbital_freq'][index] = worb / worb_scale
        result['spin_freq'][index] = scaled_wspin * worb / worb_scale

    return result
#pylint: enable=too-many-locals


def plot_leading_term_vs_e(eval_e,
                           frequencies_fname,
                           multipliers_fname,
                           *,
                           wspin='ps',
                           semimajor='const angmom',
                           expansion_order=100):
    """Make a plot of the frequency of the most important term vs eccentr."""

    leading_frequencies = get_leading_frequencies_vs_e(eval_e,
                                                       wspin,
                                                       semimajor,
                                                       expansion_order)
    if frequencies_fname is not None:
        for quantity, label in [('eccentricity', '$P_{tide}: e$'),
                                #('semimajor', '$P_{tide}: a$'),
                                #('spin', r'$P_{tide}: \Omega_\star$'),
                                ('orbital_freq', '$P_{orb}$'),
                                ('spin_freq', r'$P_\star$')]:
            pyplot.semilogy(
                eval_e,
                1.0 / numpy.absolute(leading_frequencies[quantity]),
                '-',
                label=label,
                linewidth=5
            )
        pyplot.legend()
        pyplot.ylim(0.01, 10)
        pyplot.xlim(min(eval_e), max(eval_e))
        pyplot.grid(which='both')
        pyplot.xlabel('Eccentricity')
        pyplot.ylabel('Scaled Period')
        if frequencies_fname:
            pyplot.savefig(frequencies_fname)
        else:
            pyplot.show()

    if multipliers_fname:
        for quantity, label_end in [('eccentricity', 'e')]:#,
#                                ('semimajor', 'a'),
#                                ('spin', r'\Omega_\star')]:
            sign = numpy.empty(eval_e.shape)
            sign[leading_frequencies[quantity + '_wspin_mult'] >= 0] = 1
            sign[leading_frequencies[quantity + '_wspin_mult'] < 0] = -1
            pyplot.plot(
                eval_e,
                sign * leading_frequencies[quantity + '_worb_mult'],
                'o',
                label='$m_{orb}: ' + label_end + '$',
                linewidth=5
            )
            pyplot.plot(
                eval_e,
                sign * leading_frequencies[quantity + '_wspin_mult'],
                'o',
                label='$m_{spin}: ' + label_end + '$',
                linewidth=5
            )

        pyplot.legend()
        pyplot.ylim(0, 30)
        pyplot.xlim(min(eval_e), max(eval_e))
        pyplot.grid(which='both')
        pyplot.xlabel('Eccentricity')
        pyplot.ylabel('Dominant multipliers')
        if multipliers_fname:
            pyplot.savefig(multipliers_fname)
        else:
            pyplot.show()

#pylint: disable=too-many-branches
def main(config):
    """Avoid polluting global namespace."""

    orbital_evolution_library.prepare_eccentricity_expansion(
        config.eccentricity_expansion_fname.encode('ascii'),
        1e-4,
        True,
        True
    )

    if config.frequencies_vs_e is not None or config.top_terms_vs_e is not None:
        plot_leading_term_vs_e(
            eval_e=numpy.linspace(0, 0.8, 8001),
            frequencies_fname=config.frequencies_vs_e,
            multipliers_fname=config.top_terms_vs_e,
            wspin=config.spin_frequency,
            semimajor=config.semimajor_axis
        )

    if config.single_term_rates_vs_spin is not None:
        plot_single_term_rates_vs_spin()
        if config.plot_single_term_rates_vs_spin:
            pyplot.savefig(config.plot_single_term_rates_vs_spin)
        else:
            pyplot.show()

    for eccentricity in config.rate_spectra_eccentricities:
        plot_rate_spectra(
            eccentricity,
            max_phase_lag=calc_phase_lag(config.min_lgQ),
            tidal_powerlaw=config.Q_period_powerlaw,
            break_frequency=2.0 * numpy.pi / config.Q_break_period
        )
    if config.rate_spectra_eccentricities:
        if config.rate_spectra:
            pyplot.savefig(config.rate_spectra)
        else:
            pyplot.show()

    if config.pseudo_synchronization is not None:
        pyplot.subplot(221)
        plot_pseudosynchronization(spin_scale='orbit')

        pyplot.subplot(222)
        plot_pseudosynchronization(spin_scale='ps')

        pyplot.subplot(223)
        plot_pseudosynchronization(spin_scale='ps', log=True)

        pyplot.subplot(224)
        plot_pseudosynchronization(spin_scale='orbit', log=True)

        if config.pseudo_synhcronization:
            pyplot.savefig(config.pseudo_synhcronization)
        else:
            pyplot.show()
#pylint: enable=too-many-branches
