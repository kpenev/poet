"""Create a plot showing the contribution of each tidal term to evolution."""

from os import path

from matplotlib import pyplot
import numpy
from astropy import units as u

from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution import SingleTermBody, LockedPlanet

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


def plot_single_term_synchronization():
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
    pyplot.show()



def main():
    """Avoid polluting global namespace."""

    e_expansion_fname = path.join(
        path.dirname(path.dirname(path.abspath(__file__))),
        'eccentricity_expansion_coef_O400.sqlite'
    )
    assert path.exists(e_expansion_fname)
    orbital_evolution_library.prepare_eccentricity_expansion(
        e_expansion_fname.encode('ascii'),
        1e-4,
        True,
        True
    )

    plot_single_term_synchronization()

    plot_rate_spectra(0.8,
                      max_phase_lag=1e-6,
                      tidal_powerlaw=1.0,
                      break_frequency=2e3 * numpy.pi)
    plot_rate_spectra(0.5,
                      max_phase_lag=1e-6,
                      tidal_powerlaw=1.0,
                      break_frequency=2e3 * numpy.pi)
    plot_rate_spectra(0.3,
                      max_phase_lag=1e-6,
                      tidal_powerlaw=1.0,
                      break_frequency=2e3 * numpy.pi)
    plot_rate_spectra(0.1,
                      max_phase_lag=1e-6,
                      tidal_powerlaw=1.0,
                      break_frequency=2e3 * numpy.pi)
    pyplot.show()

    pyplot.subplot(221)
    plot_pseudosynchronization(spin_scale='orbit')

    pyplot.subplot(222)
    plot_pseudosynchronization(spin_scale='ps')

    pyplot.subplot(223)
    plot_pseudosynchronization(spin_scale='ps', log=True)

    pyplot.subplot(224)
    plot_pseudosynchronization(spin_scale='orbit', log=True)

    pyplot.show()


if __name__ == '__main__':
    main()
