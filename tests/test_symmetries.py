#!/usr/bin/env python3

"""Test some evolution symmetries (i.e. exchangi primary<->secondary) etc."""

import matplotlib
matplotlib.use('TkAgg')

#Importing matplotlib and setting backend needs to be first.
#pylint: disable=wrong-import-position
#pylint: disable=wrong-import-order
from matplotlib import pyplot
import os.path
import unittest

import scipy
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy import constants

from add_poet_to_path import poet_root
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.binary import Binary
from orbital_evolution.planet_interface import LockedPlanet
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.transformations import phase_lag
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
#pylint: enable=wrong-import-position
#pylint: enable=wrong-import-order

class TestSymmetries(unittest.TestCase):
    """Unit tests based on expected symmetries in the evolution."""

    #pylint false positive, constants members exist but are not detected.
    #pylint: disable=no-member
    @staticmethod
    def _create_planet(mass=(constants.M_jup / constants.M_sun).to('')):
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(
            mass=mass,
            radius=(constants.R_jup / constants.R_sun).to('')
        )
        return planet
    #pylint: enable=no-member

    def _create_star(self, convective_phase_lag, mass, wind=True):
        """Create the star to use in the evolution."""

        star = EvolvingStar(mass=mass,
                            metallicity=0.0,
                            wind_strength=(0.17 if wind else 0.0),
                            wind_saturation_frequency=2.78,
                            diff_rot_coupling_timescale=5.0e-3,
                            interpolator=self.interp)
        star.set_dissipation(zone_index=0,
                             tidal_frequency_breaks=None,
                             spin_frequency_breaks=None,
                             tidal_frequency_powers=scipy.array([0.0]),
                             spin_frequency_powers=scipy.array([0.0]),
                             reference_phase_lag=convective_phase_lag)
        star.set_dissipation(zone_index=1,
                             tidal_frequency_breaks=None,
                             spin_frequency_breaks=None,
                             tidal_frequency_powers=scipy.array([0.0]),
                             spin_frequency_powers=scipy.array([0.0]),
                             reference_phase_lag=0.0)
        return star

    @staticmethod
    def _create_binary_system(primary,
                              secondary,
                              *,
                              disk_lock_frequency,
                              initial_semimajor,
                              disk_dissipation_age,
                              secondary_angmom=None):
        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet):
            secondary_config = dict(spin_angmom=scipy.array([0.0]),
                                    inclination=None,
                                    periapsis=None)
        else:
            secondary.select_interpolation_region(disk_dissipation_age)
            secondary_config = dict(spin_angmom=secondary_angmom,
                                    inclination=scipy.array([0.0]),
                                    periapsis=scipy.array([0.0]))
        secondary.configure(age=disk_dissipation_age,
                            companion_mass=primary.mass,
                            semimajor=initial_semimajor,
                            eccentricity=0.0,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)
        if isinstance(secondary, EvolvingStar):
            secondary.detect_stellar_wind_saturation()

        primary.select_interpolation_region(primary.core_formation_age())
        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_semimajor=initial_semimajor,
                        initial_eccentricity=0.0,
                        initial_inclination=0.0,
                        disk_lock_frequency=disk_lock_frequency,
                        disk_dissipation_age=disk_dissipation_age,
                        secondary_formation_age=disk_dissipation_age)

        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=scipy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        primary.detect_stellar_wind_saturation()
        return binary

    def _plot_evolution(self,
                        binary,
                        evolution,
                        *,
                        orbit='-r',
                        primary_core='-b',
                        primary_envelope='-g',
                        secondary_core='-c',
                        secondary_envelope='-m'):
        """
        Plot the given evolution using several subplots.

        Args:
            binary:    The binary system which was used to calculate the
                evolution.

            evolution:    The calculated evolution.

            orbit:    The line style for plotting the orbital frequency.

            primary_core:    The lyne style for plotting the spin frequency of
                the core of the primary.

            primary_envelope:    The lyne style for plotting the spin frequency of
                the envelope of the primary.

            secondary_core:    The lyne style for plotting the spin frequency of
                the core of the secondary.

            secondary_envelope:    The lyne style for plotting the spin frequency of
                the envelope of the secondary.

        Returns:
            None
        """

        self._plotted = True
        wsun = 2.0 * scipy.pi / 25.34


        worb = (2.0 * scipy.pi / binary.orbital_period(evolution.semimajor)
                /
                wsun)
        wenv = (getattr(evolution, 'envelope_angmom',
                        getattr(evolution, 'primary_envelope_angmom', None))
                /
                binary.primary.envelope_inertia(evolution.age)) / wsun
        wcore = (getattr(evolution, 'core_angmom',
                         getattr(evolution, 'primary_core_angmom', None))
                 /
                 binary.primary.core_inertia(evolution.age)) / wsun

        pyplot.loglog(evolution.age, worb, orbit)
        pyplot.loglog(evolution.age, wenv, primary_envelope)
        pyplot.loglog(evolution.age, wcore, primary_core)
        pyplot.axhline(y=binary.primary.wind_saturation_frequency / wsun,
                       color='black')

        if isinstance(binary.secondary, EvolvingStar):
            wenv = (evolution.secondary_envelope_angmom
                    /
                    binary.secondary.envelope_inertia(evolution.age)) / wsun
            wcore = (evolution.secondary_core_angmom
                     /
                     binary.secondary.core_inertia(evolution.age)) / wsun

            pyplot.loglog(evolution.age, wenv, secondary_envelope)
            pyplot.loglog(evolution.age, wcore, secondary_core)

    #As simplified as I could make it.
    #pylint: disable=too-many-locals
    def _compare_evolutions(self,
                            evolution1,
                            evolution2,
                            *,
                            primary_to_secondary=False,
                            flip=True,
                            min_age=-scipy.inf,
                            max_age=scipy.inf):
        """
        Require the two evolutions match for a, and primary Lconv & Lrad.

        Args:
            evolution1:    The first of the two evolutions to compare.

            evolution2:    The second of the two evolutions to compare.

            primary_to_secondary:    If False compares the evolution of the
                primary object in both evolutions, otherwise compares the
                evolution of the primary in evolution1 to the evolution of the
                secondary in evolution2.

            flip:    If true, also calls itself with the two evolutions swapped.

            min_age:    Discrepancies at ages before this are ignored.

            min_age:    Discrepancies at ages after this are ignored.
        """

        def output_failing(data):
            """Output to stdout the discrepant evolutions."""

            for i in range(max(d.shape[0] for d in data)):
                print(
                    (
                        '%25.16e %25.16e' % tuple(data[0][i])
                        if i < data[0].shape[0] else
                        (51 * ' ')
                    )
                    +
                    ' '
                    +
                    (
                        '%25.16e %25.16e' % tuple(data[1][i])
                        if i < data[1].shape[0] else
                        (51 * ' ')
                    )
                )


        for quantity_name in [
                ('semimajor', 'semimajor'),
                ('envelope_angmom', 'primary_envelope_angmom'),
                ('core_angmom', 'primary_core_angmom')
        ]:
            with self.subTest(quantity=quantity_name[0], flipped=not flip):
                if primary_to_secondary:
                    quantity = [
                        getattr(
                            evolution1, quantity_name[0],
                            getattr(evolution1, quantity_name[1], None)
                        ),
                        getattr(
                            evolution2,
                            (
                                (
                                    '' if quantity_name[0] == 'semimajor'
                                    else 'secondary_'
                                )
                                +
                                quantity_name[0]
                            )
                        )
                    ]
                else:
                    quantity = [
                        getattr(evol, quantity_name[0],
                                getattr(evol, quantity_name[1], None))
                        for evol in [evolution1, evolution2]
                    ]
                age_within_range = [
                    scipy.logical_and(evol.age[:-1] > min_age,
                                      evol.age[:-1] < max_age)
                    for evol in [evolution1, evolution2]
                ]

                acceptable_ages = scipy.logical_and(
                    age_within_range[0],
                    scipy.logical_and(
                        (evolution1.age[:-1] - evolution1.age[1:]) < 0,
                        scipy.isfinite(quantity[0][:-1])
                    )
                )
                interp_quantity = InterpolatedUnivariateSpline(
                    evolution1.age[:-1][acceptable_ages],
                    quantity[0][:-1][acceptable_ages]
                )

                max_difference = scipy.nanmax(
                    scipy.fabs(
                        quantity[1][:-1][age_within_range[1]]
                        -
                        interp_quantity(evolution2.age[:-1][age_within_range[1]])
                    )
                )
                max_error = max(
                    (
                        scipy.nanmax(quantity[1][:-1][age_within_range[1]])
                        -
                        scipy.nanmin(quantity[1][:-1][age_within_range[1]])
                    ) * 5e-3,
                    1e-10 * scipy.nanmean(quantity[1][:-1][age_within_range[1]])
                )
                if max_difference > max_error:
                    output_failing(
                        [
                            scipy.dstack(
                                (
                                    evolution1.age[:-1][age_within_range[0]],
                                    quantity[0][:-1][age_within_range[0]]
                                )
                            )[0],
                            scipy.dstack(
                                (
                                    evolution2.age[:-1][age_within_range[1]],
                                    quantity[1][:-1][age_within_range[1]]
                                )
                            )[0]
                        ]
                    )
                self.assertLessEqual(max_difference, max_error)

        if flip:
            self._compare_evolutions(evolution2,
                                     evolution1,
                                     primary_to_secondary=primary_to_secondary,
                                     min_age=min_age,
                                     max_age=max_age,
                                     flip=False)
    #pylint: enable=too-many-locals

    def test_star_planet_swap(self):
        """Compare evolutions with secondary planet or non-dissipative star."""

        tdisk = 5e-3

        for primary_mass in [
                1.0,
                0.8
        ]:
            for secondary_mass in [
                    1.0,
                    0.8
            ]:
                for star_phase_lag in [
                        0.0,
                        phase_lag(6.0)
                ]:
                    for wind in [
                            False,
                            True
                    ]:
                        with self.subTest(wind=wind,
                                          phase_lag=star_phase_lag,
                                          primary_mass=primary_mass,
                                          secondary_mass=secondary_mass):

                            star = self._create_star(
                                convective_phase_lag=star_phase_lag,
                                mass=primary_mass,
                                wind=wind
                            )
                            planet = self._create_planet(secondary_mass)
                            binary = self._create_binary_system(
                                star,
                                planet,
                                disk_lock_frequency=2.0 * scipy.pi / 3.0,
                                initial_semimajor=10.0,
                                disk_dissipation_age=tdisk
                            )
                            binary.evolve(10.0, 1e-2, 1e-6, None)
                            star_planet_evolution = binary.get_evolution()

                            planet.delete()
                            star.delete()
                            binary.delete()

                            #pylint false positive, members added after
                            #construction.
                            #pylint: disable=no-member
                            tdisk_index = (
                                star_planet_evolution.age.searchsorted(tdisk)
                            )
                            #pylint: enable=no-member

                            primary = self._create_star(
                                convective_phase_lag=star_phase_lag,
                                mass=primary_mass,
                                wind=wind
                            )
                            secondary = self._create_star(
                                convective_phase_lag=0.0,
                                mass=secondary_mass,
                                wind=False
                            )
                            binary = self._create_binary_system(
                                primary,
                                secondary,
                                disk_lock_frequency=2.0 * scipy.pi / 3.0,
                                initial_semimajor=10.0,
                                disk_dissipation_age=tdisk,
                                secondary_angmom=scipy.array([
                                    #pylint false positive
                                    #pylint: disable=no-member
                                    star_planet_evolution.envelope_angmom[
                                        tdisk_index
                                    ] * 1e-3,
                                    star_planet_evolution.core_angmom[
                                        tdisk_index
                                    ] * 1e-3
                                    #pylint: enable=no-member
                                ])
                            )
                            binary.evolve(10.0, 1e-2, 1e-6, None)
                            star_star_evolution = binary.get_evolution()

                            primary.delete()
                            secondary.delete()
                            binary.delete()

                            self._compare_evolutions(
                                star_planet_evolution,
                                star_star_evolution,
                                max_age=(8.49 if (star_phase_lag and wind)
                                         else scipy.inf)
                            )

    def test_primary_secondary_swap(self):
        """Compare evolutions of two stars, swapping wich one is primary."""

        def get_disk_angmom(star_mass, star_phase_lag, wind):
            """Return core & envelope angular momenta at disk dissipation."""

            star = self._create_star(convective_phase_lag=star_phase_lag,
                                     mass=star_mass,
                                     wind=wind)
            planet = self._create_planet(1.0)
            binary = self._create_binary_system(
                star,
                planet,
                disk_lock_frequency=2.0 * scipy.pi / 3.0,
                initial_semimajor=10.0,
                disk_dissipation_age=config['tdisk']
            )
            binary.evolve(config['tdisk'], 1e-2, 1e-6, None)
            disk_dissipation_state = binary.final_state()

            planet.delete()
            star.delete()
            binary.delete()

            #Pylint false positive, members added after construction.
            #pylint:disable=no-member
            return scipy.array([disk_dissipation_state.envelope_angmom,
                                disk_dissipation_state.core_angmom])
            #pylint:enable=no-member

        def get_evolution(primary_mass, secondary_mass, tdisk, star_phase_lag, wind):
            """Return the evolution for the given parameters."""

            primary = self._create_star(
                convective_phase_lag=star_phase_lag,
                mass=primary_mass,
                wind=wind
            )
            secondary = self._create_star(
                convective_phase_lag=star_phase_lag,
                mass=secondary_mass,
                wind=wind
            )
            binary = self._create_binary_system(
                primary,
                secondary,
                disk_lock_frequency=2.0 * scipy.pi / 3.0,
                initial_semimajor=10.0,
                disk_dissipation_age=tdisk,
                secondary_angmom=get_disk_angmom(secondary_mass,
                                                 star_phase_lag,
                                                 wind)
            )
            binary.evolve(10.0, 1e-2, 1e-6, None)
            evolution = binary.get_evolution()

            primary.delete()
            secondary.delete()
            binary.delete()

            return evolution

        config = dict(mass1=1.0, mass2=0.8, tdisk=5e-3)

        for star_phase_lag in [
                0.0,
                phase_lag(6.0)
        ]:
            for wind in [
                    False,
                    True
            ]:
                with self.subTest(wind=wind, phase_lag=star_phase_lag):
                    evolution1 = get_evolution(
                        primary_mass=config['mass1'],
                        secondary_mass=config['mass2'],
                        tdisk=config['tdisk'],
                        star_phase_lag=star_phase_lag,
                        wind=wind
                    )
                    evolution2 = get_evolution(
                        primary_mass=config['mass2'],
                        secondary_mass=config['mass1'],
                        tdisk=config['tdisk'],
                        star_phase_lag=star_phase_lag,
                        wind=wind
                    )
                    self._compare_evolutions(
                        evolution1,
                        evolution2,
                        primary_to_secondary=True,
                        min_age=config['tdisk'],
                        max_age=(4.59 if (star_phase_lag and wind)
                                 else scipy.inf)
                    )

    @classmethod
    def setUpClass(cls):
        """Read a serialized interpolator in self.interp."""

        cls.interp = StellarEvolutionManager(
            os.path.join(poet_root, 'stellar_evolution_interpolators')
        ).get_interpolator_by_name('default')

    def setUp(self):
        """Reset plottin flag."""

        self._plotted = False

    def tearDown(self):
        """If plotting flag was set, show the plots."""

        if self._plotted:
            pyplot.show()

if __name__ == '__main__':
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        os.path.join(poet_root,
                     'eccentricity_expansion_coef.txt').encode('ascii')
    )
    unittest.main()
