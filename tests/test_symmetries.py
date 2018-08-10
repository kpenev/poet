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

    def _compare_evolutions(self, evolution1, evolution2, flip=True):
        """Require the two evolutions match for a, and primary Lconv & Lrad."""

        for quantity_name, quantity_name_backup in [
                ('semimajor', 'semimajor'),
                ('envelope_angmom', 'primary_envelope_angmom'),
                ('core_angmom', 'primary_core_angmom')
        ]:
            with self.subTest(quantity=quantity_name, flipped=not flip):
                quantity = [
                    getattr(evol, quantity_name,
                            getattr(evol, quantity_name_backup, None))
                    for evol in [evolution1, evolution2]
                ]
                acceptable_ages = scipy.logical_and(
                    (evolution1.age[:-1] - evolution1.age[1:]) < 0,
                    scipy.isfinite(quantity[0][:-1])
                )
                interp_quantity = InterpolatedUnivariateSpline(
                    evolution1.age[:-1][acceptable_ages],
                    quantity[0][:-1][acceptable_ages]
                )

                max_difference = scipy.nanmax(
                    scipy.fabs(
                        quantity[1][:-1] - interp_quantity(evolution2.age[:-1])
                    )
                )
                max_error = max(
                    (
                        scipy.nanmax(quantity[1][:-1])
                        -
                        scipy.nanmin(quantity[1][:-1])
                    ) * 5e-3,
                    1e-10 * scipy.nanmean(quantity[1])
                )
                self.assertLessEqual(max_difference, max_error)
        if flip:
            self._compare_evolutions(evolution2, evolution1, False)

    def test_star_planet_swap(self):
        """Compare planet-non dissipative star to 2 non-dissipative stars."""

        tdisk = 5e-3

        for star_phase_lag in [
                0.0,
                phase_lag(6.0)
        ]:
            for wind in [
                    False,
                    True
            ]:
                with self.subTest(wind=wind, phase_lag=star_phase_lag):

                    max_age = 8.7 if (star_phase_lag and wind) else 10.0

                    star = self._create_star(
                        convective_phase_lag=star_phase_lag,
                        mass=1.0,
                        wind=wind
                    )
                    planet = self._create_planet(0.8)
                    binary = self._create_binary_system(
                        star,
                        planet,
                        disk_lock_frequency=2.0 * scipy.pi / 3.0,
                        initial_semimajor=10.0,
                        disk_dissipation_age=tdisk
                    )
                    binary.evolve(max_age, 1e-3, 1e-6, None)
                    star_planet_evolution = binary.get_evolution()

                    planet.delete()
                    star.delete()
                    binary.delete()

                    #pylint false positive, members added after construction.
                    #pylint: disable=no-member
                    tdisk_index = star_planet_evolution.age.searchsorted(tdisk)
                    #pylint: enable=no-member

                    primary = self._create_star(
                        convective_phase_lag=star_phase_lag,
                        mass=1.0,
                        wind=wind
                    )
                    secondary = self._create_star(convective_phase_lag=0.0,
                                                  mass=0.8,
                                                  wind=False)
                    binary = self._create_binary_system(
                        primary,
                        secondary,
                        disk_lock_frequency=2.0 * scipy.pi / 3.0,
                        initial_semimajor=10.0,
                        disk_dissipation_age=tdisk,
                        secondary_angmom=scipy.array([
                            #pylint false positive
                            #pylint: disable=no-member
                            star_planet_evolution.envelope_angmom[tdisk_index],
                            star_planet_evolution.core_angmom[tdisk_index]
                            #pylint: enable=no-member
                        ])
                    )
                    binary.evolve(max_age, 1e-3, 1e-6, None)
                    star_star_evolution = binary.get_evolution()

                    primary.delete()
                    secondary.delete()
                    binary.delete()

                    self._compare_evolutions(star_planet_evolution,
                                             star_star_evolution)

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
