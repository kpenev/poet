"""Native three-zone smoke test. Set POET_ECCENTRICITY_DB to a populated database."""

import csv
from itertools import product
import os
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
import unittest

import add_poet_to_path
import numpy as np

from basic_utils import calc_orbital_angular_momentum
from orbital_evolution.binary import Binary
from orbital_evolution.evolve_interface import library
from orbital_evolution.planet_interface import LockedPlanet
from orbital_evolution.star_interface import EvolvingStar
from stellar_evolution.library_interface import MESAInterpolator, library_track_fname


class ThreeZoneStarTest(unittest.TestCase):
    def test_orbital_evolution(self):
        coefficients = os.environ.get("POET_ECCENTRICITY_DB")
        if not coefficients:
            self.skipTest("Set POET_ECCENTRICITY_DB to run native evolution")
        with sqlite3.connect(Path(coefficients).resolve().as_uri() + "?mode=ro",
                             uri=True) as database:
            order, total, unique = database.execute(
                "SELECT MAX(s), COUNT(*), COUNT(DISTINCT m || ',' || s) "
                "FROM interpolations WHERE m IN (-2, 0, 2) AND s >= 0"
            ).fetchone()
        self.assertEqual((total, unique), (3 * (order + 1),) * 2,
                         "Use a complete eccentricity table; (m,s) rows are missing or duplicated")
        library.prepare_eccentricity_expansion(os.fsencode(coefficients),
                                               1e-4, True, True)
        cases = ((0.03, 0.04, 0.01), (0.03, 0.04, 1e-8), (0.03, 0.04, 0.0),
                 (0.03, 0.0, 0.01), (0.0, 0.04, 0.01))
        for eccentricity, inertias in product((0.0, 0.1), cases):
            with self.subTest(eccentricity=eccentricity, inertias=inertias):
                with TemporaryDirectory() as directory:
                    self.check_evolution(Path(directory), inertias, eccentricity)

    @staticmethod
    def track(directory, surface, interior, mass, radius):
        directory.mkdir()
        with (directory / library_track_fname(1.0, 0.0)).open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(("age", "R_star", "L_star", "M_star",
                             "R_tachocline", "M_ini", "M_conv", "M_rad",
                             "I_conv", "I_rad"))
            for age in np.linspace(0.8e9, 1.2e9, 21):
                writer.writerow((age, 1, 1, 1, radius, 1, 1-mass, mass,
                                 surface, interior))
        exact = dict(smoothing=float('nan'), nodes=0,
                     vs_log_age=False, log_quantity=False)
        interpolator = MESAInterpolator(
            mesa_dir=str(directory), num_threads=1,
            **MESAInterpolator.get_create_interpolator_config(
                **dict.fromkeys(MESAInterpolator.quantity_list, exact)
            )
        )
        interpolator.mass_range = lambda: (1.0, 1.0)
        interpolator.feh_range = lambda: (0.0, 0.0)
        return interpolator

    def check_evolution(self, directory, inertias, eccentricity):
        surface, shell, core = inertias
        core_mass, core_radius = (0.2, 0.3) if core else (0.0, 0.0)
        shell_mass, shell_radius = (
            (1.0, 1.0) if not surface else (0.8, 0.7) if shell
            else (core_mass, core_radius)
        )
        shell_track = self.track(directory / "shell", surface, shell,
                                 shell_mass, shell_radius)
        core_track = self.track(directory / "core", surface, core,
                                 core_mass, core_radius)
        star = EvolvingStar(
            mass=1.0, metallicity=0.0, wind_strength=0.0,
            wind_saturation_frequency=2.78, diff_rot_coupling_timescale=0.005,
            interpolator=shell_track, convective_core_interpolator=core_track
        )
        planet = LockedPlanet(mass=0.001, radius=0.1)
        binary = Binary(primary=star, secondary=planet, initial_orbital_period=4.0,
                        initial_eccentricity=eccentricity, initial_inclination=0.0,
                        disk_lock_frequency=0.5, disk_dissipation_age=0.9)
        try:
            for index, inertia in enumerate(inertias):
                self.assertAlmostEqual(star.zone_inertia(index, 1.0), inertia)
                star.set_dissipation(
                    zone_index=index, tidal_frequency_breaks=np.array([1.1]),
                    spin_frequency_breaks=None, tidal_frequency_powers=np.array([1., 0.]),
                    spin_frequency_powers=np.array([0.]), age_breaks=None,
                    reference_phase_lags=np.array([
                        (index + 1) * 1e-6 * min(1.0, inertia / 0.01)
                        if inertia > 0 else (index + 1) * 1e-6
                    ])
                )
            binary.configure(age=1.0, semimajor=binary.semimajor(4.0),
                             eccentricity=eccentricity,
                             spin_angmom=np.r_[np.array(inertias) * 0.5, 0.],
                             inclination=np.zeros(4), periapsis=np.zeros(3),
                             evolution_mode="BINARY")
            star.detect_stellar_wind_saturation()
            binary.evolve(1.001, 0.0001, 1e-8, np.array([1.0, 1.001]), timeout=10)
            result = binary.get_evolution()
            spins = star.zone_angular_momentum_history(result.age.size)
            self.assertAlmostEqual(result.age[-1], 1.001)
            self.assertTrue(np.isfinite(spins).all())
            self.assertTrue(np.all(spins[-1, np.array(inertias) > 0]
                                   != spins[0, np.array(inertias) > 0]))
            np.testing.assert_array_equal(spins[:, np.array(inertias) == 0], 0)
            total = (calc_orbital_angular_momentum(1.0, 0.001, result.semimajor,
                                                  result.eccentricity)
                     + spins.sum(axis=1) + result.planet_angmom)
            np.testing.assert_allclose(total, total[0], atol=1e-9, rtol=1e-7)
        finally:
            binary.delete()
            planet.delete()
            star.delete()
            core_track.delete()
            shell_track.delete()


if __name__ == "__main__":
    unittest.main()
