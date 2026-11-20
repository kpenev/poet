"""Run the exploratory 1.65 + 1.65 Msun, [Fe/H] = +0.5 three-zone binary."""

import argparse
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd
from mesa_reader import MesaData

from basic_utils import calc_orbital_angular_momentum
from orbital_evolution.binary import Binary
from orbital_evolution.evolve_interface import library
from orbital_evolution.star_interface import EvolvingStar
from stellar_evolution.library_interface import MESAInterpolator, library_track_fname
from stellar_evolution.MIST_realtime.build_model import convert_history_to_track


MASS = 1.65
FEH = 0.5
FLOOR = 1e-4
PRELIMINARY_AGE = 1e-3
START_WINDOW_ROWS = 10
START_WINDOW_RANGE = 0.5
SMOOTHING = 1e-8
END_AGE = 2.25
ORBITAL_PERIOD = 6.0
ECCENTRICITY = 0.1
INITIAL_SPIN_PERIODS = (8.0, 12.0)
LAGS = (1e-9, 1e-8, 1e-9)
COUPLING_TIME = 0.005
ZONES = ("convective_surface", "radiative_shell", "convective_core")
LABELS = ("Convective surface", "Radiative shell", "Convective core")
COLORS = ("#097969", "#3558a2", "#bc5629")


def grid_start_age(history):
    """Apply the default stability-window rule from more_massive_stars_grid_poet.py."""
    ages = history.star_age / 1e9
    use = ages >= PRELIMINARY_AGE
    inertias = np.column_stack(
        (history.core_inertia[use], history.env_inertia[use], history.rad_inertia[use])
    )
    for i in range(len(inertias) - START_WINDOW_ROWS + 1):
        if np.all(np.ptp(inertias[i:i + START_WINDOW_ROWS], axis=0) < START_WINDOW_RANGE):
            return float(ages[use][i])
    raise ValueError("No stable starting-age window in MESA history.")


def build_tracks(root, history, track_start):
    for quantity in ("ICONV", "IRAD"):
        MESAInterpolator.set_quantity_lower_limit(quantity, FLOOR)
    config = MESAInterpolator.get_create_interpolator_config(
        **dict.fromkeys(("ICONV", "IRAD", "MRAD", "RRAD"), {"smoothing": SMOOTHING})
    )
    tracks = []
    for zone in ZONES[1:]:
        directory = root / zone
        directory.mkdir()
        convert_history_to_track(
            MASS, str(history), directory / library_track_fname(MASS, FEH),
            min_age=track_start, zone_separated=zone
        )
        track = MESAInterpolator(mesa_dir=str(directory), num_threads=1, **config)
        track.mass_range = lambda: (MASS, MASS)
        track.feh_range = lambda: (FEH, FEH)
        tracks.append(track)
    return tracks


def make_star(tracks):
    star = EvolvingStar(
        mass=MASS, metallicity=FEH, wind_strength=0.0, wind_saturation_frequency=2.78,
        diff_rot_coupling_timescale=COUPLING_TIME, interpolator=tracks[0],
        convective_core_interpolator=tracks[1]
    )
    for index, lag in enumerate(LAGS):
        star.set_dissipation(
            zone_index=index, tidal_frequency_breaks=np.array([1.1]),
            spin_frequency_breaks=None, tidal_frequency_powers=np.array([1., 0.]),
            spin_frequency_powers=np.array([0.]), age_breaks=None,
            reference_phase_lags=np.array([lag])
        )
    return star


def evolve(tracks, start_age, timeout, end_age=END_AGE):
    stars = [make_star(tracks) for _ in INITIAL_SPIN_PERIODS]
    binary = Binary(
        primary=stars[0], secondary=stars[1], initial_orbital_period=ORBITAL_PERIOD,
        initial_eccentricity=ECCENTRICITY, initial_inclination=0.0,
        disk_lock_frequency=2 * np.pi / INITIAL_SPIN_PERIODS[0],
        disk_dissipation_age=start_age,
        secondary_formation_age=start_age
    )
    try:
        initial_momenta = []
        for star, period in zip(stars, INITIAL_SPIN_PERIODS):
            star.select_interpolation_region(start_age)
            inertia = np.array([star.zone_inertia(i, start_age) for i in range(3)])
            if not np.all(inertia >= FLOOR):
                raise RuntimeError(f"Inertia floor failed at the starting age: {inertia}")
            initial_momenta.extend(inertia * (2 * np.pi / period))
        binary.configure(
            age=start_age, semimajor=binary.semimajor(ORBITAL_PERIOD),
            eccentricity=ECCENTRICITY, spin_angmom=np.array(initial_momenta),
            inclination=np.zeros(6), periapsis=np.zeros(5), evolution_mode="BINARY"
        )
        for star in stars:
            star.detect_stellar_wind_saturation()
        ages = np.geomspace(start_age, end_age, 250)
        binary.evolve(end_age, 0.003, 1e-8, ages, timeout=timeout)
        result = binary.get_evolution(quantities=("age", "semimajor", "eccentricity"))
        period = binary.orbital_period(result.semimajor)
        series = dict(age_gyr=result.age, semimajor_rsun=result.semimajor,
                      eccentricity=result.eccentricity, orbital_period_day=period)
        total = calc_orbital_angular_momentum(MASS, MASS, result.semimajor,
                                               result.eccentricity)
        for name, star in zip(("primary", "secondary"), stars):
            angular = star.zone_angular_momentum_history(result.age.size)
            total += angular.sum(axis=1)
            for i, zone in enumerate(ZONES):
                prefix = f"{name}_{zone}"
                inertia = np.array([star.zone_inertia(i, age) for age in result.age])
                series[f"{prefix}_inertia_msun_rsun2"] = inertia
                series[f"{prefix}_spin_angmom"] = angular[:, i]
                series[f"{prefix}_spin_period_day"] = np.divide(
                    2 * np.pi * inertia, angular[:, i],
                    out=np.full_like(inertia, np.nan), where=angular[:, i] != 0
                )
        series["total_angular_momentum"] = total
    finally:
        binary.delete()
        for star in stars:
            star.delete()
    return pd.DataFrame(series)


def plot_results(data, destination):
    age = data.age_gyr
    ticks = [tick for tick in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)
             if age.iloc[0] <= tick <= age.iloc[-1]]
    title = r"$1.65+1.65\,M_\odot$ binary, $[\mathrm{Fe/H}]=+0.5$"
    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(7.4, 7.4),
                               layout="constrained")
    for ax, field, label in zip(
        axes, ("orbital_period_day", "eccentricity", "semimajor_rsun"),
        ("Orbital period (days)", "Eccentricity", r"Semimajor axis ($R_\odot$)")
    ):
        ax.plot(age, data[field], color="#3558a2", lw=1.8)
        ax.set_ylabel(label)
        ax.ticklabel_format(axis="y", useOffset=False)
    axes[-1].set_xlabel("Stellar age (Gyr)")
    figure.suptitle(title + "\nOrbital evolution")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(ticks, labels=[f"{tick:g}" for tick in ticks])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=0.2)
    figure.savefig(destination / "orbit.png", dpi=200)
    plt.close(figure)

    for field, label, filename in (
        ("inertia_msun_rsun2", r"Moment of inertia ($M_\odot R_\odot^2$)", "inertia.png"),
        ("spin_period_day", "Spin period (days)", "spin_period.png")
    ):
        figure, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10.8, 4.4),
                                   layout="constrained")
        for ax, star, initial in zip(axes, ("primary", "secondary"), INITIAL_SPIN_PERIODS):
            for zone, zone_label, color in zip(ZONES, LABELS, COLORS):
                ax.plot(age, data[f"{star}_{zone}_{field}"], label=zone_label,
                        color=color, lw=1.5)
            ax.set_title(f"{star.title()} (initial spin {initial:g} d)")
            ax.set_xlabel("Stellar age (Gyr)")
            ax.set_xscale("log")
            ax.set_xticks(ticks, labels=[f"{tick:g}" for tick in ticks])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.grid(alpha=0.2)
            if field.startswith("inertia"):
                ax.set_yscale("log")
                ax.axhline(FLOOR, color="0.35", ls=":", lw=1.0, label="Inertia floor")
            else:
                ax.plot(age, data.orbital_period_day, color="0.4", ls="--", lw=1.0,
                        label="Orbital period")
        axes[0].set_ylabel(label)
        axes[1].legend(frameon=False, fontsize=8)
        figure.suptitle(title)
        figure.savefig(destination / filename, dpi=200)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", type=Path)
    parser.add_argument("coefficients", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--timeout", type=float, default=180)
    args = parser.parse_args()
    history = MesaData(str(args.history))
    start_age = grid_start_age(history)
    earlier_ages = history.star_age[history.star_age < PRELIMINARY_AGE * 1e9]
    if not earlier_ages.size:
        raise ValueError("No MESA row before the preliminary track age.")
    previous = earlier_ages[-2] if earlier_ages.size > 1 else 0.0
    track_start = float(0.5 * (previous + earlier_ages[-1])) / 1e9
    print(f"Grid stability window starts at {start_age:.12g} Gyr.", flush=True)
    library.prepare_eccentricity_expansion(os.fsencode(args.coefficients),
                                           1e-4, True, True)
    with TemporaryDirectory() as directory:
        tracks = build_tracks(Path(directory), args.history, track_start)
        try:
            data = evolve(tracks, start_age, args.timeout)
        finally:
            for track in tracks:
                track.delete()
    print(f"POET returned {len(data)} steps ending at {data.age_gyr.iloc[-1]:.12g} Gyr.",
          flush=True)
    # An event can record the state at the same age before and after handling it.
    if len(data) < 2 or np.any(np.diff(data.age_gyr) < 0):
        raise RuntimeError("POET did not return a chronological trajectory.")
    data = data.drop_duplicates(subset="age_gyr", keep="last")
    physical_columns = data.columns.difference(data.filter(like="_spin_period_day").columns)
    invalid = ~np.isfinite(data[physical_columns]).all()
    if invalid.any():
        raise RuntimeError(f"Nonfinite evolution quantities: {invalid[invalid].index.tolist()}")
    minimum_inertia = data.filter(like="_inertia_msun_rsun2").min().min()
    if minimum_inertia < FLOOR * (1 - 1e-12):
        raise RuntimeError(f"Solver inertia dropped below the floor: {minimum_inertia}")
    total = data.total_angular_momentum
    drift = float(np.max(np.abs(total / total.iloc[0] - 1)))
    if drift > 1e-6:
        raise RuntimeError(f"Relative angular-momentum drift exceeded tolerance: {drift:g}")
    args.output.mkdir(parents=True, exist_ok=False)
    data.to_csv(args.output / "evolution.csv", index=False)
    plot_results(data, args.output)
    metadata = dict(
        history=str(args.history), coefficient_database=str(args.coefficients),
        mass_msun=MASS, feh=FEH, inertia_floor_msun_rsun2=FLOOR, start_age_gyr=start_age,
        track_min_age_gyr=track_start, start_window_rows=START_WINDOW_ROWS,
        start_window_range_msun_rsun2=START_WINDOW_RANGE,
        smoothing=dict.fromkeys(("ICONV", "IRAD", "MRAD", "RRAD"), SMOOTHING),
        requested_end_age_gyr=END_AGE, reached_age_gyr=float(data.age_gyr.iloc[-1]),
        initial_orbital_period_day=ORBITAL_PERIOD, initial_eccentricity=ECCENTRICITY,
        initial_spin_period_day=INITIAL_SPIN_PERIODS, reference_phase_lag=dict(zip(ZONES, LAGS)),
        wind_strength=0.0, differential_rotation_coupling_gyr=COUPLING_TIME,
        maximum_relative_total_angular_momentum_drift=drift, output_steps=len(data)
    )
    (args.output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
