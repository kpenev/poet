"""Useful functions for calculating orbital evolution."""

import numpy
from astropy import constants

from orbital_evolution.binary import Binary
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

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
        print('Setting planet dissipation')
        try:
            planet.set_dissipation(tidal_frequency_breaks=None,
                                   spin_frequency_breaks=None,
                                   tidal_frequency_powers=numpy.array([0.0]),
                                   spin_frequency_powers=numpy.array([0.0]),
                                   reference_phase_lag=float(phase_lag))
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

    star = EvolvingStar(mass=mass,
                        metallicity=metallicity,
                        wind_strength=wind_strength,
                        wind_saturation_frequency=wind_saturation_frequency,
                        diff_rot_coupling_timescale=diff_rot_coupling_timescale,
                        interpolator=interpolator)
    if convective_phase_lag:
        try:
            star.set_dissipation(
                zone_index=0,
                tidal_frequency_breaks=None,
                spin_frequency_breaks=None,
                tidal_frequency_powers=numpy.array([0.0]),
                spin_frequency_powers=numpy.array([0.0]),
                reference_phase_lag=float(convective_phase_lag)
            )
        except TypeError:
            star.set_dissipation(zone_index=0,
                                 **convective_phase_lag)
    star.select_interpolation_region(star.core_formation_age()
                                     if interp_age is None else
                                     interp_age)
    return star

def create_system(primary,
                  secondary,
                  disk_lock_frequency,
                  *,
                  initial_eccentricity=0.0,
                  porb_initial=3.5,
                  disk_dissipation_age=4e-3,
                  initial_inclination=0.0,
                  secondary_formation_age=None):
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
    binary.configure(age=primary.core_formation_age(),
                     semimajor=float('nan'),
                     eccentricity=float('nan'),
                     spin_angmom=numpy.array([0.0]),
                     inclination=None,
                     periapsis=None,
                     evolution_mode='LOCKED_SURFACE_SPIN')

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
                            numpy.array([0.01, 0.01])
                            if isinstance(secondary, EvolvingStar) else
                            numpy.array([0.0])
                        ),
                        inclination=initial_obliquity,
                        periapsis=initial_periapsis,
                        locked_surface=False,
                        zero_outer_inclination=True,
                        zero_outer_periapsis=True)

    primary.detect_stellar_wind_saturation()
    if isinstance(secondary, EvolvingStar):
        secondary.detect_stellar_wind_saturation()

    return binary
