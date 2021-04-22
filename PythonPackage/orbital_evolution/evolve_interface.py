#!/usr/bin/env python3

"""An interface to the POET orbital evolution library."""

import sys
sys.path.append('..')

#pylint: disable=wrong-import-position
from ctypes import\
    cdll,\
    c_int,\
    c_double,\
    c_void_p,\
    c_uint,\
    c_bool,\
    c_char_p,\
    POINTER
from ctypes.util import find_library
import numpy

from orbital_evolution.c_interface_util import ndpointer_or_null
#pylint: enable=wrong-import-position

#pylint: disable=invalid-name
#pylint: disable=too-few-public-methods
class c_binary_p(c_void_p):
    """Place holder type for binary systems from the C-library."""

class c_solver_p(c_void_p):
    """Place holder type for orbit evolution solver from the C-library."""

class c_dissipating_body_p(c_void_p):
    """Dummy class only used for type checking."""

class c_dissipating_zone_p(c_void_p):
    """Dummy class only used for type checking."""

class c_broken_powerlaw_phase_lag_zone_p(c_void_p):
    """Place holder type for BrokenPowerlawPhaseLagZone from C-library."""

#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

#Define interface to c library, does not make sense to split
#pylint: disable=too-many-statements
def initialize_library():
    """Prepare the orbital evolution library for use."""

    library_fname = find_library('evolve')
    if library_fname is None:
        raise OSError('Unable to find POET\'s evolve library.')
    result = cdll.LoadLibrary(library_fname)

    result.read_eccentricity_expansion_coefficients.argtypes = [c_char_p]
    result.read_eccentricity_expansion_coefficients.restype = None

    result.create_star_planet_system.argtypes = [c_dissipating_body_p,
                                                 c_dissipating_body_p,
                                                 c_double,
                                                 c_double,
                                                 c_double,
                                                 c_double,
                                                 c_double,
                                                 c_double]
    result.create_star_planet_system.restype = c_binary_p

    result.create_star_star_system.argtypes = (
        result.create_star_planet_system.argtypes
    )
    result.create_star_star_system.restype = (
        result.create_star_planet_system.restype
    )

    result.destroy_binary.argtypes = [
        result.create_star_planet_system.restype
    ]
    result.destroy_binary.restype = None

    result.configure_planet.argtypes = [
        c_dissipating_body_p,
        c_double,
        c_double,
        c_double,
        c_double,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        c_bool,
        c_bool,
        c_bool
    ]
    result.configure_planet.restype = None

    result.configure_star.argtypes = result.configure_planet.argtypes
    result.configure_star.restype = result.configure_planet.restype

    result.configure_system.argtypes = [
        c_binary_p,
        c_double,
        c_double,
        c_double,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        c_int
    ]
    result.configure_system.restype = None

    result.evolve_system.argtypes = [
        c_binary_p,
        c_double,
        c_double,
        c_double,
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        c_uint,
        c_bool,
        c_double
    ]
    result.evolve_system.restype = c_solver_p

    result.destroy_solver.argtypes = [result.evolve_system.restype]
    result.destroy_solver.restype = None

    result.num_evolution_steps.argtypes = [result.evolve_system.restype]
    result.num_evolution_steps.restype = c_uint

    result.get_star_planet_evolution.argtypes = [
        result.evolve_system.restype, # solver
        result.create_star_planet_system.restype, #system
        c_dissipating_body_p, #star
        c_dissipating_body_p, #planet
        ndpointer_or_null(dtype=c_double, #age
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #semimajor
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #eccentricity
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #envelope inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #core_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #envelope_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #core_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #envelope_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #core_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #planet_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #planet_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #planet_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_int, #evolution_mode
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool, #wind_saturation
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #semimajor rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #eccentricity rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #envelope inclination rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #core_inclination rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #envelope_periapsis rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #core_periapsis rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #envelope_angmom rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #core_angmom rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #planet_inclination rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #planet_periapsis rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double, #planet_angmom rate
                          ndim=1,
                          flags='C_CONTIGUOUS')

    ]
    result.get_star_planet_evolution.restype = None

    result.get_star_star_evolution.argtypes = [
        result.evolve_system.restype,
        result.create_star_planet_system.restype,
        c_dissipating_body_p,
        c_dissipating_body_p,
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_int,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS')

    ]
    result.get_star_star_evolution.restype = None

    result.get_star_planet_final_state.argtypes = [
        result.evolve_system.restype, #solver
        result.create_star_planet_system.restype, #system
        c_dissipating_body_p, #star
        c_dissipating_body_p, #planet
        POINTER(c_double), #age
        POINTER(c_double), #semimajor
        POINTER(c_double), #eccentricity
        POINTER(c_double), #envelope_inclination
        POINTER(c_double), #core_inclination
        POINTER(c_double), #envelovpe_peripasis
        POINTER(c_double), #core_periapsis
        POINTER(c_double), #envolepo_angmom
        POINTER(c_double), #core_angmom
        POINTER(c_double), #planet_inclination
        POINTER(c_double), #planet_periapsis
        POINTER(c_double), #planet_angmom
        POINTER(c_int), #evolution_mode
        POINTER(c_bool) #wind_saturation
    ]
    result.get_star_planet_final_state.restype = None

    result.get_star_star_final_state.argtypes = [
        result.evolve_system.restype,
        result.create_star_planet_system.restype,
        c_dissipating_body_p,
        c_dissipating_body_p,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_int),
        POINTER(c_bool),
        POINTER(c_bool)
    ]
    result.get_star_star_final_state.restype = None

    result.orbital_angular_velocity.argtypes = [
        c_double, #m1 in Msun
        c_double, #m2 in Msun
        c_double, #semimajor in Rsun
        c_bool    #deriv w.r.t. semimajor?
    ]
    result.orbital_angular_velocity.restype = c_double

    result.orbital_energy.argtypes = [
        c_double, #m1 in Msun
        c_double, #m2 in Msun
        c_double, #semimajor in Rsun
        c_uint    #deriv order w.r.t. semimajor
    ]
    result.orbital_energy.restype = c_double

    result.orbital_angular_momentum.argtypes = [
        c_double, #m1 in Msun
        c_double, #m2 in Msun
        c_double, #semimajor in Rsun
        c_double  #eccentricity
    ]
    result.orbital_angular_momentum.restype = c_double

    result.semimajor_from_period.argtypes = [
        c_double, #m1 in Msun
        c_double, #m2 in Msun
        c_double  #period in days
    ]
    result.semimajor_from_period.restype = c_double

    result.set_zone_dissipation.argtypes = [
        c_broken_powerlaw_phase_lag_zone_p,

        c_uint, #num_tidal_frequency_breaks,

        c_uint, #num_spin_frequency_breaks,

        ndpointer_or_null(dtype=c_double, #tidal_frequency_breaks
                          ndim=1,
                          flags='C_CONTIGUOUS'),

        ndpointer_or_null(dtype=c_double, #spin_frequency_breaks
                          ndim=1,
                          flags='C_CONTIGUOUS'),

        numpy.ctypeslib.ndpointer(dtype=c_double, #tidal_frequency_powers
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),

        numpy.ctypeslib.ndpointer(dtype=c_double, #spin_frequency_powers
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),

        c_double, #reference phase lag
        c_double, #inertial mode enhancement
        c_double  #inertial mode sharpness
    ]
    result.set_zone_dissipation.restype = None

    result.configure_zone.argtypes = [
        c_dissipating_zone_p, #zone,
        c_bool,     #initialize,
        c_double,   #age,
        c_double,   #orbital_frequency,
        c_double,   #eccentricity,
        c_double,   #orbital_angmom,
        c_double,   #spin,
        c_double,   #inclination,
        c_double,   #periapsis,
        c_bool,     #spin_is_frequency,
        numpy.ctypeslib.ndpointer(dtype=c_int, #single_term
                                  shape=(2,),
                                  flags='C_CONTIGUOUS')
    ]
    result.configure_zone.restype = None

    result.get_envelope.argtypes = [c_dissipating_body_p]
    result.get_envelope.restype = result.set_zone_dissipation.argtypes[0]

    result.get_zone_tidal_power.argtypes = [c_dissipating_zone_p]
    result.get_zone_tidal_power.restype = c_double

    result.get_zone_tidal_torque_x.argtypes = [c_dissipating_zone_p]
    result.get_zone_tidal_torque_x.restype = c_double

    result.get_zone_tidal_torque_y.argtypes = [c_dissipating_zone_p]
    result.get_zone_tidal_torque_y.restype = c_double

    result.get_zone_tidal_torque_z.argtypes = [c_dissipating_zone_p]
    result.get_zone_tidal_torque_z.restype = c_double

    return result
#pylint: enable=too-many-statements

library = initialize_library()
