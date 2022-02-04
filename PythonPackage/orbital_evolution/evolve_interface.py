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
    c_wchar_p,\
    c_longdouble,\
    POINTER
from ctypes.util import find_library
import numpy

from orbital_evolution.dissipating_body import c_dissipating_body_p
from orbital_evolution.c_interface_util import ndpointer_or_null
#pylint: enable=wrong-import-position

#pylint: disable=invalid-name
#pylint: disable=too-few-public-methods
class c_binary_p(c_void_p):
    """Place holder type for binary systems from the C-library."""

class c_solver_p(c_void_p):
    """Place holder type for orbit evolution solver from the C-library."""

#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

def initialize_library(library_fname=None):
    """Prepare the orbital evolution library for use."""

    if library_fname is None:
        library_fname = find_library('evolve')
    if library_fname is None:
        raise OSError('Unable to find POET\'s evolve library.')

    result = cdll.LoadLibrary(library_fname)

    result.prepare_eccentricity_expansion.argtypes = [
        c_char_p,   #filename
        c_double,   #precision
        c_bool,     #pre_load
        c_bool      #disable_precision_fali
    ]
    result.prepare_eccentricity_expansion.restype = None

    result.create_star_planet_system.argtypes = [
        c_dissipating_body_p,   #star
        c_dissipating_body_p,   #planet
        c_double,               #initial_semimajor
        c_double,               #initial_eccentricity
        c_double,               #initial_inclination
        c_double,               #disk_lock_frequency
        c_double,               #disk_dissipation_age
        c_double                #secondary_formation_age
    ]
    result.create_star_planet_system.restype = c_binary_p

    result.create_star_star_system.argtypes = (
        result.create_star_planet_system.argtypes
    )
    result.create_star_star_system.restype = (
        result.create_star_planet_system.restype
    )
    result.create_planet_planet_system.argtypes = (
        result.create_star_planet_system.argtypes[:-1]
    )
    result.create_planet_planet_system.restype = (
        result.create_star_planet_system.restype
    )

    result.destroy_binary.argtypes = [
        result.create_star_planet_system.restype
    ]
    result.destroy_binary.restype = None

    result.configure_planet.argtypes = [
        c_dissipating_body_p,                               #planet
        c_double,                                           #age
        c_double,                                           #companion_mass
        c_double,                                           #semimajor
        c_double,                                           #eccentricity
        numpy.ctypeslib.ndpointer(dtype=c_double,           #spin_angmom
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,                   #inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,                   #periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        c_bool,                                             #locked_surface
        c_bool,                                             #zero_outer_inclin.
        c_bool                                              #zero_outer_periaps.
    ]
    result.configure_planet.restype = None

    result.configure_star.argtypes = result.configure_planet.argtypes
    result.configure_star.restype = result.configure_planet.restype

    result.configure_system.argtypes = [
        c_binary_p,                                     #system
        c_double,                                       #age
        c_double,                                       #semimajor
        c_double,                                       #eccentricity
        numpy.ctypeslib.ndpointer(dtype=c_double,       #spin_angmom
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,               #inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,               #periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        c_int                                           #evolution_mode
    ]
    result.configure_system.restype = None

    result.evolve_system.argtypes = [
        c_binary_p,                             #system
        c_double,                               #final_age
        c_double,                               #max_time_step
        c_double,                               #precision
        ndpointer_or_null(dtype=c_double,       #required_ages
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        c_uint,                                 #num_required_ages
        c_bool,                                 #print_progress
        c_double                                #max_runtime
    ]
    result.evolve_system.restype = c_solver_p

    result.destroy_solver.argtypes = [result.evolve_system.restype]
    result.destroy_solver.restype = None

    result.num_evolution_steps.argtypes = [result.evolve_system.restype]
    result.num_evolution_steps.restype = c_uint

    result.get_star_planet_evolution.argtypes = [
        result.evolve_system.restype,               #solver
        result.create_star_planet_system.restype,   #system
        c_dissipating_body_p,                       #star
        c_dissipating_body_p,                       #planet
        ndpointer_or_null(dtype=c_double,           #age
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #semimajor
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #eccentricity
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #envelope inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #core_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #envelope_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #core_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #envelope_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #core_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #planet_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #planet_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #planet_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_int,              #evolution_mode
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool,             #wind_saturation
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #semimajor rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #eccentricity rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #envelope inclination rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #core_inclination rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #envelope_periapsis rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #core_periapsis rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #envelope_angmom rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #core_angmom rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #planet_inclination rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #planet_periapsis rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #planet_angmom rate
                          ndim=1,
                          flags='C_CONTIGUOUS')

    ]
    result.get_star_planet_evolution.restype = None

    result.get_star_star_evolution.argtypes = [
        result.evolve_system.restype,               #solver
        result.create_star_planet_system.restype,   #system
        c_dissipating_body_p,                       #primary
        c_dissipating_body_p,                       #secondary
        ndpointer_or_null(dtype=c_double,           #age
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #semimajor
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #eccentricity
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_envelope_inclin
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_core_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_envelope_perapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_core_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_envelope_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_core_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_envelope_inclin
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_core_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_envelope_periaps
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_core_periapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_envelope_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_core_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_int,              #evolution_mode
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool,             #primary_wind_saturation
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool,             #secondary_wind_saturation
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #semimajor_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #eccentricity_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_env_incl_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_core_incl_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_env_periapsis_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_core_periapsis_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_env_angmom_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_core_angmom_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_env_incl_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_core_incl_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_env_peri_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_core_peri_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_env_angmom_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_core_angmom_rate
                          ndim=1,
                          flags='C_CONTIGUOUS')

    ]
    result.get_star_star_evolution.restype = None

    result.get_planet_planet_evolution.argtypes = [
        result.evolve_system.restype,               #solver
        result.create_star_planet_system.restype,   #system
        c_dissipating_body_p,                       #primary
        c_dissipating_body_p,                       #secondary
        ndpointer_or_null(dtype=c_double,           #age
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #semimajor
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #eccentricity
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_perapsis
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_inclination
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_periaps
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_angmom
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_int,              #evolution_mode
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #semimajor_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #eccentricity_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_incl_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_periapsis_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #primary_angmom_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_incl_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_peri_rate
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,           #secondary_angmom_rate
                          ndim=1,
                          flags='C_CONTIGUOUS')
    ]
    result.get_star_star_evolution.restype = None


    result.get_star_planet_final_state.argtypes = [
        result.evolve_system.restype,               #solver
        result.create_star_planet_system.restype,   #system
        c_dissipating_body_p,                       #star
        c_dissipating_body_p,                       #planet
        POINTER(c_double),                          #age
        POINTER(c_double),                          #semimajor
        POINTER(c_double),                          #eccentricity
        POINTER(c_double),                          #envelope_inclination
        POINTER(c_double),                          #core_inclination
        POINTER(c_double),                          #envelovpe_peripasis
        POINTER(c_double),                          #core_periapsis
        POINTER(c_double),                          #envolepo_angmom
        POINTER(c_double),                          #core_angmom
        POINTER(c_double),                          #planet_inclination
        POINTER(c_double),                          #planet_periapsis
        POINTER(c_double),                          #planet_angmom
        POINTER(c_int),                             #evolution_mode
        POINTER(c_bool)                             #wind_saturation
    ]
    result.get_star_planet_final_state.restype = None

    result.get_star_star_final_state.argtypes = [
        result.evolve_system.restype,               #solver
        result.create_star_planet_system.restype,   #system
        c_dissipating_body_p,                       #primary
        c_dissipating_body_p,                       #secondary
        POINTER(c_double),                          #age
        POINTER(c_double),                          #semimajor
        POINTER(c_double),                          #eccentricity
        POINTER(c_double),                          #primary_env_inclination
        POINTER(c_double),                          #primary_core_inclination
        POINTER(c_double),                          #primary_env_periapsis
        POINTER(c_double),                          #primary_core_periapsis
        POINTER(c_double),                          #primary_env_angmom
        POINTER(c_double),                          #primary_core_angmom
        POINTER(c_double),                          #secondary_env_inclination
        POINTER(c_double),                          #secondary_core_inclination
        POINTER(c_double),                          #secondary_env_periapsis
        POINTER(c_double),                          #secondary_core_periapsis
        POINTER(c_double),                          #secondary_env_angmom
        POINTER(c_double),                          #secondary_core_angmom
        POINTER(c_int),                             #evolution_mode
        POINTER(c_bool),                            #primary_wind_saturation
        POINTER(c_bool)                             #secondary_wind_saturation
    ]
    result.get_star_star_final_state.restype = None

    result.get_expansion_coeff_precision.argtypes = [
        c_int,  #m
        c_int   #s
    ]
    result.get_expansion_coeff_precision.restype = c_double

    result.evaluate_expansion_coeff.argtypes = [
        c_int,      #m
        c_int,      #s
        c_double,   #e
        c_bool      #deriv
    ]
    result.evaluate_expansion_coeff.restype = c_double

    return result

library = initialize_library(
    '/home/kpenev/projects/git/poet/build/libs/evolve/shared/release/libevolve.so'
)
