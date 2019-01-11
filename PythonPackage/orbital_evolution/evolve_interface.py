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
        c_bool
    ]
    result.evolve_system.restype = c_solver_p

    result.destroy_solver.argtypes = [result.evolve_system.restype]
    result.destroy_solver.restype = None

    result.num_evolution_steps.argtypes = [result.evolve_system.restype]
    result.num_evolution_steps.restype = c_uint

    result.get_star_planet_evolution.argtypes = [
        result.evolve_system.restype,
        result.create_star_planet_system.restype,
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
        ndpointer_or_null(dtype=c_int,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_bool,
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
                          flags='C_CONTIGUOUS')
    ]
    result.get_star_star_evolution.restype = None

    result.get_star_planet_final_state.argtypes = [
        result.evolve_system.restype,
        result.create_star_planet_system.restype,
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
        POINTER(c_int),
        POINTER(c_bool)
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

    return result

library = initialize_library()
