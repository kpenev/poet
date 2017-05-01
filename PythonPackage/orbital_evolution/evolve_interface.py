#!/usr/bin/env python3

"""An interface to the POET orbital evolution library."""

import sys
sys.path.append('..')

from orbital_evolution.dissipating_body import c_dissipating_body_p
from orbital_evolution.c_interface_util import ndpointer_or_null
from ctypes import\
    cdll,\
    c_int,\
    c_double,\
    c_void_p,\
    c_uint,\
    c_bool,\
    c_char_p,\
    POINTER,\
    byref
from ctypes.util import find_library
import numpy
from astropy import units, constants

class c_binary_p(c_void_p) : pass

class c_solver_p(c_void_p) : pass

def initialize_library() :
    """Prepare the orbital evolution library for use."""

    library_fname = find_library('evolve')
    if(library_fname is None) :
        raise OSError('Unable to find POET\'s evolve library.') 
    library = cdll.LoadLibrary(library_fname)

    library.read_eccentricity_expansion_coefficients.argtypes = [c_char_p]
    library.read_eccentricity_expansion_coefficients.restype = None

    library.create_star_planet_system.argtypes = [c_dissipating_body_p,
                                                  c_dissipating_body_p,
                                                  c_double,
                                                  c_double,
                                                  c_double,
                                                  c_double,
                                                  c_double,
                                                  c_double]
    library.create_star_planet_system.restype = c_binary_p

    library.create_star_star_system.argtypes = (
        library.create_star_planet_system.argtypes
    )
    library.create_star_star_system.restype = (
        library.create_star_planet_system.restype
    )

    library.destroy_binary.argtypes = [
        library.create_star_planet_system.restype
    ]
    library.destroy_binary.restype = None

    library.configure_planet.argtypes = [
        c_dissipating_body_p,
        c_double,
        c_double,
        c_double,
        c_double,
        numpy.ctypeslib.ndpointer(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        c_bool,
        c_bool,
        c_bool
    ]
    library.configure_planet.restype = None

    library.configure_star.argtypes = library.configure_planet.argtypes
    library.configure_star.restype = library.configure_planet.restype

    library.configure_system.argtypes = [
        c_binary_p,
        c_double,
        c_double,
        c_double,
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        c_int
    ]
    library.configure_system.restype = None

    library.evolve_system.argtypes = [
        c_binary_p,
        c_double,
        c_double,
        c_double,
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        c_uint
    ]
    library.evolve_system.restype = c_solver_p

    library.destroy_solver.argtypes = [library.evolve_system.restype]
    library.destroy_solver.restype = None

    library.num_evolution_steps.argtypes = [library.evolve_system.restype]
    library.num_evolution_steps.restype = c_uint

    library.get_star_planet_evolution.argtypes = [
        library.evolve_system.restype,
        library.create_star_planet_system.restype,
        c_dissipating_body_p,
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_int,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_bool,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS')
    ]
    library.get_star_planet_evolution.restype = None

    library.get_star_star_evolution.argtypes = [
        library.evolve_system.restype,
        library.create_star_planet_system.restype,
        c_dissipating_body_p,
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_double,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_int,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_bool,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS'),
        ndpointer_or_null(dtype = c_bool,
                          ndim = 1,
                          flags = 'C_CONTIGUOUS')
    ]
    library.get_star_star_evolution.restype = None

    library.get_star_planet_final_state.argtypes = [
        library.evolve_system.restype,
        library.create_star_planet_system.restype,
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
    library.get_star_planet_final_state.restype = None

    library.get_star_star_final_state.argtypes = [
        library.evolve_system.restype,
        library.create_star_planet_system.restype,
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
    library.get_star_star_final_state.restype = None

    return library

library = initialize_library()
