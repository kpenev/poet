#!/usr/bin/env python3

"""An interface to the POET orbital evolution library."""

import sys
sys.path.append('..')

from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

from orbital_evolution.c_interface_util import ndpointer_or_null
from basic_utils import Structure, semimajor, orbital_frequency
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

class c_dissipating_body_p(c_void_p) : pass

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

class DissipatingBody :
    """A base class for any body in a POET system."""

    def configure(self,
                  age,
                  companion_mass,
                  semimajor,
                  eccentricity,
                  spin_angmom,
                  inclination,
                  periapsis,
                  locked_surface,
                  zero_outer_inclination,
                  zero_outer_periapsis) :
        """
        Tell the body what orbit it is in.

        Args:
            - age:
                The age to set the body to.

            - companion_mass:
                The mass of the second body in the system.

            - semimajor:
                The semimajor axis of the orbit in solar radii.

            - eccentricity:
                The eccentricity of the orbit.

            - spin_angmom:
                The spin angular momenta of the non-locked zones of the body
                (outermost zone to innermost).

            - inclination:
                The inclinations of the zones of the body (same order as 
                spin_angmom).

            - periapsis:
                The arguments of periapsis of the zones of the bodies (same
                order as spin_angmom).

            - locked_surface:
                If true, the outermost zone's spin is assumed locked to a
                disk and spin_angmom is assumed to start from the next zone.

            - zero_outer_inclination:
                If true, the outermost zone's inclination is assumed to be
                zero and the inclination argument is assumed to start from
                the next zone.

            - zero_outer_periapsis:
                If true, the outermost zone's periapsis is assumed to be zero
                and the inclination argument is assumed to start from the
                next zone.

        Returns: None.
        """

        self.lib_configure_body(self.c_body,
                                age,
                                companion_mass,
                                semimajor,
                                eccentricity,
                                spin_angmom,
                                inclination,
                                periapsis,
                                locked_surface,
                                zero_outer_inclination,
                                zero_outer_periapsis)

class Binary :
    """A class for binaries POET can evolve."""

    evolution_modes = ['LOCKED_SURFACE_SPIN',
                       'BINARY',
                       'SINGLE',
                       'TABULATION']
    evolution_mode_ids = {
        mode: c_int.in_dll(library, mode + '_EVOL_MODE').value
        for mode in evolution_modes
    }

    evolution_quantities = ['age',
                            'semimajor',
                            'eccentricity',
                            'envelope_inclination',
                            'core_inclination',
                            'envelope_periapsis',
                            'core_periapsis',
                            'envelope_angmom',
                            'core_angmom',
                            'evolution_mode',
                            'wind_saturation']

    @staticmethod
    def evolution_quantity_c_type(quantity) :
        """Return the ctypes type of the given evolution quantity."""

        if quantity == 'evolution_mode' : return c_int
        elif quantity == 'wind_saturation' : return c_bool
        else : return c_double

    def __init__(self,
                 primary,
                 secondary,
                 disk_lock_frequency,
                 disk_dissipation_age,
                 initial_semimajor = None,
                 initial_orbital_period = None,
                 initial_eccentricity = 0.0,
                 initial_inclination = 0.0,
                 secondary_formation_age = None) :
        """
        Create a binary system out of two bodies.

        Args:
            - primary:
                The first body in the system. Assumed to always be there, so
                for a star-planet system this should be the star.

            - secondary:
                The second body in the system, initially may not be there and
                later may be engulfed by the first body.

            - disk_lock_frequency:
                Frequency of the surface spin of the primary when disk is
                present in rad/day.

            - disk_dissipation_age:
                Age when disk dissipates in Gyrs.

            - initial_semimajor:
                The semimajor axis of the orbit at which the secondary forms
                in solar radii. If omitted, initial_orbital_period must be
                specified.

            - initial_orbital_period:
                Alternative to specifying the initial semimajor axis.

            - initial_eccentricity:
                The eccentricity of the orbit at which the secondary forms.

            - initial_inclination:
                Inclination between surface zone of primary and initial orbit
                in radians.

            - secondary_formation_age:
                Age when the secondary forms.

        Returns: None
        """

        assert(isinstance(primary, EvolvingStar))
        self.primary = primary
        self.secondary = secondary
        if initial_semimajor is None :
            initial_semimajor = self.semimajor(initial_orbital_period)
        if secondary_formation_age is None :
            secondary_formation_age = disk_dissipation_age
        if isinstance(secondary, LockedPlanet) :
            c_create_func = library.create_star_planet_system
        else :
            assert(isinstance(secondary, EvolvingStar))
            c_create_func = library.create_star_star_system
           
        self.c_binary = c_create_func(
            primary.c_body,
            secondary.c_body,
            initial_semimajor,
            initial_eccentricity,
            initial_inclination,
            disk_lock_frequency,
            disk_dissipation_age,
            secondary_formation_age
        )

    def delete(self) :
        """Destroy the binary created at construction."""

        library.destroy_binary(self.c_binary)
        if hasattr(self, 'c_solver') :
            library.destroy_solver(self.c_solver)

    def configure(self,
                  age,
                  semimajor,
                  eccentricity,
                  spin_angmom,
                  inclination,
                  periapsis,
                  evolution_mode) :
        """
        Set the current state (orbit) of a system.

        Args:
            - age:
                The age to set the system to.

            - semimajor:
                The semimajor axis of the orbit in solar radii.

            - eccentricity:
                The eccentricity of the orbit.

            - spin_angmom:
                The spin angular momenta of the zones of the bodies (body 1
                first, outermost zone to innermost, followed by body 2).

            - inclination:
                The inclinations of the zones of the bodies (same order as
                spin_angmom). The surface zone inclination must be omitted
                for single body systems.

            - periapsis:
                The arguments of periapsis of the zones of the bodies (same
                order as spin_angmom, but not including the surface zone of
                the first body).

            - evolution_mode:
                The evolution mode to assume. Must be one of the constants
                defined.

        Returns: None
        """

        library.configure_system(self.c_binary,
                                 age,
                                 semimajor,
                                 eccentricity,
                                 spin_angmom,
                                 inclination,
                                 periapsis,
                                 self.evolution_mode_ids[evolution_mode])

    def evolve(self,
               final_age,
               max_time_step,
               precision,
               required_ages) :
        """
        Evolve the system forward from its current state.

        Args:
            - final_age:
                The age at which to stop the evolution in Gyrs. The starting
                age must be already set for the system through configure.

            - max_time_step:
                The maximum size of the time step allowed in Gyrs.

            - precision:
                The precision to require of the solution.

            - required_ages:
                Ages at which the evolution must stop precisely.

        Returnns: None
        """

        if hasattr(self, 'c_solver') :
            library.destroy_solver(self.c_solver)
        self.c_solver = library.evolve_system(
            self.c_binary,
            final_age,
            max_time_step,
            precision,
            required_ages,
            (0 if required_ages is None else required_ages.size)
        )
        self.num_evolution_steps = library.num_evolution_steps(self.c_solver)

    def get_evolution(self, quantities = evolution_quantities) :
        """
        Return the last calculated evolution.
        
        Args:
            - quantities:
                An iterable of quantities to read the evolution of. The
                evolution of omitted quantities can still be obtained later
                by subsequent calls to this method. The allowed entries are
                in the evolution_quantities property.
        
        Returns: 
            A structure with mebers named the same way as the input list of
            quantities containing the values of the corresponding quantity at
            each evolution step. The order is always in increasing age.
        """

        result = Structure()
        for q in quantities :
            setattr(result,
                    q, 
                    numpy.empty(self.num_evolution_steps,
                                dtype = self.evolution_quantity_c_type(q)))

        library.get_evolution(
            self.c_solver,
            self.c_binary,
            self.primary.c_body,
            *[getattr(result, quantity, None)
              for quantity in self.evolution_quantities]
        )
        return result

    def final_state(self) :
        """Return the final evolution state of a system (all quantities)."""

        result = Structure()
        library_final_state = [self.evolution_quantity_c_type(q)()
                               for q in self.evolution_quantities]
        library.get_final_state(self.c_solver,
                                self.c_binary,
                                self.primary.c_body,
                                *library_final_state)
        for quantity, library_value in zip(self.evolution_quantities,
                                           library_final_state) :
            setattr(result, quantity, library_value.value)
        return result

    def orbital_frequency(self, semimajor) :
        """
        The orbital frequency of the system for the given semimajor axis.

        Args:
            - semimajor:
                The semimajor axis at which the system's orbital period is
                required in solar radii. 

        Returns:
            The orbital period in days if the two bodies of this system are
            in an orbit with the given semimajor axis.
        """

        return orbital_frequency(self.primary.mass,
                                 self.secondary.mass,
                                 semimajor)

    def orbital_period(self, semimajor) :
        """
        The orbital period of the system for the given semimajor axis.

        Args:
            - semimajor:
                The semimajor axis at which the system's orbital period is
                required in solar radii.

        Returns:
            The orbital period in days if the two bodies of this system are
            in an orbit with the given semimajor axis.
        """

        return 2.0 * numpy.pi / self.orbital_frequency(semimajor)

    def semimajor(self, orbital_period) :
        """
        The semimajor axis of the system for the given orbital period.

        Args:
            - orbital_period:
                The orbital period at which the system's orbital period is
                required in days.

        Returns:
            The semimajor axis in solar radii if the two bodies of this
            system are in an orbit with the given period.
        """

        return semimajor(self.primary.mass,
                         self.secondary.mass,
                         orbital_period)
