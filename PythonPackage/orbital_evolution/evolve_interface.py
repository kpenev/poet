#!/usr/bin/env python3

"""An interface to the POET orbital evolution library."""

from ctypes import cdll, c_int, c_double, c_void_p, c_uint, c_bool
from ctypes.util import find_library
import numpy

class c_dissipating_body_p(c_void_p) : pass

class c_binary_p(c_void_p) : pass

class c_solver_p(c_void_p) : pass

def initialize_library() :
    """Prepare the orbital evolution library for use."""

    library_fname = find_library('evolve')
    if(library_fname is None) :
        raise OSError('Unable to find POET\'s evolve library.') 
    library = cdll.LoadLibrary(library_fname)

    library.create_disk_planet_system.argtypes = [c_dissipating_body_p,
                                                  c_dissipating_body_p,
                                                  c_double,
                                                  c_double,
                                                  c_double,
                                                  c_double,
                                                  c_double,
                                                  c_double]
    library.create_disk_planet_system.restype = c_binary_p

    library.destroy_binary.argtypes = [
        library.create_disk_planet_system.restype
    ]
    library.destroy_binary.restype = None

    library.configure_body.argtypes = [
        c_dissipating_body_p,
        c_double,
        c_double,
        c_double,
        c_double,
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        c_bool,
        c_bool,
        c_bool
    ]
    library.configure_body.restype = None

    library.configure_system.argtypes = [
        c_binary_p,
        c_double,
        c_double,
        c_double,
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype = c_double,
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
        numpy.ctypeslib.ndpointer(dtype = c_double,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        c_uint
    ]
    library.evolve_system.restype = c_solver_p

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

        library.configure_body(self.body,
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

    def __init__(self,
                 primary,
                 secondary,
                 initial_semimajor,
                 initial_eccentricity,
                 initial_inclination,
                 disk_lock_frequency,
                 disk_dissipation_age,
                 secondary_formation_age) :
        """
        Create a binary system out of two bodies.

        Args:
            - primary:
                The first body in the system. Assumed to always be there, so
                for a star-planet system this should be the star.

            - secondary:
                The second body in the system, initially may not be there and
                later may be engulfed by the first body.

            - initial_semimajor:
                The semimajor axis of the orbit at which the secondary forms
                in solar radii.

            - initial_eccentricity:
                The eccentricity of the orbit at which the secondary forms.

            - initial_inclination:
                Inclination between surface zone of primary and initial orbit
                in radians.

            - disk_lock_frequency:
                Frequency of the surface spin of the primary when disk is
                present in rad/day.

            - disk_dissipation_age:
                Age when disk dissipates in Gyrs.

            - secondary_formation_age:
                Age when the secondary forms.

        Returns: None
        """

        self.primary = primary
        self.secondary = secondary
        self.binary = library.create_disk_planet_system(
            primary.body,
            secondary.body,
            initial_semimajor,
            initial_eccentricity,
            initial_inclination,
            disk_lock_frequency,
            disk_dissipation_age,
            secondary_formation_age
        )

    def delete(self) :
        """Destroy the binary created at construction."""

        library.destroy_binary(self.binary)

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

        library.configure_system(self.binary,
                                 age,
                                 semimajor,
                                 eccentricity,
                                 spin_angmom,
                                 inclination,
                                 periapsis,
                                 evolution_mode)

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

        self.solver = library.evolve_system(self.binary,
                                            final_age,
                                            max_time_step,
                                            precision,
                                            required_ages,
                                            required_ages.size)

def phase_lag(lgQ) :
    """Return the phase lag corresponding to the given Q value."""

    return 15.0 / (16.0 * numpy.pi * 10.0**lgQ);
