"""Define base class for all bodies participating in evolution calc."""

from ctypes import c_void_p
from abc import ABC, abstractmethod

#Naming convention follows ctypes.
#pylint: disable=invalid-name
#Only used for type checking so no need for public methods
#pylint: disable=too-few-public-methods
class c_dissipating_body_p(c_void_p):
    """Dummy class only used for type checking."""
#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

class DissipatingBody(ABC):
    """A base class for any body in a POET system."""

    @abstractmethod
    def lib_configure_body(self, *args, **kwargs):
        """Ensure a function for configuring the body is defined."""

    def __init__(self):
        """Initialize some attributes."""

        self.dissipation = dict()
        self.spin_angmom = None

    def configure(self,
                  *,
                  age,
                  companion_mass,
                  semimajor,
                  eccentricity,
                  spin_angmom,
                  inclination,
                  periapsis,
                  locked_surface,
                  zero_outer_inclination,
                  zero_outer_periapsis):
        """
        Tell the body what orbit it is in.

        Args:
            age:    The age to set the body to.

            companion_mass:    The mass of the other body in the system.

            semimajor:    The semimajor axis of the orbit in solar radii.

            eccentricity:    The eccentricity of the orbit.

            spin_angmom:    The spin angular momenta of the non-locked zones of
                the body (outermost zone to innermost).

            inclination:    The inclinations of the zones of the body (same
                order as spin_angmom).

            periapsis:    The arguments of periapsis of the zones of the bodies
                (same order as spin_angmom).

            locked_surface:    If true, the outermost zone's spin is assumed
                locked to a disk and spin_angmom is assumed to start from the
                next zone.

            zero_outer_inclination:    If true, the outermost zone's inclination
                is assumed to be zero and the inclination argument is assumed to
                start from the next zone.

            zero_outer_periapsis:    If true, the outermost zone's periapsis is
                assumed to be zero and the inclination argument is assumed to
                start from the next zone.

        Returns:
            None
        """

        self.spin_angmom = spin_angmom
        #Objects inheriting from this class are expected to define self.c_body
        #before invoking this method
        #pylint: disable=no-member
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
        #pylint: enable=no-member

    def set_dissipation(self, *, zone_index, **dissipation_params):
        """Record the dissipation that was defined for a zone of the body."""

        self.dissipation[zone_index] = dict(dissipation_params)
