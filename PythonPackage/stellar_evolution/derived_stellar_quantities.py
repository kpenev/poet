"""Define more quantities directly calculable from POET stellar evolution."""

from astropy import units, constants
import numpy

class TeffK :
    """Stellar effective temperature in Kelvin."""

    def __init__(self, radius, luminosity, reference_Teff = 0.0) :
        """
        Effective temperature from radius and luminosity.

        Args:
            - radius:
                A library_interface Quantity instance, giving the stellar
                radius as a function of age.
            - luminosity:
                A library_interface Quantity instance, giving the stellar
                luminosity as a function of age.
            - reference_Teff:
                The returned value is the deviation from this. Useful when
                passing to solvers or minimizers.

        Returns: None
        """

        self.radius = radius
        self.luminosity = luminosity
        self.reference = reference_Teff
        self.min_age = max(radius.min_age, luminosity.min_age)
        self.max_age = min(radius.max_age, luminosity.max_age)

    def __call__(self, age) :
        """Return the effective temperature at the given age."""

        return (
            (
                self.luminosity(age) / (4.0 * numpy.pi * self.radius(age)**2)
                *
                (
                    constants.L_sun
                    /
                    (constants.R_sun**2 * constants.sigma_sb)
                )
            )**0.25
        ).to(units.K).value - self.reference

class LogGCGS :
    """Log10 of stellar surface gravity it cgs units."""

    def __init__(self, mass, radius, reference_logg = 0.0) :
        """
        Log10 of the gravitational acceleration from mass and radius.

        Args:
            - mass:
                The mass of the star whose gravity we are interpolating.
            - radius:
                A library_interface Quantity instance, giving the stellar
                radius as a function of age.               
            - reference_logg:
                The returned value is the deviation from this. Useful when
                passing to solvers or minimizers.

        Returns: None
        """

        self.mass = mass
        self.radius = radius
        self.reference = reference_logg
        self.min_age, self.max_age = radius.min_age, radius.max_age

    def __call__(self, age) :
        """Return the log10(g) at the given age."""

        return log10(
            (
                constants.G * self.mass * constants.M_sun 
                /
                (self.radius(age) * constants.R_sun)**2
            ).to(units.cm / units.s**2).value
        ) - self.reference

class RhoCGS :
    """Mean stellar density in cgs units."""

    def __init__(self, mass, radius, reference_rho = 0.0) :
        """
        Stellar density (cgs) from mass and radius.

        Args:
            - mass:
                The mass of the star whose gravity we are interpolating.
            - radius:
                A library_interface Quantity instance, giving the stellar
                radius as a function of age.               
            - reference_rho:
                The returned value is the deviation from this. Useful when
                passing to solvers or minimizers.
        """

        self.mass = mass
        self.radius = radius
        self.reference = reference_rho
        self.min_age, self.max_age = radius.min_age, radius.max_age

    def __call__(self, age) :
        """Return the density (g/cm^3) at the given age."""

        return (
            3.0 * self.mass * constants.M_sun
            /
            (4.0 * numpy.pi * (self.radius(age) * constants.R_sun)**3)
        ).to(units.g / units.cm**3).value - self.reference
                                                                   


