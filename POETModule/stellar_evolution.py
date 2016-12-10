#!/usr/bin/python
from ctypes import cdll, c_int, c_double

library = cdll.LoadLibrary('libstellarEvolution.so')

class MESAInterpolator :
    """A class for interpolating among a set of MESA tracks."""

    __quantity_ids = {
        q: c_int.in_dll(library, q).value
        for q in ['RADIUS', 'ICONV', 'LUM', 'IRAD', 'MRAD', 'RRAD']
    }

    def __init__(self, mesa_dir) :
        """
        Prepare a MESA based interpolation.
        
        Args:
            - mesa_dir: a directory contaning a grid (mass and metallicity)
                        of MESA tracks to base the interpolation on.
        Returns: None.
        """

        self.interpolator = library.create_interpolator(bytes(mesa_dir))

    def delete() :
        """
        Free the resources allocated at construction.
        """

        library.destroy_interpolator(self.interpolator)

    def __call__(self, quantity, mass, metallicity) :
        """
        Return a stellar quantity interpolated to the given mass and [Fe/H].

        Args:
            - quantity: A string identifying the quantity to interpolate.
                        Valid values: 'radius', 'iconv', 'lum', 'irad',
                        'mrad', 'rrad'. Case insensitive.
            - mass: The mass of the star for which this quantity should be
                    defined in solar masses.
            - metallicity: The metallicity of the star for which this quantity 
                           should be defined as [Fe/H].

        Returns: A Quantity instance, callable with an age parameter.
        """

        return Quantity(
            library.create_quantity(self.interpolator,
                                    self.__quantity_ids[quantity.upper()],
                                    c_double(mass),
                                    c_double(metallicity))
        )

class Quantity :

    def __init__(self, underlying_quantity) :
        """Wrap the underlying EvolvingStellarQuantity into a callable."""

        self.underlying_quantity = underlying_quantity

    def delete() :
        """Destroy the underlying quantity."""

        library.destroy_quantity(self.underlying_quantity)

    def __call__(self, age) :
        """Evaluate the underlying quantity at the given age (in Gyr)."""

        return library.evaluate_quantity(self.underlying_quantity,
                                         c_double(age))

if __name__ == '__main__' :
    interpolator = MESAInterpolator('../poet_src/StellarEvolution/MESA_sub')
    radius = interpolator('radius', 1.0, 0.0)
    print('R(1.0) = ' + repr(radius(1.0)))
