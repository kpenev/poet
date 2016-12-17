#!/usr/bin/python
from ctypes import cdll, c_int, c_double, c_void_p, c_char_p, c_uint
import numpy

def initialize_library() :
    """Prepare the stellarEvolution library for use."""

    library = cdll.LoadLibrary('libstellarEvolution.so')

    library.create_interpolator.argtypes = [c_char_p]
    library.create_interpolator.restype = c_void_p

    library.destroy_interpolator.argtypes = [
        library.create_interpolator.restype
    ]
    library.destroy_interpolator.restype = None

    library.create_quantity.argtypes = [c_void_p, c_int, c_double, c_double]
    library.create_quantity.restype = c_void_p

    library.destroy_quantity.argtypes = [c_void_p]
    library.destroy_quantity.restype = None

    library.evaluate_quantity.argtypes = [c_void_p, c_double]
    library.evaluate_quantity.restype = c_double

    library.evaluate_quantity_array.argtypes = [
        c_void_p,
        numpy.ctypeslib.ndpointer(dtype = numpy.float64,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        c_uint
    ]

    library.quantity_min_age.restype = c_double

    library.quantity_max_age.restype = c_double

    library.save_interpolator.argtypes = [c_void_p, c_char_p]
    library.save_interpolator.restype = None

    library.load_interpolator.argtypes = [c_char_p]
    library.load_interpolator.restype = c_void_p

    library.differentiate_quantity.argtypes = [c_void_p, c_double]
    library.differentiate_quantity.restype = numpy.ctypeslib.ndpointer(
        dtype = numpy.float64,
        shape = (3,),
        flags = 'C_CONTIGUOUS'
    )

    library.differentiate_quantity_array.argtypes = [
        c_void_p,
        numpy.ctypeslib.ndpointer(dtype = numpy.float64,
                                  ndim = 1,
                                  flags = 'C_CONTIGUOUS'),
        c_uint
    ]

    return library

library = initialize_library()

class MESAInterpolator :
    """A class for interpolating among a set of MESA tracks."""

    quantity_list = ['RADIUS', 'ICONV', 'LUM', 'IRAD', 'MRAD', 'RRAD']

    quantity_ids = {q: c_int.in_dll(library, q).value for q in quantity_list}

    def __init__(self, **kwargs) :
        """
        Prepare a MESA based interpolation.
        
        Kewyrd only arguments:
            - mesa_dir: A directory contaning a grid (mass and metallicity)
                        of MESA tracks to base the interpolation on. Must not
                        be specified if interpolator_fname is.
            - interpolator_fname: The filename of a previously saved
                                  interpolator state. Must not be specified
                                  together with mesa_dir.
        Returns: None.
        """

        if 'mesa_dir' in kwargs :
            self.interpolator = library.create_interpolator(
                kwargs['mesa_dir']
            )
        else :
            assert('interpolator_fname' in kwargs)
            self.interpolator = library.load_interpolator(
                kwargs['interpolator_fname']
            )

    def delete() :
        """
        Free the resources allocated at construction.
        """

        library.destroy_interpolator(self.interpolator)

    def save(self, filename) :
        """
        Save the interpolator created to the given file for faster creation.

        Args:
            - filename: The name of the file to use for saving the state.
                        Overwritten if exists.

        Returns: None
        """

        library.save_interpolator(self.interpolator, filename)

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
                                    self.quantity_ids[quantity.upper()],
                                    c_double(mass),
                                    c_double(metallicity))
        )

class Quantity :

    def __init__(self, underlying_quantity) :
        """Wrap the underlying EvolvingStellarQuantity into a callable."""

        self.underlying_quantity = underlying_quantity
        self.min_age = library.quantity_min_age(underlying_quantity)
        self.max_age = library.quantity_max_age(underlying_quantity)

    def delete() :
        """Destroy the underlying quantity."""

        library.destroy_quantity(self.underlying_quantity)

    def __call__(self, age) :
        """
        Evaluate the underlying quantity at the given age(s) (in Gyr).
        
        Args:
            - age: Either a single float or a numpy array of floats giving
                   the ages at which to evaluate the quantity.
        
        Returns: The value(s) of the quantity in the same format as age.
        """

        if type(age) is numpy.ndarray :
            library.evaluate_quantity_array.restype = (
                numpy.ctypeslib.ndpointer(dtype = numpy.float64,
                                          shape = (age.size,),
                                          flags = 'C_CONTIGUOUS')
            )
            return library.evaluate_quantity_array(self.underlying_quantity,
                                                   age,
                                                   age.size)
        else :
            return library.evaluate_quantity(self.underlying_quantity,
                                             c_double(age))

    def deriv(self, age) :
        """
        Return the 0-th, 1-st and 2-nd order derivatives of the quantity.

        Args:
            - age: Either a single float or a numpy array of floats giving
                   the ages at which to evaluate the quantity.

        Returns: A numpy array either 1-D (if age is a single float) or 2-D
                 if age is a numpy array where the outside (or only) index is
                 the derivative order.
        """

        if type(age) is numpy.ndarray :
            library.differentiate_quantity_array.restype = (
                numpy.ctypeslib.ndpointer(dtype = numpy.float64,
                                          shape = (3, age.size),
                                          flags = 'C_CONTIGUOUS')
            )
            return library.differentiate_quantity_array(
                self.underlying_quantity,
                age,
                age.size
            )
        else :
            return library.differentiate_quantity(self.underlying_quantity,
                                                  c_double(age))

if __name__ == '__main__' :
    mesa_dir = '../poet_src/StellarEvolution/MESA'
#    interpolator = MESAInterpolator(mesa_dir = mesa_dir)
#    interpolator.save(mesa_dir + '/saved_interpolator')
#    for quantity_name in MESAInterpolator.quantity_list :
#        quantity = interpolator(quantity_name, 1.0, 0.0)
#        print(quantity_name + '(1.0) = ' + repr(quantity(1.0)))
#        print(quantity_name + '(4.6) = ' + repr(quantity(4.6)))
    loaded_interpolator = MESAInterpolator(
        interpolator_fname = ('saved_interpolator')
    )
    for quantity_name in MESAInterpolator.quantity_list :
        quantity = loaded_interpolator(quantity_name, 1.0, 0.25)
        ages = numpy.exp(
            numpy.linspace(numpy.log(max(quantity.min_age, 1e-5)),
                           numpy.log(quantity.max_age),
                           5)[1:-1]
        )

        print(quantity_name
              +
              '('
              +
              repr(ages)
              +
              ') = '
              +
              repr(quantity(ages)))
