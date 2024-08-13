#!/usr/bin/env python3

"""An interface to the POET stellar evolution interpolation utilities."""

from ctypes import\
    cdll,\
    c_int, c_double, c_void_p, c_char_p, c_uint, c_bool,\
    byref
from ctypes.util import find_library
import re
from collections import namedtuple

import numpy

#class naming convention mimicks ctypes naming
#This is just a placeholder, so no public methods.
#pylint: disable=invalid-name
#pylint: disable=too-few-public-methods
class c_interpolator_p(c_void_p):
    """Type corresponding to pointer to intepolator in the POET library."""

class c_quantity_p(c_void_p):
    """Type corresponding to pointer to evolution quantities in POET library."""
#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

def initialize_library():
    """Prepare the stellarEvolution library for use."""

    library_fname = find_library('stellarEvolution')
    if library_fname is None:
        raise OSError('Unable to find POET\'s stellarEvolution library.')
    result = cdll.LoadLibrary(library_fname)

    num_quantities = c_int.in_dll(result, 'NUM_QUANTITIES').value

    result.set_interp_quantity_lower_limit.argtypes = [c_int, c_double]
    result.set_interp_quantity_lower_limit.restype = None

    result.create_interpolator.argtypes = [
        c_char_p,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  shape=(num_quantities,),
                                  flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_int,
                                  ndim=1,
                                  shape=(num_quantities,),
                                  flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_bool,
                                  ndim=1,
                                  shape=(num_quantities,),
                                  flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_bool,
                                  ndim=1,
                                  shape=(num_quantities,),
                                  flags='C_CONTIGUOUS'),
        c_uint
    ]
    result.create_interpolator.restype = c_interpolator_p

    result.destroy_interpolator.argtypes = [
        result.create_interpolator.restype
    ]
    result.destroy_interpolator.restype = None

    result.create_quantity.argtypes = [result.create_interpolator.restype,
                                       c_int,
                                       c_double,
                                       c_double]
    result.create_quantity.restype = c_quantity_p

    result.destroy_quantity.argtypes = [c_quantity_p]
    result.destroy_quantity.restype = None

    result.evaluate_quantity.argtypes = [c_quantity_p, c_double]
    result.evaluate_quantity.restype = c_double

    result.evaluate_quantity_array.argtypes = [
        c_quantity_p,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
    ]

    result.quantity_min_age.restype = c_double

    result.quantity_max_age.restype = c_double

    result.quantity_continuous_range.restype = None

    result.save_interpolator.argtypes = [
        result.create_interpolator.restype,
        c_char_p
    ]
    result.save_interpolator.restype = None

    result.load_interpolator.argtypes = [c_char_p]
    result.load_interpolator.restype = result.create_interpolator.restype

    result.differentiate_quantity.argtypes = [
        c_quantity_p,
        c_double,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]

    result.differentiate_quantity_array.argtypes = [
        c_quantity_p,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS')
    ]

    result.default_smoothing.argtypes = [c_int]
    result.default_smoothing.restype = c_double

    result.default_vs_log_age.argtypes = [c_int]
    result.default_vs_log_age.restype = c_bool

    result.default_nodes.argtypes = [c_int]
    result.default_nodes.restype = c_int

    result.metallicity_from_feh.argtypes = [c_double]
    result.metallicity_from_feh.restype = c_double

    result.feh_from_metallicity.argtypes = [c_double]
    result.feh_from_metallicity.restype = c_double

    result.feh_from_z.argtypes = [c_double]
    result.feh_from_z.restype = c_double

    result.z_from_feh.argtypes = [c_double]
    result.z_from_feh.restype = c_double

    lib_constants = ['Yprimordial', 'Yprotosun', 'Zprotosun']
    result.constants = namedtuple(
        'StellarEvolutionConstants',
        lib_constants
    )(
        *(
            c_double.in_dll(result, const_name).value
            for const_name in lib_constants
        )
    )

    return result

library = initialize_library()

library_track_fname_rex = re.compile(
    'M(?P<MASS>[0-9.E+-]+)_Z(?P<Z>[0-9.E+-]+).csv'
)

def library_track_fname(mass, feh):
    """
    Returns the base name expected by library for a track.

    Args:
        mass:    The mass of the star whose evolution is stored in the track.

        feh:    The [Fe/H] value of the star whose evolution is stored in the
            track.

    Returns:
        str:
            The base filename the stellar evolution library expects to be used
            for the given track.
    """

    return 'M%s_Z%s.csv' % (repr(float(mass)),
                            repr(library.z_from_feh(feh)))

class MESAInterpolator:
    """A class for interpolating among a set of MESA tracks."""

    quantity_list = ['RADIUS', 'ICONV', 'LUM', 'IRAD', 'MRAD', 'RRAD']

    quantity_ids = {q: c_int.in_dll(library, q).value for q in quantity_list}
    quantity_names = {c_int.in_dll(library, q).value: q
                      for q in quantity_list}

    default_smoothing = {q_name: library.default_smoothing(q_id)
                         for q_name, q_id in quantity_ids.items()}

    default_nodes = {q_name: library.default_nodes(q_id)
                     for q_name, q_id in quantity_ids.items()}

    default_vs_log_age = {q_name: library.default_vs_log_age(q_id)
                          for q_name, q_id in quantity_ids.items()}

    default_log_quantity = {q_name: library.default_log_quantity(q_id)
                            for q_name, q_id in quantity_ids.items()}

    @classmethod
    def set_quantity_lower_limit(cls, quantity, limit):
        """
        Set up a min value when interpolating of the given quantity.

        Must be called before creating an interpolator object.
        """

        library.set_interp_quantity_lower_limit(
            cls.quantity_ids[quantity.upper()],
            1e-6
        )


    def __init__(self, **kwargs):
        """
        Prepare a MESA based interpolation.

        Args:
            mesa_dir:    A directory contaning a grid (mass and metallicity) of
                MESA tracks to base the interpolation on. Must not be specified
                if interpolator_fname is.

            smoothing:    A numpy float array of the smoothing arguments to use
                for the interpolation of each quantity. Should be in the order
                defined by quantity_ids.

            nodes:    A numpy integer array of the nodes to use for the
                interpolation of each quantity. Same order as smoothing.

            vs_log_age:    A numpy boolean array indicating whether the
                interpolation for each quantity should be done vs log(age)
                instead of age.

            log_quantity:    A numpy boolean array indicating whether the
                interpolation for each quantity should be of log(quantity)
                instead of quantity.

            interpolator_fname:    The filename of a previously saved
                interpolator state. Must not be specified together with
                mesa_dir. If passed, the smoothing and nodes arguments are
                ignored.

            num_threads:    The number of simultaneous threads to use when
                constructing the interpolation.

        Returns: None.
        """

        if 'mesa_dir' in kwargs:
            self.interpolator = library.create_interpolator(
                kwargs['mesa_dir'].encode('ascii'),
                kwargs['smoothing'],
                kwargs['nodes'],
                kwargs['vs_log_age'],
                kwargs['log_quantity'],
                kwargs['num_threads']
            )
        else:
            assert 'interpolator_fname' in kwargs
            self.filename = kwargs['interpolator_fname']
            self.interpolator = library.load_interpolator(
                kwargs['interpolator_fname'].encode('ascii')
            )

    @classmethod
    def get_create_interpolator_config(cls, **custom_config):
        """
        Return args for creating new interpolator, filling defaults as needed.

        Args:
            custom_config:    Configuration for how to interpolate the POET
                relevant quantities. May have an entry for everything in
                `MESAInterpolator.quantity_list`, omitted quantities use the
                defaults from MESAInterpolator. For each qunatity user may
                specify: ``smoothing``, ``nodes``, ``vs_log_age``,
                ``log_quantity``.

        Returns:
            dict:
                Entries for ``smoothing``, ``nodes``, ``vs_log_age``, and
                ``log_quantity`` to pass to __init__ to create an interpolator.
        """

        num_quantities = len(cls.quantity_list)

        kwargs = dict(
            smoothing=numpy.empty(num_quantities, dtype=c_double),
            nodes=numpy.empty(num_quantities, dtype=c_int),
            vs_log_age=numpy.empty(num_quantities, dtype=c_bool),
            log_quantity=numpy.empty(num_quantities, dtype=c_bool)
        )

        for q_name, q_index in cls.quantity_ids.items():
            q_config = custom_config.get(q_name, dict())

            for param in kwargs:
                kwargs[param][q_index] = q_config.get(
                    param,
                    getattr(cls, 'default_' + param)[q_name]
                )
        return kwargs


    def delete(self):
        """Free the resources allocated at construction."""

        library.destroy_interpolator(self.interpolator)

    def save(self, filename):
        """
        Save the interpolator created to the given file for faster creation.

        Args:
            filename:    The name of the file to use for saving the state.
                Overwritten if exists.

        Returns:
            None
        """

        self.filename = filename
        library.save_interpolator(self.interpolator,
                                  filename.encode('ascii'))

    def __call__(self, quantity, mass, feh):
        """
        Return a stellar quantity interpolated to the given mass and [Fe/H].

        Args:
            quantity:    A string identifying the quantity to interpolate. The
                following values are allowed: 'radius', 'iconv', 'lum',
                'irad', 'mrad', 'rrad'. This is a case insensitive argument.

            mass:    The mass of the star for which this quantity should be
                defined in solar masses.

            feh:    The [Fe/H] of the star for which this  quantity  should
                be defined.

        Returns:
            Quantity:
                callable with an age parameter evaluating to the quantity at the
                given age.
        """

        return Quantity(
            library.create_quantity(self.interpolator,
                                    self.quantity_ids[quantity.upper()],
                                    c_double(mass),
                                    c_double(feh))
        )

class Quantity:
    """Callable that evaluates to the value of the quantity at a given age."""

    def __init__(self, underlying_quantity):
        """Wrap the underlying EvolvingStellarQuantity into a callable."""

        self.underlying_quantity = underlying_quantity
        self.min_age = library.quantity_min_age(underlying_quantity)
        self.max_age = library.quantity_max_age(underlying_quantity)

    def delete(self):
        """Destroy the underlying quantity."""

        library.destroy_quantity(self.underlying_quantity)

    def __call__(self, age):
        """
        Evaluate the underlying quantity at the given age(s) (in Gyr).

        Args:
            age:    Either a single float or a numpy array of floats giving
               the ages at which to evaluate the quantity.

        Returns:
            type(age):
                The value(s) of the quantity in the same format as age.
        """

        if isinstance(age, numpy.ndarray):
            assert (age >= self.min_age).all()
            assert (age <= self.max_age).all()
            result = numpy.empty(dtype=c_double,
                                 shape=(age.size,),
                                 order='C')
            library.evaluate_quantity_array(self.underlying_quantity,
                                            age,
                                            age.size,
                                            result)
            return result

        assert age > self.min_age
        assert age < self.max_age
        return library.evaluate_quantity(self.underlying_quantity,
                                         c_double(age))

    def deriv(self, age):
        """
        Return the 0-th, 1-st and 2-nd order derivatives of the quantity.

        Args:
            age(float or numpy array):    Either a single float or a numpy array
                of floats giving the ages at which to evaluate the quantity.

        Returns:
            numpy array:
                Either 1-D (if age is a single float) or 2-D array if age is a
                numpy array where the outside (or only) index is the derivative
                order.
        """

        if isinstance(age, numpy.ndarray):
            result = numpy.empty(dtype=c_double,
                                 shape=(3, age.size),
                                 order='C')
            library.differentiate_quantity_array(
                self.underlying_quantity,
                age,
                age.size,
                result.reshape(3 * age.size)
            )
        else:
            result = numpy.empty(dtype=c_double,
                                 shape=(3,),
                                 order='C')
            library.differentiate_quantity(self.underlying_quantity,
                                           c_double(age),
                                           result)
        return result

    def continuous_range(self, age):
        """
        Return the range around age over which the quantity is continuous.

        Args:
            - age: The age around which the continuous region is required.

        Returns: A 2-tuple of the minimum and maximum ages surrounding age
                 over thich the quantity is guaranteed continuous
        """

        min_age, max_age = c_double(), c_double()
        library.quantity_continuous_range(self.underlying_quantity,
                                          c_double(age),
                                          byref(min_age),
                                          byref(max_age))
        return min_age.value, max_age.value

def example():
    """Example of the usage of the interface."""

#    mesa_dir = '../poet_src/StellarEvolution/MESA'
#    interpolator = MESAInterpolator(mesa_dir = mesa_dir)
#    interpolator.save(mesa_dir + '/saved_interpolator')
#    for quantity_name in MESAInterpolator.quantity_list :
#        quantity = interpolator(quantity_name, 1.0, 0.0)
#        print(quantity_name + '(1.0) = ' + repr(quantity(1.0)))
#        print(quantity_name + '(4.6) = ' + repr(quantity(4.6)))
    loaded_interpolator = MESAInterpolator(
        interpolator_fname=('saved_interpolator')
    )
    for quantity_name in MESAInterpolator.quantity_list:
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

if __name__ == '__main__':
    example()
