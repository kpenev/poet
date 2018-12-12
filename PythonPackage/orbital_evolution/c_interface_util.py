"""Utilities for creating the interface to the shared C/C++ libraries."""

import numpy

def ndpointer_or_null(*args, **kwargs):
    """
    Allow None (->NULL) to be passed for c-style array function arguments.

    Modified from:
    http://stackoverflow.com/questions/32120178/how-can-i-pass-null-to-an-external-library-using-ctypes-with-an-argument-decla
    """

    base = numpy.ctypeslib.ndpointer(*args, **kwargs)

    #Call signature dictated by numpy.ctypeslib
    #pylint: disable=unused-argument
    def from_param(cls, obj):
        """Construct numpy.ndpointer from the given object."""

        if obj is None:
            return obj

        return base.from_param(obj)
    #pylint: enable=unused-argument

    return type(base.__name__, (base,),
                {'from_param': classmethod(from_param)})
