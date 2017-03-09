import numpy

def ndpointer_or_null(*args, **kwargs) :
    """
    Allow None (->NULL) to be passed for c-style array function arguments.

    Modified from:
    http://stackoverflow.com/questions/32120178/how-can-i-pass-null-to-an-external-library-using-ctypes-with-an-argument-decla
    """

    base = numpy.ctypeslib.ndpointer(*args, **kwargs)

    def from_param(cls, obj) :
        """Construct numpy.ndpointer from the given object."""

        if obj is None :
            return obj
        else :
            return base.from_param(obj)

    return type(base.__name__, (base,), 
                {'from_param': classmethod(from_param)})


