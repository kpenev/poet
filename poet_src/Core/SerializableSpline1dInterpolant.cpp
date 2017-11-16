#include "SerializableSpline1dInterpolant.h"

namespace Core {
    std::string alglib_dtype_name(alglib_impl::ae_datatype dtype)
    {
        switch(dtype) {
            case alglib_impl::DT_BOOL: return "bool";
            case alglib_impl::DT_INT: return "int";
            case alglib_impl::DT_REAL: return "real";
            case alglib_impl::DT_COMPLEX : return "complex";
            default : return "unknown: " + std::to_string(dtype);
        }
    }

    std::ostream &operator<<(std::ostream &os,
                             const SerializableSpline1dInterpolant &interp)
    {
        alglib_impl::spline1dinterpolant *interp_impl = 
            const_cast<alglib_impl::spline1dinterpolant *>(interp.c_ptr());

        os << "\tperiodic: " << interp_impl->periodic << std::endl
           << "\tn: " << interp_impl->n << std::endl
           << "\tk: " << interp_impl->k << std::endl;

        os << "\tc.cnt: " << interp_impl->c.cnt << std::endl
           << "\tc.datatype: "
           << alglib_dtype_name(interp_impl->c.datatype)
           << std::endl
           << "\tc.is_attached: " << interp_impl->c.is_attached << std::endl
           << "\tc.ptr.p_double (@ " << interp_impl->c.ptr.p_double << "):";
        for(int i = 0; i < 4; ++i)
            os << " " << interp_impl->c.ptr.p_double[i];
        os << std::endl;

        os << "\tx.cnt: " << interp_impl->x.cnt << std::endl
           << "\tx.datatype: "
           << alglib_dtype_name(interp_impl->x.datatype)
           << std::endl
           << "\tx.is_attached: " << interp_impl->x.is_attached << std::endl
           << "\tx.ptr.p_double (@ " << interp_impl->x.ptr.p_double << "):";
        for(int i = 0; i < 2; ++i)
            os << " " << interp_impl->x.ptr.p_double[i];
        os << std::endl;
        return os;
    }

    std::ostream &operator<<(std::ostream &os,
                             const alglib::real_1d_array &array)
    {
        os << "len: " << array.length() 
           << ", @: " << array.getcontent()
           << ":";
        for(int i = 0; i<3; ++i) os << " " << array[i];
        return os;
    }

} //End Core namespace.
