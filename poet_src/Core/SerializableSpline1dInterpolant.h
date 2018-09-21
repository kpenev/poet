#ifndef __SERIALIZABLE_SPLINE_1D_INTERPOLANT_H
#define __SERIALIZABLE_SPLINE_1D_INTERPOLANT_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../third_party_libs/alglib/alglib/src/interpolation.h"

#include <boost/serialization/base_object.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>

namespace Core {

    ///A serializable (using boost serialization) alglib 1D interpolant.
    class LIB_LOCAL SerializableSpline1dInterpolant : 
        public alglib::spline1dinterpolant
    {
    private:
        friend class boost::serialization::access;

        ///Serialize (see boost serialization library) the interpolation.
        ///
        ///Serializes everything EXCEPT p_struct->x.data and
        ///p_struct->y.data, because I have no idea what they are. Also 
        ///assumes that the spline data are all stored as doubles.
        template<class Archive> void serialize(
            ///The archive to serialize to.
            Archive & ar,

            ///Version number. Ignored!
            const unsigned int
        );
    public:
        ///\brief Human readable output for interpolants (data is currently 
        ///not output).
        friend std::ostream &operator<<(
            ///The stream to output to.
            std::ostream &os,

            ///The interpolant to output.
            const SerializableSpline1dInterpolant &interp
        );
    }; //End SerializableSpline1dInterpolant class.

    ///Human readable output for ALGLIB arrays (no data is output).
    std::ostream &operator<<(
        ///The stream to output to.
        std::ostream &os,

        ///The vector to output.
        const alglib::real_1d_array &array);

    template<class Archive> void SerializableSpline1dInterpolant::serialize(
        Archive & ar,
        const unsigned int
    )
    {
        ar & p_struct->periodic;
        ar & p_struct->n;
        ar & p_struct->k;
        ar & p_struct->c.cnt;
        ar & p_struct->c.datatype;
        ar & p_struct->x.cnt;
        ar & p_struct->x.datatype;
        if (Archive::is_loading::value) {
            using namespace alglib_impl;
            ae_state state;
            ae_vector_init(&p_struct->c,
                           (p_struct->n - 1) * 4 + 2,
                           p_struct->c.datatype,
                           &state);
            ae_vector_init(&p_struct->x,
                           p_struct->n,
                           p_struct->x.datatype,
                           &state);
        }
        for (
            int node_index=0;
            node_index < p_struct->n - 1;
            node_index++
        ) {
            double *coef=(p_struct->c.ptr.p_double + 4 * node_index);
            ar & coef[0] & coef[1] & coef[2] & coef[3];
            ar & p_struct->x.ptr.p_double[node_index];
        }
        ar & p_struct->x.ptr.p_double[p_struct->n - 1];
    } //End SerializableSpline1dInterpolant::serialize definition.

} //End Core namespace.

#endif
