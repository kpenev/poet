#ifndef __INTERPOLATING_FUNCTION_ALGLIB_H
#define __INTERPOLATING_FUNCTION_ALGLIB_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Functions.h"
#include "SerializableSpline1dInterpolant.h"
#include "InterpSolutionIterator.h"

namespace Core {

    ///Function which interpolates, with possible smoothing, between points.
    class LIB_LOCAL InterpolatingFunctionALGLIB :
        public OneArgumentDiffFunction {
    private:
#ifndef NO_SERIALIZE
        ///Needed for serialization to work.
        friend class boost::serialization::access;

        ///Serialize this function.
        //The second parameter should eventually be version
        template<class Archive>
            void serialize(Archive & ar, const unsigned int) {
                ar & boost::serialization::base_object<OneArgumentDiffFunction>(
                    *this
                );
                ar & __spline;
                ar & __min_x & __max_x;
            }
#endif

        ///\brief The interpolating function information necessary to
        ///evaluate at any given point.
        SerializableSpline1dInterpolant __spline;

        ///The smallest abscissa covered by the spline points.
        double __min_x,

               ///The largest abscissa covered by the spline points.
               __max_x;
    public:
        ///Needed by the Boost serializer.
        InterpolatingFunctionALGLIB() {};

        ///\brief Constuct an interpolating function.
        ///
        ///Based on the ALGLIB (smoothing) cubic spline interpolation.
        InterpolatingFunctionALGLIB(
            ///The abscissas of tabulated points to fit.
            const std::valarray<double> &x,

            ///The ordinates of tabulated points to fit.
            const std::valarray<double> &y, 

            ///The values of the derivatives to impose on the
            ///nodes.
            const std::valarray<double> &yprime=
            std::valarray<double>(),

            ///How much smoothing to apply. Omit for no
            ///smoothnig, i.e. the interpolating curve passes
            ///through all the points).
            double smoothing=NaN,

            ///How many degrees of freedom to use for smoothing
            ///interpolation. Ignored for non-smoothing
            ///interpolation. If omitted it is set to 3 times the
            ///number of points being fitted.
            int degrees_of_freedom=-1
        )
        : InterpolatingFunctionALGLIB(
            &(x[0]),
            &(y[0]),
            x.size(),
            (yprime.size() ? &(yprime[0]) : NULL),
            smoothing,
            degrees_of_freedom
        )
        {
            assert(x.size() == y.size());
            assert(yprime.size() == 0 || yprime.size() == x.size());
        }

        ///\brief Constuct an interpolating function.
        ///
        ///Based on the ALGLIB (smoothing) cubic spline interpolation.
        InterpolatingFunctionALGLIB(
            ///The abscissas of tabulated points to fit.
            const double *x,

            ///The ordinates of tabulated points to fit.
            const double *y, 

            ///The number of input points.
            size_t num_points,

            ///The values of the derivatives to impose on the
            ///nodes.
            const double *yprime = NULL,

            ///How much smoothing to apply. Omit for no
            ///smoothnig, i.e. the interpolating curve passes
            ///through all the points).
            double smoothing=NaN,

            ///How many degrees of freedom to use for smoothing
            ///interpolation. Ignored for non-smoothing
            ///interpolation. If omitted it is set to 3 times the
            ///number of points being fitted.
            int degrees_of_freedom=-1);

        ///\brief Returns the value of the interpolating function at the
        ///given abscissa.
        inline double operator()(double x) const
        {return alglib::spline1dcalc(__spline, x);}

        ///\brief The derivatives of the interpolating function at the given
        ///abscissa.
        ///
        ///Returns a newly allocated structure, which must be destroyed
        ///when no longer needed.
        inline const CubicSplineDerivatives *deriv(double x) const
        {
            double v, dv, d2v;
            alglib::spline1ddiff(__spline, x, v, dv, d2v); 
            return new CubicSplineDerivatives(v, dv, d2v);
        }

        ///The maximum abscissa at which the function is defined.
        inline double range_high() const {return __max_x;}

        ///The minimum abscissa at which the function is defined.
        inline double range_low() const {return __min_x;}

        ///\brief Iterator over the abscissas where the function takes the
        ///given y value.
        InterpSolutionIterator crossings(double y=0) const;
    }; //End InterpolatingFunctionALGLIB class
    //BOOST_CLASS_EXPORT_GUID(OneArgumentDiffFunction, 
    //                        "OneArgumentDiffFunction")
 
} //End Core namespace

#ifndef NO_SERIALIZE
    BOOST_CLASS_EXPORT_KEY(Core::InterpolatingFunctionALGLIB)
    BOOST_CLASS_EXPORT_KEY(Core::ZeroFunction)
#endif

        //BOOST_CLASS_EXPORT_GUID(Core::InterpolatingFunctionALGLIB,
        //	"InterpolatingFunctionALGLIB")

#endif
