#include "InterpolatingFunctionALGLIB.h"

#ifndef NO_SERIALIZE
    BOOST_CLASS_EXPORT_IMPLEMENT(Core::InterpolatingFunctionALGLIB)
#endif

namespace Core {

    InterpolatingFunctionALGLIB::InterpolatingFunctionALGLIB(
        const double *x,
        const double *y,
        size_t num_points,
        const double *yprime,
        double smoothing,
        int degrees_of_freedom
    )
    {
        if(degrees_of_freedom<0) degrees_of_freedom=std::min(
            -degrees_of_freedom, static_cast<int>(3 * num_points));

        alglib::real_1d_array alglib_x;
        alglib_x.setcontent(num_points, x);
        alglib::real_1d_array alglib_y;
        alglib_y.setcontent(num_points, y);

        alglib::real_1d_array alglib_yprime;
        if(yprime) alglib_yprime.setcontent(num_points, yprime);

        auto x_range = std::minmax_element(x, x + num_points);
        __min_x = *x_range.first;
        __max_x = *x_range.second;

        if(std::isnan(smoothing)) {
            if(yprime)
                alglib::spline1dbuildhermite(alglib_x,
                                             alglib_y,
                                             alglib_yprime,
                                             __spline);
            else 
                alglib::spline1dbuildcubic(alglib_x,
                                           alglib_y,
                                           __spline);
        }
        else {
            if(yprime)
                throw Error::BadFunctionArguments(
                    "Smoothing not supported when derivatives "
                    "are specified in InterpolatingFunctionALGLIB "
                    "constructor.");
            alglib::ae_int_t fit_info;
            alglib::spline1dfitreport fit_report;
            alglib::spline1dfitpenalized(alglib_x,
                                         alglib_y, 
                                         degrees_of_freedom,
                                         smoothing,
                                         fit_info, 
                                         __spline,
                                         fit_report);
            if(fit_info<=0) throw Error::ALGLIB("Spline fit failed.");
        }
    }

    InterpSolutionIterator InterpolatingFunctionALGLIB::crossings(
        double y
    )  const
    {
        return InterpSolutionIterator(__spline, y);
    }



} //End Core namespace.
