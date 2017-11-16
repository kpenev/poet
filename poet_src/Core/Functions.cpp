/**\file
 * 
 * \brief The definition of some of the methods of the various function
 * classes.
 *
 * \ingroup Utilities_group
 */

#include "Functions.h"
#include "Common.h"
#include "InterpSolutionIterator.h"
#include <iostream>

#ifndef NO_SERIALIZE
    BOOST_CLASS_EXPORT_IMPLEMENT(Core::ZeroFunction)
#endif

namespace Core {

    CubicSplineDerivatives::CubicSplineDerivatives(double func_value, 
                                                   double first_deriv,
                                                   double second_deriv) :
        zeroth(func_value), first(first_deriv), second(second_deriv)
    {}

    double CubicSplineDerivatives::order(unsigned deriv_order) const
    {
        if(deriv_order==0) return zeroth;
        else if(deriv_order==1) return first;
        else if(deriv_order==2) return second;
        else if(deriv_order==3) return NaN;
        else return 0.0;
    }

    InterpSolutionIterator ZeroFunction::crossings(double) const
    {
        throw Error::Runtime("Called ZeroQuantity::crossings, "
                             "which are ill defined.");
    }


} //End Core namespace.
