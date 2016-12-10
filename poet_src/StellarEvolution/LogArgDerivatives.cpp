/**\file
 *
 * \brief Defines some of the methods of the EvolvingStellarQuantity class
 * used for interpolating among stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#include "LogArgDerivatives.h"

namespace StellarEvolution {

    double LogArgDerivatives::transform_log_arg_deriv(unsigned order) const
    {
        if(order==1) return underlying_deriv_values[0]/x;
        else if(order==2)
            return (underlying_deriv_values[1]-underlying_deriv_values[0])/(x*x);
        else throw Core::Error::BadFunctionArguments(
            "Transforming log-derivatives of order higher than 2 is not "
            "implemented."
        );
    }

    double LogArgDerivatives::order(unsigned deriv_order) const
    {
        if(deriv_order==0) {
            if(std::isnan(value)) value=calc_deriv(0);
            return value;
        }
        if(deriv_values.size()<deriv_order)
            for(unsigned i=deriv_values.size(); i<deriv_order; ++i) {
                if(correct_log_arg) {
                    underlying_deriv_values.push_back(calc_deriv(i+1));
                    deriv_values.push_back(transform_log_arg_deriv(i+1));
                } else deriv_values.push_back(calc_deriv(i+1));
            }
        return deriv_values[deriv_order-1];
    }

} //End of StellarEvolution namespace.
