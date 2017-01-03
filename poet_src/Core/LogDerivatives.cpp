/**\file
 *
 * \brief Defines some of the methods of the EvolvingStellarQuantity class
 * used for interpolating among stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#include "LogDerivatives.h"

namespace StellarEvolution {

    double LogDerivatives::transform_log_x_deriv(unsigned order) const
    {
        switch(order) {
            case 0:
                return __underlying_deriv_values[0];
            case 1:
                return __underlying_deriv_values[1] / __x;
            case 2:
                return (
                    __underlying_deriv_values[2]
                    -
                    __underlying_deriv_values[1]
                )/(__x * __x);
            default: 
                throw Core::Error::BadFunctionArguments(
                    "Transforming log-derivatives of order higher than 2 is"
                    " not implemented."
                );
        }
    }

    double LogDerivatives::transform_log_y_deriv(
        double uncorrected_derivative,
        unsigned order
    ) const
    {
        switch(order) {
            case 0: 
                return std::exp(uncorrected_derivative);
            case 1: 
                return __deriv_values[0] * uncorrected_derivative;
            case 2:
                return (
                    __deriv_values[0] * uncorrected_derivative
                    +
                    __deriv_values[1] * __deriv_values[1] / __deriv_values[0]
                );
        }
    }

    double LogDerivatives::order(unsigned deriv_order) const
    {
        if(__deriv_values.size() <= deriv_order) {
            for(
                unsigned order = __deriv_values.size();
                order < deriv_order; 
                ++order
            ) {
                __underlying_deriv_values.push_back(calc_deriv(order));
                double corrected_deriv_value = (
                    __log_x
                    ? transform_log_x_deriv(order)
                    : calc_deriv(order)
                );
                if(__log_y)
                    corrected_deriv_value = transform_log_y_deriv(
                        corrected_deriv_value,
                        order
                    );
                __deriv_values.push_back(corrected_deriv_value);
            }
        }
        return __deriv_values[deriv_order];
    }

} //End of StellarEvolution namespace.
