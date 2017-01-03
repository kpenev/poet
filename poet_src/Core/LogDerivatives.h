/**\file
 *
 * \brief Declares a class for differentiating functions of log(arg) w.r.t.
 * arg.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __LOG_ARG_DERIVATIVES_H
#define __LOG_ARG_DERIVATIVES_H

#include "Functions.h"
#include "Error.h"
#include <vector>

namespace StellarEvolution {

    using Core::FunctionDerivatives;
    using Core::NaN;

    ///\brief Calculate dy/dx given dy/dx, dy/dln(x), dln(y)/dx or 
    ///dln(y)/dln(x).
    ///
    ///\ingroup StellarSystem_group
    class LogDerivatives : public FunctionDerivatives {
    private:
        ///The value of the argument at which derivatives are calculated.
        double __x;

        bool 
            ///Is the underlying derivative w.r.t. log(argument)?
            __log_x,
            
            ///Is the underlying derivative of log(quantity)
            __log_y;

        mutable std::vector<double> 
            ///Cache for previously computed underlying derivative values.
            __underlying_deriv_values,

            ///Cache for previously computed corrected derivative values.
            __deriv_values;

        ///\brief Correct for differentiating w.r.t. log(arg) instead of arg.
        ///
        ///Uses the pre-computed array of the derivatives up to the given
        ///order in __underlying_deriv_values with respect to ln(x) to return
        ///the order-th derivative with respect to x.
        double transform_log_x_deriv(unsigned order) const;

        ///\brief Correct for differentiating of log(quantity) instead of
        ///quantity.
        ///
        ///Uses all lower order pre-computed derivatives in __deriv_values.
        double transform_log_y_deriv(
            ///The derivative of log(quantity) w.r.t. argument of the same
            ///order. Must already be corrected for log(argument) if
            ///necessary.
            double uncorrected_derivative,

            ///The order of the derivative being calculated.
            unsigned order
        ) const;

    protected:
        ///\brief Should be overwritten to calculate the derivatives with respect
        ///to either arg or log(arg) as specified on construction.
        virtual double calc_deriv(unsigned deriv_order) const =0;
    public:
        ///\brief Create a derivative for possibly log(functions) of possibly
        ///log(arg).
        ///
        ///The created object corrects for the fact that the underlying
        ///derivative (defined by the calc_deriv method) may be of the
        ///logarithm of a quantity and/or with respect to the
        ///logarithm of the argument.
        LogDerivatives(
            ///If not NaN, the underlying derivatives are assumed to be
            ///w.r.t. log(x) and should be converted to derivatives w.r.t. x.
            double x = NaN,

            ///Should the underlying derivatives of log(y) be converted to
            ///derivatives of y. If false, the underlying derivatives are
            ///assumed to be of y itself.
            bool log_y = false
        ) : __x(x), __log_x(!std::isnan(x)), __log_y(log_y) {}

        ///Returns the deriv_order-th derivative of the quantity
        double order(unsigned deriv_order=1) const;
    }; //End of LogDerivatives class

} //End of StellarEvolution namespace.

#endif
