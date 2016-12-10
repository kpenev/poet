/**\file
 *
 * \brief Declares a class for differentiating functions of log(arg) w.r.t.
 * arg.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __LOG_ARG_DERIVATIVES_H
#define __LOG_ARG_DERIVATIVES_H

#include "../Core/Functions.h"
#include "../Core/Error.h"
#include <vector>

namespace StellarEvolution {

    using Core::FunctionDerivatives;
    using Core::NaN;

    ///\brief A class that calculates derivatives with respect to an argument for
    ///functions of the log(argument).
    ///
    ///\ingroup StellarSystem_group
    class LogArgDerivatives : public FunctionDerivatives {
    private:
        ///The value of the argument at which derivatives are calculated.
        double x;

        ///The value (zeroth derivative).
        mutable double value;

        ///\brief Was the interpolation done against the logarithm of the
        ///argument and hence the derivative needs to be corrcted.
        bool correct_log_arg;

        ///\brief The currently computed derivatives.
        ///
        ///With respect to log(arg) if correct_log_arg is true, or with respect
        ///to arg if not.
        mutable std::vector<double> underlying_deriv_values,

                ///Previously calculated values of the derivatives.
                ///
                ///These are reused if requested multiple times.
                deriv_values;

        ///\brief Actually corrects for differentiating w.r.t. log(arg) instead
        ///of arg.
        ///
        ///Uses the pre-computed array of the derivatives up to the given order
        ///in underlying_deriv_values with respect to ln(x) to return the
        ///order-th derivative with respect to x.
        double transform_log_arg_deriv(unsigned order) const;
    protected:
        ///\brief Should be overwritten to calculate the derivatives with respect
        ///to either arg or log(arg) as specified on construction.
        virtual double calc_deriv(unsigned deriv_order) const =0;
    public:
        ///\brief Create a derivative for functions of possibly log(arg).
        ///
        ///The created object corrects for the fact that the underlying
        ///derivative (defined by the calc_deriv method) is with respect to the
        ///logarithm of the argument if arg_val is not NaN. No correction if it
        ///is Nan.
        LogArgDerivatives(double arg_val=NaN) :
            x(arg_val), value(NaN), correct_log_arg(!std::isnan(arg_val)) {}

        ///Returns the deriv_order-th derivative of the quantity
        double order(unsigned deriv_order=1) const;
    }; //End of LogArgDerivatives class

} //End of StellarEvolution namespace.

#endif
