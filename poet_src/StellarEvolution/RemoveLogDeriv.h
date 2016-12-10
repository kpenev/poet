/**\file
 *
 * \brief Declares a class to switch from logarithmic to linear derivative.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __REMOVE_LOG_DERIV_H
#define __REMOVE_LOG_DERIV_H

#include "LogArgDerivatives.h"

namespace StellarEvolution {

    ///\brief Makes a derivative with respect to linear argument from a
    ///derivative with respect to log(argument).
    class RemoveLogDeriv : public LogArgDerivatives {
    private:
        ///The original logarithmic derivative 
        const FunctionDerivatives *__log_deriv;

        ///Whether to delete the underlying log-derivative on destruction.
        bool __delete_deriv;
    protected:
        ///Returns the deriv_order-th derivative of the quantity
        double calc_deriv(unsigned deriv_order) const
        {return __log_deriv->order(deriv_order);}
    public:
        ///Create a linear derivative from a log one.
        RemoveLogDeriv(double age, const FunctionDerivatives *log_deriv,
                       bool delete_deriv) :
            LogArgDerivatives(age), __log_deriv(log_deriv), 
            __delete_deriv(delete_deriv) {}

        ///\brief Deletes the input logarithmic derivative if so specified on 
        ///creation.
        ~RemoveLogDeriv()
        {if(__delete_deriv) delete __log_deriv;}
    }; //End RemoveLogDeriv class declaration.

} //End StellarEvolution namespace.

#endif
