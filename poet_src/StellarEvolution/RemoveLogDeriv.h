/**\file
 *
 * \brief Declares a class to switch from logarithmic to linear derivative.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __REMOVE_LOG_DERIV_H
#define __REMOVE_LOG_DERIV_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/LogDerivatives.h"

namespace StellarEvolution {

    ///\brief Return dy/dx given dy/dln(x), dln(y)/dx or dln(y)/dln(x).
    class LIB_LOCAL RemoveLogDeriv : public LogDerivatives {
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
        RemoveLogDeriv(
            double age, 
            bool log_quantity,
            const FunctionDerivatives *log_deriv,
            bool delete_deriv
        ) :
            LogDerivatives(age, log_quantity),
            __log_deriv(log_deriv), 
            __delete_deriv(delete_deriv)
        {}

        ///\brief Deletes the input logarithmic derivative if so specified on
        ///creation.
        ~RemoveLogDeriv()
        {if(__delete_deriv) delete __log_deriv;}
    }; //End RemoveLogDeriv class declaration.

} //End StellarEvolution namespace.

#endif
