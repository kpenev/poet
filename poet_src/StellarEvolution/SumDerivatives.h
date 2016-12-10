/**\file
 *
 * \brief Declares & defines a class for the derivative of a quantity which
 * is the sum of two other quantities.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __SUM_DERIVATIVES_H
#define __SUM_DERIVATIVES_H

#include "../../Core/Functions.h"

namespace StellarEvolution {

    ///\brief Derivative class for a quantity that is the sum of two other
    ///quantities.
    ///
    ///\ingroup StellarSystem_group
    class SumDerivatives : public FunctionDerivatives {
    private:
        ///The derivatives of the first quantity in the sum.
        const FunctionDerivatives *q1_deriv,

              ///The derivatives of the second quantity in the sum.  
              *q2_deriv;

        ///\brief Whether to delete the input derivative when the object is 
        ///destroyed.
        bool destroy_derivs;
    public:
        ///Create a derivative object for a sum of two quantities: q1+q2.
        SumDerivatives(
            ///Pointer to the derivative of the first quantity (q1).
            const FunctionDerivatives *derivative1,

            ///Pointer to the derivative of the second quantity (q2).
            const FunctionDerivatives *derivative2,

            ///Delete the input derivatives on destruction?
            bool delete_inputs=true)
            : q1_deriv(derivative1), q2_deriv(derivative2),
            destroy_derivs(delete_inputs) {}

        ///The deriv_order-th derivative.
        double order(unsigned deriv_order=1) const
        {return q1_deriv->order(deriv_order)+q2_deriv->order(deriv_order);}

        ///Clean up.
        ~SumDerivatives()
        {if(destroy_derivs) {delete q1_deriv; delete q2_deriv;}}
    }; //End SumDerivatives class declaration.

} //End StellarEvolution namespace.

#endif
