/**\file
 *
 * \brief Declare a class for a stellar evolution quantity which is the sum
 * of two other quantities.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __SUM_QUANTITY_H
#define __SUM_QUANTITY_H

#include "../Core/SharedLibraryExportMacros.h"
#include "EvolvingStellarQuantity.h"

namespace StellarEvolution {

    ///\brief A clas for stellar quantities that are the sum of two other
    ///quantities.
    ///
    ///\ingroup StellarSystem_group
    class LIB_LOCAL SumQuantity : public EvolvingStellarQuantity {
    private:
        ///This quantity will be q1+q2
        EvolvingStellarQuantity *q1, *q2;

        ///Whether to destroy the input quantities on destruction
        bool destroy_qs;
    public:
        ///Create a quantity that is (*quantity1)-(*quantity2)
        SumQuantity(EvolvingStellarQuantity *quantity1,
                    EvolvingStellarQuantity *quantity2,
                    bool delete_inputs=false)
            : q1(quantity1), q2(quantity2), destroy_qs(delete_inputs) {}

        ///See EvolvingStellarQuantity::select_interpolation_region.
        virtual void select_interpolation_region(double age) const
        {
            q1->select_interpolation_region(age);
            q2->select_interpolation_region(age);
        }

        ///Return the value the quantity takes at the given age.
        double operator()(double age) const
        {return (*q1)(age)+(*q2)(age);}

        ///Return the age derivative of the quantity at the given age.
        const FunctionDerivatives *deriv(double age) const
        {return new SumDerivatives(q1->deriv(age), q2->deriv(age), true);}

        ///The largest age for which the quantity can be interpolated
        double range_high() const
        {return std::min(q1->range_high(), q2->range_high());}

        ///The smallest age for which the quantity can be interpolated.
        double range_low() const
        {return std::max(q1->range_low(), q2->range_low());}

        ///See EvolvingStellarQuantity::next_discontinuity.
        double next_discontinuity() const
        {return std::min(q1->next_discontinuity(), q2->next_discontinuity());}

        ///See EvolvingStellarQuantity::enable_next_interpolation_region.
        void enable_next_interpolation_region() const
        {
            if(q1->next_discontinuity() == q2->next_discontinuity()) {
                q1->enable_next_interpolation_region();
                q2->enable_next_interpolation_region();
            } else if(q1->next_discontinuity() < q2->next_discontinuity())
                q1->enable_next_interpolation_region();
            else 
                q2->enable_next_interpolation_region();
        }

        ///An iterator over the ages where the quantity takes the given y value.
        InterpSolutionIterator crossings(double =0) const
        {
            throw Core::Error::Runtime(
                "Called EvolvingStellarQuantity::crossings, which are ill "
                "defined."
            );
        }

        ///Clean up.
        ~SumQuantity()
        {if(destroy_qs) {delete q1; delete q2;}}
    }; //End SumQuantity declaration.


} //End StellarEvolution namespace.

#endif
