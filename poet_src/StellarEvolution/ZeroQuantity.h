/**\file
 * 
 * \brief Declare a class for an identically zero stellar evolution quantity.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __ZERO_QUANTITY_H
#define __ZERO_QUANTITY_H

#include "../Core/SharedLibraryExportMacros.h"
#include "EvolvingStellarQuantity.h"

namespace StellarEvolution {

    class LIB_LOCAL ZeroQuantity : public EvolvingStellarQuantity {
    public:
        ///Do nothing. See EvolvingStellarQuantity::select_interpolation_region.
        void select_interpolation_region(double age) {}

        ///Return the value the quantity takes at the given age.
        double operator()(double) const {return 0;}

        ///Return the age derivative of the quantity at the given age.
        const FunctionDerivatives *deriv(double) const
        {return new Core::ZeroDerivatives;}

        ///The largest age for which the quantity can be interpolated
        double range_high() const {return Inf;}

        ///The smallest age for which the quantity can be interpolated.
        double range_low() const {return -Inf;}

        ///No discontinuities. See EvolvingStellarQuantity::next_discontinuity.
        double next_discontinuity() const {return Inf;}

        ///\brief Do nothing.
        ///See EvolvingStellarQuantity::enable_next_interpolation_region.
        double enable_next_interpolation_region();

        ///An iterator over the ages where the quantity takes the given y value.
        InterpSolutionIterator crossings(double =0) const
        {throw Error::Runtime("Called ZeroQuantity::crossings, "
                              "which are ill defined.");}
    }; //End ZeroQuantity declaration.

} //End StellarEvolution namespace.

#endif
