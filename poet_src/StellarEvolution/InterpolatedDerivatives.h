#ifndef __INTERPOLATED_DERIVATIVES_H
#define __INTERPOLATED_DERIVATIVES_H

#include "mass_metallicity_interp.h"
#include "../Core/LogDerivatives.h"
#include "../Core/Functions.h"
#include <vector>

namespace StellarEvolution {

    ///\brief Derivative class for stellar quantities which are interpolated
    ///age, mass and metallicity.
    ///
    ///\ingroup StellarSystem_group
    class InterpolatedDerivatives : public LogDerivatives {
    private:
        double 
            ///The mass to interpolate to in \f$M_\odot\f$.
            __stellar_mass,

            ///The metallicity to interpolate to in Solar metallicities;
            __stellar_metallicity;

        ///The age derivatives for each stellar model.
        std::vector<const FunctionDerivatives *> *__interp_deriv;


        const alglib::real_1d_array 
            ///The masses of the stelar models in \f$M_\odot\f$
            &__interp_masses,

            ///The metallicities of the stellar models in Solar metallicities.
            &__interp_metallicities;

        ///Whether to delete the derivatives it was created with
        bool __delete_derivatives;
    protected:
        ///Returns the deriv_order-th derivative of the quantity
        double calc_deriv(unsigned deriv_order) const;
    public:
        ///\brief Create an object that interpolates derivatives from
        ///evolution tracks.
        ///
        ///The input grid of derivatives may be of quantity of log(quantity)
        ///vs. age or log(age). The returned derivatives are always of
        ///quantity vs age (no log of anything).
        InterpolatedDerivatives(
            ///The stellar mass at which to evaluate the interpolated
            ///derivatives.
            double mass,

            ///The stellar metallicity at which to evaluate the interpolated
            ///derivatives.
            double metallicity,

            ///The derivatives at each grid intersection.
            std::vector<const FunctionDerivatives*> *derivatives,

            ///The masses of the grid intersections.
            const alglib::real_1d_array &interp_masses,

            ///The sorted metallicities of the grid intersections.
            const alglib::real_1d_array &interp_metallicities,

            ///If not NaN, \p derivatives are assumed to calculate
            ///derivatives w.r.t. ln(age), whereas this object always returns
            ///derivatives w.r.t. age.
            double age = NaN,

            ///Are \p derivatives of log(quantity) instead of quantity?
            bool log_quantity = false,

            ///Should \p derivatives be deleted when this object is
            ///destroyed?
            bool delete_derivatives = false
        );

        ///Deletes the interpolation data if so specified on creation.
        ~InterpolatedDerivatives()
        {
            if(__delete_derivatives) {
                for(size_t i = 0; i < __interp_deriv->size(); i++)
                    delete (*__interp_deriv)[i];
                delete __interp_deriv; 
            }
        }
    }; //End of InterpolatedDerivatives class.

} //End of StellarEvolution namespace.

#endif
