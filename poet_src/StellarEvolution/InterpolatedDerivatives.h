#ifndef __INTERPOLATED_DERIVATIVES_H
#define __INTERPOLATED_DERIVATIVES_H

#include "LogArgDerivatives.h"
#include "mass_metallicity_interp.h"
#include "../Core/Functions.h"
#include <vector>

namespace StellarEvolution {

    ///\brief Derivative class for stellar quantities which are interpolated
    ///age, mass and metallicity.
    ///
    ///\ingroup StellarSystem_group
    class InterpolatedDerivatives : public LogArgDerivatives {
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
        ///\brief Create an object that interpolates derivatives from evolution
        ///tracks.
        ///
        ///If age is specified the input derivatives are assumed to be with
        ///respect to ln(age), while derivatives always with respect to age are
        ///output.
        InterpolatedDerivatives(
            double mass,
            double metallicity,
            std::vector<const FunctionDerivatives*> *derivatives,
            const alglib::real_1d_array &interp_masses,
            const alglib::real_1d_array &interp_metallicities,
            double age = NaN,
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
