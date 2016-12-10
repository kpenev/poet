#include "InterpolatedDerivatives.h"

namespace StellarEvolution {

    double InterpolatedDerivatives::calc_deriv(unsigned deriv_order) const
    {
        if(deriv_order > 2) return 0.0;

        alglib::real_1d_array interp_values;
        interp_values.setlength(__interp_deriv->size());
        for(unsigned i = 0; i < __interp_deriv->size(); ++i)
            interp_values[i] = (*__interp_deriv)[i]->order(deriv_order);

        return mass_metallicity_interp(__interp_masses,
                                       __interp_metallicities,
                                       interp_values,
                                       __stellar_mass,
                                       __stellar_metallicity);
    }

    InterpolatedDerivatives::InterpolatedDerivatives(
        double mass,
        double metallicity,
        std::vector<const FunctionDerivatives*> *derivatives,
        const alglib::real_1d_array &interp_masses,
        const alglib::real_1d_array &interp_metallicities,
        double age,
        bool delete_derivatives
    ) :
        LogArgDerivatives(age),
        __stellar_mass(mass),
        __stellar_metallicity(metallicity),
        __interp_deriv(derivatives),
        __interp_masses(interp_masses),
        __interp_metallicities(interp_metallicities),
        __delete_derivatives(delete_derivatives)
    {
        assert(derivatives->size()
               == 
               interp_masses.length() * interp_metallicities.length());
    }

} //End of StellarEvolution namespace.
