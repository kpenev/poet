#include "InterpolatedDerivatives.h"

namespace StellarEvolution {

    double InterpolatedDerivatives::calc_deriv(unsigned deriv_order) const
    {
        if(deriv_order > 2) return 0.0;

        alglib::real_1d_array interp_values;
        interp_values.setlength(__interp_deriv->size());
        for(unsigned i = 0; i < __interp_deriv->size(); ++i)
            interp_values[i] = (*__interp_deriv)[i]->order(deriv_order);

        return mass_feh_interp(__interp_masses,
                               __interp_feh,
                               interp_values,
                               __stellar_mass,
                               __stellar_feh);
    }

    InterpolatedDerivatives::InterpolatedDerivatives(
        double mass,
        double feh,
        std::vector<const FunctionDerivatives*> *derivatives,
        const alglib::real_1d_array &interp_masses,
        const alglib::real_1d_array &interp_feh,
        double age,
        bool log_quantity,
        bool delete_derivatives
    ) :
        LogDerivatives(age, log_quantity),
        __stellar_mass(mass),
        __stellar_feh(feh),
        __interp_deriv(derivatives),
        __interp_masses(interp_masses),
        __interp_feh(interp_feh),
        __delete_derivatives(delete_derivatives)
    {
        assert(static_cast<long>(derivatives->size())
               == 
               interp_masses.length() * interp_feh.length());
    }

} //End of StellarEvolution namespace.
