#include "mass_metallicity_interp.h"

namespace StellarEvolution {

    double mass_metallicity_interp(
        const alglib::real_1d_array &interp_masses,
        const alglib::real_1d_array &interp_metallicities,
        const alglib::real_1d_array &interp_values,
        double stellar_mass,
        double stellar_metallicity
    )
    {
        assert(interp_masses.length() * interp_metallicities.length()
               ==
               interp_values.length());
        alglib::spline2dinterpolant spline;
        alglib::spline2dbuildbicubicv(interp_masses, 
                                      interp_masses.length(), 
                                      interp_metallicities,
                                      interp_metallicities.length(),
                                      interp_values,
                                      1,
                                      spline);
        return alglib::spline2dcalc(spline,
                                    stellar_mass,
                                    stellar_metallicity);
    }

}
