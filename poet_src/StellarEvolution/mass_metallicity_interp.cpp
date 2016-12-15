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

        if(interp_values.length() == 1) {
            assert(stellar_mass == interp_masses[0]);
            assert(stellar_metallicity == interp_metallicities[0]);
            return interp_values[0];
        } else if(interp_masses.length() == 1) {
            alglib::spline1dinterpolant spline;
            alglib::spline1dbuildcubic(interp_metallicities,
                                       interp_values,
                                       spline);
            return alglib::spline1dcalc(spline, stellar_metallicity);
        } else if(interp_metallicities.length() == 1) {
            alglib::spline1dinterpolant spline;
            alglib::spline1dbuildcubic(interp_masses,
                                       interp_values,
                                       spline);
            return alglib::spline1dcalc(spline, stellar_mass);
        } else {
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

}
