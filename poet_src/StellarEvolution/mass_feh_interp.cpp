#define BUILDING_LIBRARY
#include "mass_feh_interp.h"

namespace StellarEvolution {

    double mass_feh_interp(
        const alglib::real_1d_array &interp_masses,
        const alglib::real_1d_array &interp_feh,
        const alglib::real_1d_array &interp_values,
        double stellar_mass,
        double stellar_feh
    )
    {
        assert(interp_masses.length() * interp_feh.length()
               ==
               interp_values.length());

        if(interp_values.length() == 1) {
            assert(stellar_mass == interp_masses[0]);
            assert(stellar_feh == interp_feh[0]);
            return interp_values[0];
        } else if(interp_masses.length() == 1) {
            alglib::spline1dinterpolant spline;
            alglib::spline1dbuildcubic(interp_feh,
                                       interp_values,
                                       spline);
            return alglib::spline1dcalc(spline, stellar_feh);
        } else if(interp_feh.length() == 1) {
            alglib::spline1dinterpolant spline;
            alglib::spline1dbuildcubic(interp_masses,
                                       interp_values,
                                       spline);
            return alglib::spline1dcalc(spline, stellar_mass);
        } else {
            alglib::spline2dinterpolant spline;
            alglib::spline2dbuildbicubicv(interp_masses, 
                                          interp_masses.length(), 
                                          interp_feh,
                                          interp_feh.length(),
                                          interp_values,
                                          1,
                                          spline);
            return alglib::spline2dcalc(spline,
                                        stellar_mass,
                                        stellar_feh);
        }
    }

}
