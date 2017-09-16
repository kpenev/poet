/**\file
 * 
 * \brief Define a single function performing mass-[Fe/H] interpolations.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __MASS_METALLICITY_INTERP
#define __MASS_METALLICITY_INTERP

#include "alglib/src/interpolation.h"
#include <cassert>

namespace StellarEvolution {

    ///Perform a bi-cubic spline interpolation of a single quantity.
    double mass_feh_interp(
        ///The masses of the stelar models on which to base the 
        ///interpolation in \f$M_\odot\f$
        const alglib::real_1d_array &interp_masses,

        ///The [Fe/H] of the stellar models on which to base the
        ///interpolation.
        const alglib::real_1d_array &interp_feh,

        ///The values of the quantity being interpolated on the grid defined 
        ///by \p interp_masses and \p interp_feh.
        const alglib::real_1d_array &interp_values,

        ///The stellar mass to which to interpolate in \f$M_\odot\f$.
        double stellar_mass,

        ///The stellar \f$[Fe/H]\f$ to which to interpolate.
        double stellar_feh
    );

}

#endif
