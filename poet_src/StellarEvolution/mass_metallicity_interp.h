/**\file
 * 
 * \brief Define a single function performing mass-metallicity 
 * interpolations.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __MASS_METALLICITY_INTERP
#define __MASS_METALLICITY_INTERP

#include "alglib/src/interpolation.h"
#include <cassert>

namespace StellarEvolution {

    ///Perform a bi-cubic spline interpolation of a single quantity.
    double mass_metallicity_interp(
        ///The masses of the stelar models on which to base the 
        ///interpolation in \f$M_\odot\f$
        const alglib::real_1d_array &interp_masses,

        ///The metallicities of the stellar models on which to base the
        ///interpolation in \f$[Fe/H]\f$.
        const alglib::real_1d_array &interp_metallicities,

        ///The values of the quantity being interpolated on the grid defined 
        ///by \p interp_masses and \p interp_metallicities.
        const alglib::real_1d_array &interp_values,

        ///The stellar mass to which to interpolate in \f$M_\odot\f$.
        double stellar_mass,

        ///The stellar metallicity to which to interpolate in \f$[Fe/H]\f$.
        double stellar_metallicity
    );

}

#endif
