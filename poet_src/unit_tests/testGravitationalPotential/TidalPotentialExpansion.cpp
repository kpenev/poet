/**\file
 *
 * \brief Define the functions that sums up the expansion series.
 *
 * \ingroup UnitTests_group
 */

#include "TidalPotentialExpansion.h"

namespace testGravitationalPotential {

    double TidalPotentialExpansion::tidal_term(int mprime,
                                               double radiual_distance,
                                               double azimuthal_angle,
                                               double polar_angle,
                                               double time)
    {
        result = 0.0;
        azimuthal_correction = 2.0 * M_PI * mprime * time / orbit.orbital_period();
        for(int m=-2; m<=2; ++m) {
            double no_deriv, inclination_deriv, eccentricity_deriv;
            __expansion_coef(__eccentricity,
                             m,
                             mprime,
                             no_deriv,
                             inclination_deriv,
                             eccentricity_deriv);
            result += (
                no_deriv
                *
                std::pow(radial_distance, 2)
                *
                spherical_harmonic_r(2,
                                     m,
                                     polar_angle,
                                     azimuthal_angle - azimuthal_correction / m)
            );
        }
        return result;
    }

    double TidalPotentialExpansion::evaluate_spherical_coords(
        double radiual_distance,
        double azimuthal_angle,
        double polar_angle,
        double time
    )
    {
        assert(radial_distance >= 0);
        assert(azimumthal_angel >= 0);
        assert(azimuthal_angle < 2 * M_PI);
        assert(polar_angle >= 0);
        assert(polar_angle <= M_PI);

        __expansion_coef.configure(__inclination);

        double potential_norm = -(
            Core::AstroConst::G
            *
            __secondary_mass * Core::AstroConst::solar_mass
            /
            std::pow(__semimajor * Core::AstroConst::solar_radius, 3)
        );
    }

}//End testGravitationalPotential namespace
