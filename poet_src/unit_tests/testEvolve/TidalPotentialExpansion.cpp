/**\file
 *
 * \brief Define the functions that sums up the expansion series.
 *
 * \ingroup UnitTests_group
 */

#include "TidalPotentialExpansion.h"

namespace Evolve {

    double TidalPotentialExpansion::tidal_term(int mprime,
                                               double radial_distance,
                                               double azimuthal_angle,
                                               double polar_angle,
                                               double orbital_phase) const
    {
        std::complex<double> result = 0.0;
        for(int m=-2; m<=2; ++m) {
            double no_deriv, inclination_deriv, eccentricity_deriv;
            __expansion_coef(__eccentricity,
                             m,
                             mprime,
                             no_deriv,
                             inclination_deriv,
                             eccentricity_deriv);
            assert(std::isfinite(no_deriv));
            result += (
                no_deriv
                *
                std::pow(radial_distance, 2)
                *
                boost::math::spherical_harmonic(
                    2,
                    m,
                    polar_angle,
                    azimuthal_angle
                )
                *
                std::complex<double>(std::cos(mprime * orbital_phase),
                                     -std::sin(mprime * orbital_phase))
            );
            assert(std::isfinite(result.real()));
            assert(std::isfinite(result.imag()));
        }
        return result.real();
    }

    double TidalPotentialExpansion::evaluate_spherical_coords(
        double radial_distance,
        double azimuthal_angle,
        double polar_angle,
        double time
    )
    {
        assert(radial_distance >= 0);
        assert(azimuthal_angle >= 0);
        assert(azimuthal_angle < 2 * M_PI);
        assert(polar_angle >= 0);
        assert(polar_angle <= M_PI);

        __expansion_coef.configure(__inclination);

        double orbital_phase = (
            2.0 * M_PI * time
            /
            EccentricOrbit(__primary_mass,
                           __secondary_mass,
                           __semimajor,
                           __eccentricity).orbital_period()
        ) + __arg_of_periapsis;

        double potential_norm = -(
            Core::AstroConst::G
            *
            __secondary_mass * Core::AstroConst::solar_mass
            /
            std::pow(__semimajor, 3)
        ) /  Core::AstroConst::solar_radius;

        double result = 0.0;
        int mprime_range = (static_cast<int>(__expansion_coef.current_e_order())
                            +
                            2);
        for(
            int mprime = -mprime_range;
            mprime <= mprime_range;
            ++mprime
        ) {
            result += tidal_term(mprime,
                                 radial_distance,
                                 azimuthal_angle,
                                 polar_angle,
                                 orbital_phase);
        }
        return potential_norm * result;
    }

}//End Evolve namespace
