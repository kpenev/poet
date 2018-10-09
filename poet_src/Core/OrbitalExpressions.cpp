/**\file
 *
 * \brief The definitions of the orbital expression functions.
 */

#include "OrbitalExpressions.h"

namespace Core {

    double factorial(unsigned n)
    {
        double result=1;
        for(unsigned i=2; i<=n; ++i) result*=i;
        return result;
    }

    double orbital_angular_velocity(double m1,
                                    double m2,
                                    double semimajor,
                                    bool deriv)
    {
        return (
            (deriv ? -1.5 : 1.0)
            *
            std::sqrt(
                Core::AstroConst::G * (m1 + m2) * Core::AstroConst::solar_mass
                /
                std::pow(
                    semimajor * Core::AstroConst::solar_radius,
                    (deriv ? 5 : 3)
                )
            )
            * 
            Core::AstroConst::day
        );
    }

    double orbital_energy(double m1,
                          double m2,
                          double semimajor,
                          unsigned deriv_order)
    {
        return (deriv_order%2 ? 1 : -1) * m1 * m2 * factorial(deriv_order)
               /
               (2.0 * std::pow(semimajor, static_cast<int>(deriv_order + 1)))
               *
               Core::AstroConst::G * Core::AstroConst::solar_mass
               *
               std::pow(Core::AstroConst::day, 2)
               /
               std::pow(Core::AstroConst::solar_radius, 3);
    }

    double orbital_angular_momentum(double m1,
                                    double m2,
                                    double semimajor,
                                    double eccentricity)
    {
        return m1 * m2
               *
               std::sqrt(
                   semimajor*(1.0 - std::pow(eccentricity, 2)) / (m1 + m2)
                   *
                   Core::AstroConst::G * Core::AstroConst::solar_mass
                   /
                   std::pow(Core::AstroConst::solar_radius, 3)
               )
               *
               Core::AstroConst::day;
    }

}//End Core namespace.
