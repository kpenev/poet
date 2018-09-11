/**\file
 *
 * \brief Define the non-inlnine methods of PolynomialEvolutionZone.
 *
 * \ingroup UnitTests_group
 */

#include "PolynomialEvolutionZone.h"

namespace Evolve {
    void PolynomialEvolutionZone::configure(bool initialize,
                                            double age,
                                            double orbital_frequency,
                                            double eccentricity,
                                            double orbital_angmom,
                                            double spin,
                                            double inclination,
                                            double periapsis,
                                            bool spin_is_frequency)
    {
        for(int deriv_order = 0; deriv_order <= 2; ++deriv_order) {
            __current_radius[deriv_order] = outer_radius(age,
                                                         deriv_order);
            __current_mass[deriv_order] = outer_mass(age,
                                                     deriv_order);
            __current_inertia[deriv_order] = moment_of_inertia(age,
                                                               deriv_order);
        }

        DissipatingZone::configure(initialize,
                                   age,
                                   orbital_frequency,
                                   eccentricity,
                                   orbital_angmom,
                                   spin,
                                   inclination,
                                   periapsis,
                                   spin_is_frequency);
    }

    double PolynomialEvolutionZone::evaluate_polynomial(
        const std::valarray<double> &coefficients,
        double age,
        int deriv_order
    ) const
    {
        std::valarray<double> gsl_result(deriv_order + 1);
        gsl_poly_eval_derivs(&(coefficients[0]),
                             coefficients.size(),
                             age,
                             &(gsl_result[0]),
                             deriv_order + 1);
        return gsl_result[deriv_order];
    }

}//End Evolve namespace.
