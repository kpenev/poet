/**\file
 *
 * \brief Declare a dissipative zone class where the mass and radius evolution
 * are given by a polynomial.
 *
 * \ingroup UnitTests_group
 */

#ifndef __POLYNOMIAL_EVOLUTION_ZONE_H
#define __POLYNOMIAL_EVOLUTION_ZONE_H

#include "../../Evolve/DissipatingZone.h"

#include <gsl/gsl_poly.h>

namespace Evolve {

    class PolynomialEvolutionZone : public DissipatingZone {
    private:
        std::valarray<double>
            ///The coefficients of the polynomial for the mass evolution.
            __mass_coefficients,

            ///The coefficients of the polynomial for the radius evolution.
            __radius_coefficients,
            
            ///\brief The coefficients of the polynomial for the moment of
            ///inertia evolution.
            __inertia_coefficients;

        std::valarray<double>
            ///\brief The outer mass ofthe zone at the age of the last call to
            ///configure() and its first and second derivatives.
            __current_mass,

            ///\brief The outer radius ofthe zone at the age of the last call to
            ///configure() and its first and second derivatives.
            __current_radius,
            
            ///\brief The moment of inertia ofthe zone at the age of the last
            ///call to configure() and its first and second derivatives.
            __current_inertia;

        ///\brief Evaluate either the given order derivative of either the mass
        ///or the radius polynomial.
        double evaluate_polynomial(
            ///The coefficienst of the polynomial to evaluate.
            const std::valarray<double> &coefficients,

            ///The age at which to evaluate the polynomial.
            double age,

            ///The derivative order desired.
            int deriv_order
        ) const;
    public:
        ///Construct a polynomial evolution zone with the given coefficients.
        PolynomialEvolutionZone(
            ///The coefficients of the polynomial for the mass evolution.
            const std::valarray<double> &mass_coefficients,

            ///The coefficients of the polynomial for the radius evolution.
            const std::valarray<double> &radius_coefficients,

            ///The coefficients of the polynomial for the moment of inertia
            ///evolution.
            const std::valarray<double> &inertia_coefficients
        ):
            __mass_coefficients(mass_coefficients),
            __radius_coefficients(radius_coefficients),
            __inertia_coefficients(inertia_coefficients),
            __current_mass(Core::NaN, 3),
            __current_radius(Core::NaN, 3),
            __current_inertia(Core::NaN, 3)
        {}

        ///See DissipatingZone::configure().
        virtual void configure(bool initialize,
                               double age,
                               double orbital_frequency,
                               double eccentricity,
                               double orbital_angmom,
                               double spin,
                               double inclination,
                               double periapsis,
                               bool spin_is_frequency);

        ///See DissipatingZone::outer_mass(int).
        double outer_mass(int deriv_order=0) const
        {return __current_mass[deriv_order];}

        ///See DissipatingZone::outer_mass(double, int).
        double outer_mass(double age, int deriv_order=0) const
        {return evaluate_polynomial(__mass_coefficients, age, deriv_order);}

        ///See DissipatingZone::outer_radius(int).
        double outer_radius(int deriv_order=0) const
        {return __current_radius[deriv_order];}

        ///See DissipatingZone::outer_radius(double, int).
        double outer_radius(double age, int deriv_order=0) const
        {return evaluate_polynomial(__radius_coefficients, age, deriv_order);}

        ///See DissipatingZone::moment_of_inertia(int).
        virtual double moment_of_inertia(int deriv_order=0) const
        {return __current_inertia[deriv_order];}

        ///See DissipatingZone::moment_of_inertia(double, int).
        virtual double moment_of_inertia(double age, int deriv_order=0) const
        {return evaluate_polynomial(__inertia_coefficients, age, deriv_order);}

        ///See DissipatingZone::love_coefficient()
        double love_coefficient(int, int, Dissipation::QuantityEntry) const
        {return 0.0;}
    };//End PolynomialEvolutionZone class.

} //End Evolve namespace.

#endif
