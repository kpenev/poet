/**\file
 *
 * \brief Declares a function of the stellar angular momentum that evolves
 * linearly with time when only the m = 2, m' = 0 term is active for orbits with
 * arbitrary initial obliquity.
 *
 * \ingroup UnitTests_group
 */

#ifndef __OBLIQUE_2_0_LINEAR_QUANTITY_H
#define __OBLIQUE_2_0_LINEAR_QUANTITY_H

#include "../../Core/Functions.h"
#include "../../Core/InterpSolutionIterator.h"

///\brief A function of the stellar angular momentum expected to evolve linearly
///with time under the m = 1, m' = 0 term.
class Oblique20LinearQuantity : public Core::OneArgumentDiffFunction {
private:
    double
        ///\brief The Magnitude of the total angular momentum in the system in
        ///units of the orbital angular momentum (conserved).
        __total_angmom,

        ///\brief The Magnitude of the initial stellar spin angular momentum in
        //units of the orbital angular momentum.
        __initial_star_angmom,
 
        ///\brief The initial orbital angular momentum (everything is scaled by
        ///this quantity).
        __angmom_scale,
 
        ///The value of indefinite integral(__initial_star_angmom)
        __initial_indefinite_integral;

    ///\brief Return the real part of the indefinite integral of the inverse of
    ///the rate of change of the stellar angular momentum.
    double indefinite_integral(double star_angmom) const;

public:
    Oblique20LinearQuantity(double total_angmom,
                            double orbital_angmom,
                            double initial_star_angmom);

    double operator()(double star_angmom) const;

    double range_high() const {return __initial_star_angmom * __angmom_scale;}
    double range_low() const {return (__total_angmom - 1.0) * __angmom_scale;}

    ///\brief An iterator over the abscissas where the function takes
    ///the given y value.
    Core::InterpSolutionIterator crossings(double = 0) const
    {
        throw Core::Error::Runtime(
            "Finding all solutinos of Oblique20LinearQuantity not supported!"
        );
    };

    ///\brief Returns a pointer to the derivative of the function.
    ///
    ///Result must be deleted when no longer needed.
    ///
    ///The use of a pointer allows avoiding potentially expensive copy
    ///opertaions.
    const Core::FunctionDerivatives *deriv(double star_angmom) const;

}; //End Oblique20LinearQuantity class.

#endif
