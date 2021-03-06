/**\file
 *
 * \brief Declares a OneArgumentDiffFunction sub-class giving the evolution of
 * the stellar convective zone angular momentum when some function of it evolves
 * linearly with time.
 *
 * \ingroup UnitTests_group
 */

#ifndef __INVERSE_LINEAR_LCONV_EVOLUTION_H
#define __INVERSE_LINEAR_LCONV_EVOLUTION_H

#include "InverseFunction.h"

template<class LINEAR_QUANTITY_TYPE>
class InverseLinearLconvEvolution : public Core::OneArgumentDiffFunction {
private:
    double 
        ///The lifetime of the protoplanetary disk.
        __disk_lifetime,

        ///\brief The rate at which the linear quantity (see
        ///for example Oblique10LinearQuantity) evolves.
        __evolution_rate;

    ///\brief A function of the convective zone angular momentum which evolves
    ///linearly with time.
    LINEAR_QUANTITY_TYPE __linear_quantity;

    ///Evaluates to Lconv given \f$ \frac{3\pi}{5} T_0 \Delta_{10}(t - tdisk)\f$
    InverseFunction __find_lconv;

public:
    InverseLinearLconvEvolution(
        ///See __disk_lifetime member.
        double disk_lifetime,

        ///The total angular momentum (does not evolve).
        double total_angmom,

        ///The orbital angular momentum (does not evolve).
        double orbital_angmom,

        ///The initial convective angular momentum.
        double initial_conv_angmom,

        ///See __evolution_rate member.
        double evolution_rate
    ) :
        __disk_lifetime(disk_lifetime),
        __evolution_rate(evolution_rate),
        __linear_quantity(total_angmom, orbital_angmom, initial_conv_angmom),
        __find_lconv(__linear_quantity,
                     (
                         (total_angmom - orbital_angmom)
                         *
                         (1.0 + 100.0 * std::numeric_limits<double>::epsilon())
                     ),
                     initial_conv_angmom)
    {}

    ///\brief Return the expected value for the stellar convective angular
    ///momentum at the given age.
    double operator()(double age) const
    {
        return __find_lconv(
            std::min(0.0, -(age - __disk_lifetime) * __evolution_rate)
        );
    }

    ///The disk lifetime specified at construction.
    double disk_lifetime() const {return __disk_lifetime;}

    double range_high() const
    {return __find_lconv.range_high() / __evolution_rate;}

    double range_low() const
    {return -Core::Inf;}

    ///\brief An iterator over the abscissas where the function takes
    ///the given y value.
    Core::InterpSolutionIterator crossings(double = 0) const
    {
        throw Core::Error::Runtime(
            "Finding all solutinos of ConservedLELconvEvolution not supported!"
        );
    };

    ///\brief Returns a pointer to the derivative of the function.
    ///
    ///Result must be deleted when no longer needed.
    ///
    ///The use of a pointer allows avoiding potentially expensive copy
    ///opertaions.
    const Core::FunctionDerivatives *deriv(double age) const
    {
        if(age <= __disk_lifetime)
            return new Core::CubicSplineDerivatives(__find_lconv(0.0),
                                                    0.0,
                                                    0.0);
        const Core::FunctionDerivatives *scaled_deriv = __find_lconv.deriv(
            -(age - __disk_lifetime) * __evolution_rate
        );
        double value = scaled_deriv->order(0),
               order1 = -scaled_deriv->order(1) * __evolution_rate,
               order2 = scaled_deriv->order(2) * std::pow(__evolution_rate, 2);
        delete scaled_deriv;
        return new Core::CubicSplineDerivatives(value, order1, order2);
    }
};


#endif
