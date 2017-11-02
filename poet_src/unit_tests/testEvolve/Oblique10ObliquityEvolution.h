/**\file
 *
 * \brief Declares a OneArgumentDiffFunction sub-class giving the evolution of
 * the obliquity under the m = 1, m' = 0 term.
 *
 * \ingroup UnitTests_group
 */

#ifndef __OBLIQUE_1_0_OBLIQUITY_EVOLUTION_H
#define __OBLIQUE_1_0_OBLIQUITY_EVOLUTION_H

#include "Oblique10LconvEvolution.h"

class Oblique10ConvObliquityEvolution : public Core::OneArgumentDiffFunction {
private:
    ///The evolution of the convective zone angular momentum.
    const Oblique10LconvEvolution &__lconv_evol;

    ///The orbital angular momentum (does not evolve).
    double __orbital_angmom;

    ///The difference between the squares of the total and orbital angular
    ///momenta.
    double __total2_minus_orbital2;

    ///The value of the obliquity for the given convective angular momentum.
    double value(double lconv) const
    {
        return std::acos((__total2_minus_orbital2 - lconv * lconv)
                         /
                         (2.0 * lconv * __orbital_angmom));
    }

public:
    Oblique10ConvObliquityEvolution(
        ///See __lconv_evol member.
        const Oblique10LconvEvolution &lconv_evol,

        ///See __orbital_angmom member.
        double orbital_angmom,

        ///The total angular momentum in the system (does not evolve).
        double total_angmom
    ) :
        __lconv_evol(lconv_evol),
        __orbital_angmom(orbital_angmom),
        __total2_minus_orbital2(std::pow(total_angmom, 2)
                                -
                                std::pow(orbital_angmom, 2))
    {}

    double operator()(double age) const
    {
        return (age  < __lconv_evol.disk_lifetime()
                ? 0.0
                : value(__lconv_evol(age)));
    }

    double range_high() const {return __lconv_evol.range_high();}
    double range_low() const {return __lconv_evol.range_low();}

    Core::InterpSolutionIterator crossings(double = 0) const
    {
        throw Core::Error::Runtime(
            "Finding all solutinos of Oblique10ObliquityEvolution not supported!"
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
        throw Core::Error::Runtime(
            "Derivatives not implemented for 1-0 obliquity evolution."
        );
        if(age < __lconv_evol.disk_lifetime())
            return new Core::CubicSplineDerivatives(0.0, 0.0, 0.0);
        const Core::FunctionDerivatives *lconv_deriv = __lconv_evol.deriv(age);
        double obliquity = value(lconv_deriv->order(0)),
               order1 =  - 1.0 / std::sin(obliquity) * (
                   __total2_minus_orbital2
                   /
                   (2.0 * __orbital_angmom * std::pow(lconv_deriv->order(0), 2))
                   -
                   1.0 / (2.0 * __orbital_angmom)
               ),
               order2;
        delete lconv_deriv;
        return new Core::CubicSplineDerivatives(obliquity, order1, order2);
    }

};

class Oblique10RadObliquityEvolution : public Core::OneArgumentDiffFunction {
private:
    ///The evolution of the convective zone angular momentum.
    const Oblique10LconvEvolution &__lconv_evol;

    ///The evolution of the convective zone obliquity.
    const Oblique10ConvObliquityEvolution &__conv_obliq_evol;

    ///The orbital angular momentum magnitude (does not evolve).
    double __orbital_angmom;

    ///The obliquity at the time the disk dissipates.
    double __initial_obliquity;

    ///\brief The initial component of the convective zone angular momentum
    ///perpendicular to the orbital angular momentum.
    double __initial_lconv_perp;

    ///The radiative zone obliquity for the given parameters.
    double value(
        ///The convective zone angular momentum magnitude.
        double lconv,

        ///The convective zone obliqity
        double conv_obliq
    ) const
    {
        double lconv_perp = lconv * std::sin(conv_obliq),
               lconv_par = lconv * std::cos(conv_obliq),
               repeated = (std::pow(__orbital_angmom, 2)
                           +
                           std::pow(lconv, 2)
                           +
                           2.0 * __orbital_angmom * lconv_par);
                
        return __initial_obliquity - std::acos(
            (
                (__orbital_angmom + lconv_par)
                *
                std::sqrt(
                    repeated
                    -
                    std::pow(__initial_lconv_perp, 2)
                )
                +
                __initial_lconv_perp * lconv_perp
            )
            /
            repeated
        );
    }

public:
    Oblique10RadObliquityEvolution(
        ///See __lconv_evol member.
        const Oblique10LconvEvolution &lconv_evol,

        ///See __conv_obliq_evol member.
        const Oblique10ConvObliquityEvolution &conv_obliq_evol,

        ///See __orbital_angmom member.
        double orbital_angmom
    ) :
        __lconv_evol(lconv_evol),
        __conv_obliq_evol(conv_obliq_evol),
        __orbital_angmom(orbital_angmom),
        __initial_obliquity(conv_obliq_evol(__lconv_evol.disk_lifetime())),
        __initial_lconv_perp(lconv_evol(__lconv_evol.disk_lifetime())
                             *
                             std::sin(__initial_obliquity))
    {
        std::cerr << "Initial obliquity = " << __initial_obliquity << std::endl;
    }

    double operator()(double age) const
    {
        return (age < __lconv_evol.disk_lifetime()
                ? 0.0
                : value(__lconv_evol(age), __conv_obliq_evol(age)));
    }

    double range_high() const {return std::min(__lconv_evol.range_high(),
                                               __conv_obliq_evol.range_high());}
    double range_low() const {return std::max(__lconv_evol.range_low(),
                                              __conv_obliq_evol.range_low());}

    ///\brief An iterator over the abscissas where the function takes
    ///the given y value.
    Core::InterpSolutionIterator crossings(double = 0) const
    {
        throw Core::Error::Runtime(
            "Solving Oblique10RadObliquityEvolution not implemented!"
        );
    };

    ///\brief Returns a pointer to the derivative of the function.
    ///
    ///Result must be deleted when no longer needed.
    ///
    ///The use of a pointer allows avoiding potentially expensive copy
    ///opertaions.
    const Core::FunctionDerivatives *deriv(double) const
    {
        throw Core::Error::Runtime(
            "Differentiating Oblique10RadObliquityEvolution not implemented!"
        );
    }
};

#endif
