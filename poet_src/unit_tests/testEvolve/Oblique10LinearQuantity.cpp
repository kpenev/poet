/**\file
 *
 * \brief Implements some methods of the Oblique10LinearQuantity class.
 *
 * \ingroup UnitTests_group
 */

#include "Oblique10LinearQuantity.h"

double Oblique10LinearQuantity::indefinite_integral(double star_angmom) const
{
    double t2minus1 = std::pow(__total_angmom, 2) - 1.0;
    return (
        (2.0 * star_angmom) / (t2minus1 - std::pow(star_angmom, 2))
        +
        (1.0 / __total_angmom - 1.0) * std::atanh((__total_angmom - 1.0)
                                                  /
                                                  star_angmom)
        -
        2.0 * atanh(star_angmom / std::sqrt(t2minus1)) / std::sqrt(t2minus1)
        +
        (0.5 / __total_angmom + 0.5) * (
            std::log((1.0 + __total_angmom) / star_angmom + 1.0)
            -
            std::log((1.0 + __total_angmom) / star_angmom - 1.0)
        )
    );
}

Oblique10LinearQuantity::Oblique10LinearQuantity(double total_angmom,
                                                 double orbital_angmom,
                                                 double initial_star_angmom) :
    __total_angmom(total_angmom / orbital_angmom),
    __initial_star_angmom(initial_star_angmom / orbital_angmom),
    __angmom_scale(orbital_angmom),
    __initial_indefinite_integral(indefinite_integral(__initial_star_angmom))
{
    assert(total_angmom > 0);
    assert(orbital_angmom > 0);
    assert(initial_star_angmom > 0);
    assert(total_angmom <= orbital_angmom + initial_star_angmom);
    assert(orbital_angmom <= total_angmom + initial_star_angmom);
    assert(initial_star_angmom <= orbital_angmom + total_angmom);
    assert(std::pow(total_angmom, 2) - std::pow(orbital_angmom, 2)
           >
           std::pow(initial_star_angmom, 2));
    assert(!std::isnan(__initial_indefinite_integral));
}

double Oblique10LinearQuantity::operator()(double star_angmom) const
{
    return (indefinite_integral(star_angmom / __angmom_scale)
            -
            __initial_indefinite_integral) * __angmom_scale;
}

const Core::FunctionDerivatives *Oblique10LinearQuantity::deriv(
    double star_angmom
) const
{
    star_angmom /= __angmom_scale;
    double s2 = std::pow(star_angmom, 2),
           t2 = std::pow(__total_angmom, 2);
    return new Core::CubicSplineDerivatives(
        operator()(star_angmom),
        (
            (__total_angmom + star_angmom - 1.0)
            *
            (1.0 + star_angmom - __total_angmom)
            *
            (1.0 + __total_angmom - star_angmom)
            *
            (1.0 + __total_angmom + star_angmom)
            *
            std::pow(1.0 + s2 - t2, 2)
            /
            std::pow(4.0 * s2, 2)
        ) * __angmom_scale,
        (
            (1.0 - s2 - t2)
            *
            (1.0 + s2 - t2)
            *
            (1.0 + s2 * s2 - 2.0 * (1.0 + s2) * t2 + t2 * t2)
            /
            (4.0 * s2 * s2 * star_angmom)
        ) * __angmom_scale
    );
}
