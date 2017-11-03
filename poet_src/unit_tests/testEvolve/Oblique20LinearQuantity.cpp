/**\file
 *
 * \brief Implements some methods of the Oblique20LinearQuantity class.
 *
 * \ingroup UnitTests_group
 */

#include "Oblique20LinearQuantity.h"

double Oblique20LinearQuantity::indefinite_integral(double star_angmom) const
{
    double s2 = std::pow(star_angmom, 2),
           t2 = std::pow(__total_angmom, 2),
           t3 = t2 * __total_angmom,
           atanh_good_arg = star_angmom / (1.0 + __total_angmom),
           atanh_bad_arg = star_angmom / (1.0 - __total_angmom);
    return (
        (
            2.0 * star_angmom * __total_angmom
            *
            (1.0 - s2 * t2 + t2 * t2 - (s2 + 2.0 * t2))
            /
            (1.0 + std::pow(s2 - t2, 2) - 2.0 * (s2 + t2))
        ) 
        +
        (t3 - 1.0) * (std::log(std::abs(atanh_bad_arg + 1.0))
                      -
                      std::log(std::abs(atanh_bad_arg - 1.0))) / 2.0
        +
        (1.0 + t3) * std::atanh(atanh_good_arg)
    ) / (2.0 * t3);
}

Oblique20LinearQuantity::Oblique20LinearQuantity(double total_angmom,
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
    assert(!std::isnan(__initial_star_angmom));
    assert(!std::isnan(__initial_indefinite_integral));
}

double Oblique20LinearQuantity::operator()(double star_angmom) const
{
    return (indefinite_integral(star_angmom / __angmom_scale)
            -
            __initial_indefinite_integral) * __angmom_scale;
}

const Core::FunctionDerivatives *Oblique20LinearQuantity::deriv(
    double
) const
{
    throw Core::Error::Runtime("Differentiating m=2, m' = 0 linearly evolving "
                               "function of the stellar angular momentum is not"
                               " supported.");
}
