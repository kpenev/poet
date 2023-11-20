/**\file
 *
 * \brief Definitions of some of the methods of EvolvingStellarZone.
 *
 * \ingroup Star_group
 */

#define BUILDING_LIBRARY
#include "EvolvingStellarZone.h"

namespace Star {

    void EvolvingStellarZone::reset_current_quantities()
    {
        for(size_t i=0; i<__current_age_quantities.size(); ++i)
            if(__current_age_quantities[i]) {
                delete __current_age_quantities[i];
                __current_age_quantities[i]=NULL;
            }
    }

    double EvolvingStellarZone::current_age_quantity(size_t quantity,
                                                     unsigned deriv_order) const
    {
        assert(quantity < __evolving_quantities.size());

        if(__current_age_quantities[quantity]==NULL) {
            __current_age_quantities[quantity] =
                __evolving_quantities[quantity]->deriv(__current_age);
        }
        return __current_age_quantities[quantity]->order(deriv_order);
    }

    double EvolvingStellarZone::any_age_quantity(size_t quantity,
                                                 double age,
                                                 unsigned deriv_order) const
    {
        if(deriv_order == 0)
            return (*(__evolving_quantities[quantity]))(age);
        else {
            const Core::FunctionDerivatives
                *deriv = __evolving_quantities[quantity]->deriv(age);
            double result = deriv->order(deriv_order);
            delete deriv;
            return result;
        }
    }

    void EvolvingStellarZone::configure(bool initialize,
                                        double age,
                                        double orbital_frequency,
                                        double eccentricity,
                                        double orbital_angmom,
                                        double spin,
                                        double inclination,
                                        double periapsis,
                                        bool spin_is_frequency)
    {
        if(__current_age != age) {
            __current_age = age;
            reset_current_quantities();
        }
        BrokenPowerlawPhaseLagZone::configure(initialize,
                                              age,
                                              orbital_frequency,
                                              eccentricity,
                                              orbital_angmom,
                                              spin,
                                              inclination,
                                              periapsis,
                                              spin_is_frequency);
    }

    EvolvingStellarZone::~EvolvingStellarZone()
    {
        reset_current_quantities();

        for(size_t i = 0; i < __evolving_quantities.size(); ++i)
            if(__evolving_quantities[i]) delete __evolving_quantities[i];
    }

    void EvolvingStellarZone::reached_critical_age(double age)
    {
#ifndef NDEBUG
        std::cerr << "Adjusting zone quantity interpolation per age = "
                  << age
                  << std::endl;
#endif
        for(size_t i = 0; i < __evolving_quantities.size(); ++i) {
#ifndef NDEBUG
            std::cerr << "Quantity " << i << " next discontinuity at: "
                      << __evolving_quantities[i]->next_discontinuity();
#endif

            while(__evolving_quantities[i]->next_discontinuity() <= age) {
                __evolving_quantities[i]->enable_next_interpolation_region();
#ifndef NDEBUG
                std::cerr << ", switching to next region up to t = "
                          << __evolving_quantities[i]->next_discontinuity();
#endif
            }
#ifndef NDEBUG
            std::cerr << std::endl;
#endif
        }
        BrokenPowerlawPhaseLagZone::reached_critical_age(age);
    }

    double EvolvingStellarZone::next_stop_age() const
    {
        double result = Evolve::BrokenPowerlawPhaseLagZone::next_stop_age();
        for(size_t i = 0; i < __evolving_quantities.size(); ++i)
            result = std::min(
                result,
                __evolving_quantities[i]->next_discontinuity()
            );
        return result;
    }

    double EvolvingStellarZone::min_interp_age() const
    {
        double result = -Core::Inf;
        for(
            std::vector<
                const StellarEvolution::EvolvingStellarQuantity*
            >::const_iterator quantity_iter = __evolving_quantities.begin();
            quantity_iter != __evolving_quantities.end();
            ++quantity_iter
        )
            result = std::max(result, (*quantity_iter)->range_low());
        return result;
    }

    void EvolvingStellarZone::select_interpolation_region(double age) const
    {
        for(
            std::vector<
                const StellarEvolution::EvolvingStellarQuantity*
            >::const_iterator quantity_iter = __evolving_quantities.begin();
            quantity_iter != __evolving_quantities.end();
            ++quantity_iter
        )
            (*quantity_iter)->select_interpolation_region(age);
    }

}//End Star namespace.
