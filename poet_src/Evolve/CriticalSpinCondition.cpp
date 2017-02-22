#include "CriticalSpinCondition.h"

namespace Evolve {

    void CriticalSpinCondition::set_num_subconditions()
    {
        if(
            __critical_below_iter == __critical_spins.end()
            ||
            __critical_above_iter == __critical_spins.end()
        )
            __num_subconditions = 1;
        else
            __num_subconditions = 2;
    }

    double CriticalSpinCondition::derivative(double surf_angmom_deriv,
                                             double wcritical) const
    {
        return (
            (
                (__zone.spin_frequency() < 0 ? -1.0 : 1.0)
                *
                surf_angmom_deriv
                - 
                __zone.moment_of_inertia(1)
                *
                std::abs(__zone.spin_frequency())
            )
            /
            (__zone.moment_of_inertia() * wcritical)
        );
    }

    void CriticalSpinCondition::fill_locked_derivs(
        Core::EvolModeType 
#ifndef NDEBUG
        evol_mode
#endif
        ,
        const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
        assert(evol_mode == Core::BINARY);
        double deriv_factor = (-1.5 * __zone.spin_frequency()
                               *
                               derivatives[0] / orbit[0]);
        unsigned subcond = 0;
        if(__critical_below_iter != __critical_spins.end())
            stop_deriv[subcond++] = deriv_factor * (*__critical_below_iter);
        if(__critical_above_iter != __critical_spins.end())
            stop_deriv[subcond] = deriv_factor * (*__critical_above_iter);
    }

    void CriticalSpinCondition::fill_unlocked_derivs(
        Core::EvolModeType evol_mode,
        const std::valarray<double> &,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
        unsigned angmom_index = 1 + __zone_index + 2 * __body.number_zones();
        if(evol_mode == Core::BINARY) 
            angmom_index += (
                2 * __other_body.number_zones()
                + 
                (__primary 
                 ? 0
                 : (__other_body.number_zones()
                    -
                    __other_body.number_locked_zones()))
            );
        else angmom_index -= 3;
        assert(angmom_index <= derivatives.size());

        double surf_angmom_deriv = derivatives[angmom_index];

        if(__critical_below_iter != __critical_spins.end())
            stop_deriv[0] = derivative(surf_angmom_deriv,
                                       *__critical_below_iter);
        if(__critical_above_iter != __critical_spins.end())
            stop_deriv[1] = derivative(surf_angmom_deriv,
                                       *__critical_above_iter);
    }

    CriticalSpinCondition::CriticalSpinCondition(
        const DissipatingBody &body,
        const DissipatingBody &other_body,
        bool primary,
        unsigned zone_index,
        std::vector<double> critical_spins
    ) : 
        __critical_spins(critical_spins),
        __zone(body.zone(zone_index)),
        __critical_above_iter(std::lower_bound(__critical_spins.begin(),
                                               __critical_spins.end(),
                                               __zone.spin_frequency())),
        __critical_below_iter(
            __critical_above_iter == __critical_spins.begin()
            ? __critical_spins.end()
            : __critical_above_iter - 1
        ),
        __body(body),
        __other_body(other_body),
        __primary(primary),
        __zone_index(zone_index)
    {
        if(__zone.spin_frequency() == *__critical_above_iter)
            throw Core::Error::BadFunctionArguments(
                "Starting evolution from exactly a critical spin frequency "
                "is not currently supported."
            );
        set_num_subconditions();
    }

    std::valarray<double> CriticalSpinCondition::operator()(
        Core::EvolModeType evol_mode,
        const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
        assert(evol_mode != Core::LOCKED_SURFACE_SPIN);
        assert(__primary || evol_mode == Core::BINARY);

        stop_deriv.resize(Core::NaN, __num_subconditions);


        if(__zone.locked())
            fill_locked_derivs(evol_mode, orbit, derivatives, stop_deriv);
        else
            fill_unlocked_derivs(evol_mode, orbit, derivatives, stop_deriv);

        std::valarray<double> result(__num_subconditions);

        unsigned subcond = 0;
        if(__critical_below_iter != __critical_spins.end())
            result[subcond++] = (
                (__zone.spin_frequency() - *__critical_below_iter)
                /
                *__critical_below_iter
            );
        if(__critical_above_iter != __critical_spins.end())
        result[subcond] = (
            (__zone.spin_frequency() - *__critical_above_iter)
            /
            *__critical_above_iter
        );
        return result;
    }

    void CriticalSpinCondition::reached(short deriv_sign, unsigned index)
    {
        assert(index < __num_subconditions);
        if(__critical_below_iter != __critical_spins.end() && index == 0) {
            assert(deriv_sign == -1);
            --__critical_above_iter;
            assert(__critical_above_iter == __critical_below_iter);
            if(__critical_below_iter == __critical_spins.begin())
                __critical_below_iter = __critical_spins.end();
            else 
                --__critical_below_iter;
        } else {
            assert(deriv_sign == 1);
            assert(__critical_above_iter != __critical_spins.end());
            if(__critical_below_iter == __critical_spins.end())
                __critical_below_iter = __critical_spins.begin();
            else
                ++__critical_below_iter;
            assert(__critical_above_iter == __critical_below_iter);
            ++__critical_above_iter;
        }
        set_num_subconditions();
    }

}//End Evolve namespace
