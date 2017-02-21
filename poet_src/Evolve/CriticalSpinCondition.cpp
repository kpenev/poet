#include "CriticalSpin.h"

namespace Evolve {

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
        Core::EvolModeType evol_mode,
        const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
        double deriv_factor = (-1.5 * __zone.spin_frequency()
                               *
                               derivatives[0] / orbit[0]);
        if(__critical_below_iter != __critical_spins.end())
            stop_deriv[0] = deriv_factor * (*critical_below_iter);
        if(__critical_above_iter != __critical_spins.end())
            stop_deriv[0] = deriv_factor * (*critical_above_iter);
    }

    void CriticalSpinCondition::fill_unlocked_derivs(
        Core::EvolModeType evol_mode,
        const std::valarray<double> &orbit,
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
        stop_deriv.resize(Core::NaN, 2);
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
        bool primary
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
        if(__zone.spin_frequency() == *critical_above_iter)
            throw Core::Error::BadFunctionArguments(
                "Starting evolution from exactly a critical spin frequency "
                "is not currently supported."
            )
    }

    std::valarray<double> CriticalSpinCondition::operator()(
        Core::EvolModeType evol_mode,
        const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
        assert(evol_mode != LOCKED_SURFACE_SPIN);
        assert(__primary || evol_mode == Core::BINARY);

        if(__zone.locked())
            fill_locked_derivs(evol_mode, orbit, derivatives, stop_deriv);
        else
            fill_unlocked_derivs(evol_mode, orbit, derivatives, stop_deriv);

        std::valarray<double> result(2);
        result[0] = (__critical_below_iter == __critical_spins.end()
                     ? Core::Inf
                     : ((__zone.spin_frequency() - *critical_below_iter)
                        /
                        *critical_below_iter));
        result[1] = (__critical_above_iter == __critical_spins.end()
                     ? Core::Inf
                     : ((__zone.spin_frequency() - *critical_above_iter)
                        /
                        *critical_above_iter));
        return result;
    }
