#include "LagSpinBreakCondition.h"
#include "BrokenPowerlawPhaseLagZone.h"

namespace Evolve {

    void LagSpinBreakCondition::set_num_subconditions()
    {

        if(
            __powerlaw_index == 0
            ||
            __powerlaw_index == __zone.__spin_frequency_breaks.size()
        )
            __num_subconditions = 1;
        else
            __num_subconditions = 2;
    }

    double LagSpinBreakCondition::derivative(double surf_angmom_deriv,
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
                __zone.spin_frequency()
            )
            /
            (__zone.moment_of_inertia() * wcritical)
        );
    }

    void LagSpinBreakCondition::fill_locked_derivs(
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
        if(__powerlaw_index > 0)
            stop_deriv[subcond++] = (
                deriv_factor
                /
                __zone.__spin_frequency_breaks[__powerlaw_index - 1]
            );
        if(__powerlaw_index < __zone.__spin_frequency_breaks.size())
            stop_deriv[subcond] = (
                deriv_factor
                /
                __zone.__spin_frequency_breaks[__powerlaw_index]
            );
    }

    void LagSpinBreakCondition::fill_unlocked_derivs(
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

        if(__powerlaw_index > 0)
            stop_deriv[0] = derivative(
                surf_angmom_deriv,
                __zone.__spin_frequency_breaks[__powerlaw_index - 1]
            );
        if(__powerlaw_index < __zone.__spin_frequency_breaks.size())
            stop_deriv[1] = derivative(
                surf_angmom_deriv,
                __zone.__spin_frequency_breaks[__powerlaw_index]
            );
    }

    LagSpinBreakCondition::LagSpinBreakCondition(
        BrokenPowerlawPhaseLagZone &zone,
        const DissipatingBody &body,
        const DissipatingBody &other_body,
        bool primary,
        unsigned zone_index
    ) : 
        __zone(zone),
        __body(body),
        __other_body(other_body),
        __primary(primary),
        __zone_index(zone_index),
        __powerlaw_index(__zone.__spin_index)
    {
        set_num_subconditions();
    }

    std::valarray<double> LagSpinBreakCondition::operator()(
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
        if(__powerlaw_index > 0) {
            double critical_frequency =
                __zone.__spin_frequency_breaks[__powerlaw_index - 1];
            result[subcond++] = (
                (__zone.spin_frequency() - critical_frequency)
                /
                critical_frequency
            );
        } if(__powerlaw_index < __zone.__spin_frequency_breaks.size()) {
            double critical_frequency =
                __zone.__spin_frequency_breaks[__powerlaw_index];
            result[subcond] = (
                (__zone.spin_frequency() - critical_frequency)
                /
                critical_frequency
            );
        }
        return result;
    }

    void LagSpinBreakCondition::reached(short deriv_sign, unsigned index)
    {
        assert(index < __num_subconditions);

        if(__powerlaw_index > 0 && index == 0) {
            assert(deriv_sign == -1);
            --__powerlaw_index;
            --__zone.__spin_index;
        } else {
            assert(deriv_sign == 1);
            assert(__powerlaw_index < __zone.__spin_frequency_breaks.size());
            ++__powerlaw_index;
            ++__zone.__spin_index;
        }

        set_num_subconditions();

    }

    short LagSpinBreakCondition::expected_crossing_deriv_sign(
        unsigned index
    ) const
    {
        if(index == 1 || __powerlaw_index == 0) {
            assert(__powerlaw_index
                   <
                   __zone.__spin_frequency_breaks.size());
            return 1;
        } else {
            return -1;
        }
    }

    std::string LagSpinBreakCondition::describe(int index) const
    {
        std::ostringstream description;
        description << (__primary ? "Primary" : "Secondary")
                    << " body, zone "
                    << __zone_index
                    << " spin frequency leaving the interval ";
        if(__powerlaw_index > 0 && index <= 0)
            description
                << __zone.__spin_frequency_breaks[__powerlaw_index - 1]
                << " < ";
        description << " |Wspin| ";
        if(
            __powerlaw_index < __zone.__spin_frequency_breaks.size()
            &&
            (
                index < 0
                ||
                index == __num_subconditions - 1
            )
        )
            description << " < "
                        << __zone.__spin_frequency_breaks[__powerlaw_index];
        description << std::endl;
        return description.str();
    }

}//End Evolve namespace
