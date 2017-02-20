#include "BrokenPowerlawPhaseLagZone.h"

namespace Evolve {

    BrokenPowerlawPhaseLagZone::BrokenPowerlawPhaseLagZone(
        std::vector<double> tidal_frequency_breaks,
        std::vector<double> spin_frequency_breaks,
        std::vector<double> tidal_frequency_powers,
        std::vector<double> spin_frequency_powers,
        double reference_phase_lag
    ) :
        __tidal_frequency_breaks(tidal_frequency_breaks),
        __spin_frequency_breaks(spin_frequency_breaks),
        __tidal_frequency_powers(tidal_frequency_powers),
        __spin_frequency_powers(spin_frequency_powers),
        __break_phase_lags(spin_frequency_breaks.size()
                           *
                           tidal_frequency_breaks.size())
    {
        assert(tidal_frequency_powers.size()
               == 
               tidal_frequency_breaks.size() + 1);
        assert(spin_frequency_powers.size()
               ==
               spin_frequency_breaks.size() + 1);
        unsigned break_lag_i = 0;
        for(
            unsigned spin_break_i = 0;
            spin_break_i < spin_frequency_breaks.size();
            ++spin_break_i
        ) {
            if(break_lag_i == 0)
                __break_phase_lags[break_lag_i] = reference_phase_lag;
            else 
                __break_phase_lags[break_lag_i] = (
                    __break_phase_lags[
                        break_lag_i - tidal_frequency_breaks.size()
                    ]
                    *
                    std::pow(
                        (
                            spin_frequency_breaks[spin_break_i]
                            /
                            spin_frequency_breaks[spin_break_i - 1]
                        ),
                        spin_frequency_powers[spin_break_i]
                    )
                );
            ++break_lag_i;
            for(
                unsigned tidal_break_i = 1;
                tidal_break_i < tidal_frequency_breaks.size();
                ++tidal_break_i
            ) {
                __break_phase_lags[break_lag_i] = (
                    __break_phase_lags[break_lag_i - 1]
                    *
                    std::pow(
                        (
                            tidal_frequency_breaks[tidal_break_i]
                            /
                            tidal_frequency_breaks[tidal_break_i - 1]
                        ),
                        tidal_frequency_powers[spin_break_i]
                    )
                );
                ++break_lag_i;
            }
        }
    }//End BrokenPowerlawPhaseLagZone::BrokenPowerlawPhaseLagZone definition.

    double BrokenPowerlawPhaseLagZone::modified_phase_lag(
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        double forcing_frequency,
        Dissipation::Derivative deriv,
        double &above_lock_value
    ) const
    {
        if(deriv == Dissipation::AGE) return 0;

        std::vector<double>::const_iterator 
            tidal_break_iter = std::lower_bound(
                __tidal_frequency_breaks.begin(),
                __tidal_frequency_breaks.end() - 1,
                forcing_frequency
            ),
            spin_break_iter = std::lower_bound(
                __spin_frequency_breaks.begin(),
                __spin_frequency_breaks.end() - 1,
                spin_frequency()
            );
        unsigned tidal_index = (tidal_break_iter
                                -
                                __tidal_frequency_breaks.begin()),
                 spin_index = (spin_break_iter
                               -
                               __spin_frequency_breaks.begin()),
                 lag_index = (spin_index * __tidal_frequency_breaks.size()
                              +
                              tidal_index);
        double 
            tidal_power = __tidal_frequency_powers[tidal_index],
            spin_power = __spin_frequency_powers[spin_index],
            result = (
                __break_phase_lags[lag_index]
                *
                std::pow(forcing_frequency / *tidal_break_iter, tidal_power)
                *
                std::pow(spin_frequency() / *spin_break_iter, spin_power)
            );
        switch(deriv) {
            case Dissipation::SPIN_FREQUENCY :
                result *= (
                    spin_power / spin_frequency()
                    -
                    spin_frequency_multiplier
                    *
                    tidal_power
                    /
                    forcing_frequency
                );
                break;
            case Dissipation::ORBITAL_FREQUENCY :
                result *= (
                    orbital_frequency_multiplier
                    * 
                    tidal_power
                    /
                    forcing_frequency
                );
                break;
            default :
                assert(deriv == Dissipation::NO_DERIV);
        }

        if(forcing_frequency == 0) above_lock_value  = -result;
        return (forcing_frequency > 0 ? result : -result);

    }//End BrokenPowerlawPhaseLagZone::modified_phase_lag definition.

    CombinedStoppingCondition *
        BrokenPowerlawPhaseLagZone::stopping_conditions(
            ///The system being evolved.
            BinarySystem &system, 

            ///Is the body this zone is part of, the primary in the system.
            bool primary,

            ///The index of the zone in the body.
            unsigned zone_index
        )
    {
        CombinedStoppingCondition  *result = 
            DissipatingZone::stopping_conditions(system,
                                                 primary,
                                                 zone_index);
    }//End BrokenPowerlawPhaseLagZone::stopping_conditions definition.

} //End BinarySystem namespace.
