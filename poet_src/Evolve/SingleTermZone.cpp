#define BUILDING_LIBRARY

#include "SingleTermZone.h"

namespace Evolve {

    void SingleTermZone::setup(int orbital_frequency_multiplier,
                               int spin_frequency_multiplier,
                               double phase_lag)
    {
        __orbital_frequency_multiplier = orbital_frequency_multiplier;
        __spin_frequency_multiplier = spin_frequency_multiplier;
        __phase_lag = phase_lag;
    }

    double SingleTermZone::modified_phase_lag(
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        double forcing_frequency,
        Dissipation::QuantityEntry entry,
        double &above_lock_value
    ) const
    {
        if(
            !(
                (
                    orbital_frequency_multiplier == __orbital_frequency_multiplier
                    &&
                    spin_frequency_multiplier == __spin_frequency_multiplier
                )
                ||
                (
                    orbital_frequency_multiplier == -__orbital_frequency_multiplier
                    &&
                    spin_frequency_multiplier == -__spin_frequency_multiplier
                )
            )
            ||
            entry != Dissipation::NO_DERIV
        )
            return 0.0;
        if(forcing_frequency == 0) {
            if(spin_frequency_multiplier >= 0) {
                above_lock_value = -__phase_lag;
                return __phase_lag;
            } else {
                above_lock_value = __phase_lag;
                return -__phase_lag;
            }

        }

        return (forcing_frequency > 0 ? __phase_lag : -__phase_lag);

    }

} //End Evolve namespace
