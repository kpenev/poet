#include "PowerlawPhaseLagZone.h"

double PowerlawPhaseLagZone::modified_phase_lag(
		int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
		double forcing_frequency,
        Dissipation::QuantityEntry entry,
		double &above_lock_value) const
{
    if(entry == Dissipation::AGE || entry == Dissipation::EXPANSION_ERROR)
        return 0;

    double result = (__phase_lag_one_day
                     *
                     std::pow(std::abs(forcing_frequency) / (2.0 * M_PI),
                              __phase_lag_powerlaw_index));

    if(result > __max_phase_lag) {
        if(entry == Dissipation::NO_DERIV) result = __max_phase_lag;
        else result = 0.0;
    } else if(entry != Dissipation::NO_DERIV) {
        result *= __phase_lag_powerlaw_index / forcing_frequency;
        if(entry == Dissipation::SPIN_FREQUENCY) 
            result *= -spin_frequency_multiplier;
        if(entry == Dissipation::ORBITAL_FREQUENCY)
            result *= orbital_frequency_multiplier;
    }

    if(forcing_frequency == 0) above_lock_value = -result;
    return (forcing_frequency > 0 ? result : -result);
}
