#ifndef __POWERLAW_PHASE_LAG_ZONE_H
#define __POWERLAW_PHASE_LAG_ZONE_H

/**\file
 *
 * \brief Declares the class that provides the phase lag function to
 * DissipatingZone objects.
 *
 * \ingroup StellarSystem_group
 */

#include "DissipatingZone.h"

///A class that only defines the phase lag function for zones.
class PowerlawPhaseLagZone : virtual public DissipatingZone {
private:
	///The phase lag at a forcing frequency of one day.
	double __phase_lag_one_day,

		   ///The powerlaw index of the phase lag period dependence.
		   __phase_lag_powerlaw_index,
           
           ///The maximum phase lag achieved for any frequency
           __max_phase_lag;
public:
	///Create the zone with the phase lag scaling.
	PowerlawPhaseLagZone(
			///The phase lag at a forcing frequency of 1 day.
			double phase_lag_one_day = 0,

			///The powerlaw index of the phase lag period dependence.
			double phase_lag_powerlaw_index = 0,
            
            ///The maximum phase lag achieved at any frequency.
            double max_phase_lag = 0) :
		__phase_lag_one_day(phase_lag_one_day),
		__phase_lag_powerlaw_index(phase_lag_powerlaw_index),
        __max_phase_lag(max_phase_lag) {}

	///Set the modified phase lag at a forcing frequency of one day.
	void set_phase_lag_one_day(double lag_1day)
	{__phase_lag_one_day = lag_1day;}

	///Set the modified pase lag period dependence powerlaw index.
	void set_phase_lag_powerlaw_index(double powerlaw_index)
	{__phase_lag_powerlaw_index = powerlaw_index;}

	///Set the maximum phase lag achieved for any frequency.
	void set_max_phase_lag(double max_lag)
	{__max_phase_lag = max_lag;}

	///\brief Should return the tidal phase lag time the love number for the
	///given tidal term (or one of its derivatives).
	///
	///In case the forcing frequency is exactly zero, it should return the
	///phase lag for the case of the spin frequency approaching the term from
	///below. The lag for spin frequency approaching from above should be
	///written to above_lock_value. If the forcing frequency is non-zero, 
	///leave above_lock_value untouched.
	virtual double modified_phase_lag(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,
			
			///The current forcing frequency in rad/day.
			double forcing_frequency,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::QuantityEntry entry,

			///If the lag of a locked term is calculated this should be set
			///to the lag assuming the spin frequency is just above the lock.
			///Otherwise, leave untouched.
			double &above_lock_value
    ) const;

	///\brief Should return the corresponding component of the love
	///coefficient (Lai 2012 Equation 24).
	virtual double love_coefficient(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::QuantityEntry
    ) const
    {return 0;}
};

#endif
