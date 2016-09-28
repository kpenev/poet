#ifndef __TWO_PHASE_LAG_ZONE_H
#define __TWO_PHASE_LAG_ZONE_H

/**\file
 *
 * \brief Declares the class that provides the phase lag function to
 * DissipatingZone objects.
 *
 * \ingroup StellarSystem_group
 */

#include "DissipatingZone.h"

///A class that only defines the phase lag function for zones.
class TwoPhaseLagZone : virtual public DissipatingZone {
private:
	///The modified pase lag outside the inertia wave frequency range.
	double __equilibrium_modified_lag,

		   ///The modified phase lag in the inertial wave frequency range.
		   __inertial_modified_lag;
public:
	///Create the zone with the given phase lags.
	TwoPhaseLagZone(
			///The modified pase lag outside the inertia wave frequency
			///range.
			double equilibrium_modified_lag=0,

			///The modified phase lag in the inertial wave frequency range.
			double inertial_modified_lag=0) :
		__equilibrium_modified_lag(equilibrium_modified_lag),
		__inertial_modified_lag(inertial_modified_lag) {}

	///Set the modified pase lag outside the inertia wave frequency range.
	void set_equilibrium_modified_lag(double lag)
	{__equilibrium_modified_lag=lag;}

	///Set the modified pase lag in the inertia wave frequency range.
	void set_inertial_modified_lag(double lag)
	{__inertial_modified_lag=lag;}

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
			Dissipation::Derivative deriv,

			///If the lag of a locked term is calculated this should be set
			///to the lag assuming the spin frequency is just above the lock.
			///Otherwise, leave untouched.
			double &above_lock_value) const;

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
			Dissipation::Derivative) const {return 0;}
};

#endif
