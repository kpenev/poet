#include "TwoPhaseLagZone.h"

double TwoPhaseLagZone::modified_phase_lag(
		int orbital_frequency_multiplier, int spin_frequency_multiplier,
		double forcing_frequency, Dissipation::Derivative deriv,
		double &above_lock_value) const
{
	if(deriv!=Dissipation::NO_DERIV) return 0;
	if(std::abs(forcing_frequency)>2.0*std::abs(spin_frequency()))
		return (forcing_frequency>0 ? __equilibrium_modified_lag
									: -__equilibrium_modified_lag);
	else if(forcing_frequency==0) {
		above_lock_value=(spin_frequency_multiplier>0
						  ? -__inertial_modified_lag
						  : __inertial_modified_lag);
		return -above_lock_value;
	} else return (forcing_frequency>0 ? __inertial_modified_lag
									: -__inertial_modified_lag);
}
