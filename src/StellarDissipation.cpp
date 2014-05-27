#include "StellarDissipation.h"

double Star::modified_phase_lag(int m,  double forcing_frequency,
		const SpinOrbitLockInfo &lock, PhaseLag::Derivative derivative) const
{
	if(derivative) return 0;
	if(forcing_frequency==0) {
		assert(false);
#ifdef DEBUG
		assert(lock.lock_direction()!=0);
#endif
		return lock.lock_direction()*(m>0 ? 1 : -1)*__mod_phase_lag*forcing_frequency;
	} else {
		if(std::abs(forcing_frequency)<2.0*spin())
			return __mod_phase_lag_inertial*forcing_frequency;
		else return __mod_phase_lag*forcing_frequency;
	}
}
