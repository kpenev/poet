#include "StellarDissipation.h"

double Star::modified_phase_lag(int m,  double forcing_frequency,
		const SpinOrbitLockInfo &lock, PhaseLag::Derivative derivative) const
{
	return (forcing_frequency>0 ? 1 : -1)*__mod_phase_lag;
}
