#include "SecondaryDeathCondition.h"

std::valarray<double> SecondaryDeathCondition::operator()(
		EvolModeType
#ifdef DEBUG
		evol_mode
#endif
		,
		const std::valarray<double> &orbit, 
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode==BINARY);
	assert(orbit.size()==1 + 3*__system.number_zones() -
						 __system.number_locked_zones());
	assert(orbit.size()==derivatives.size());
#endif
	double min_semimajor=__system.minimum_semimajor(),
		   semimajor=__system.semimajor(),
		   dsemimajor_dt=derivatives[0];
	if(__system.number_locked_zones()==0) 
		dsemimajor_dt*=semimajor/(6.5*orbit[0]);
	stop_deriv.resize(1,(dsemimajor_dt*min_semimajor
						  -semimajor*__system.minimum_semimajor(true))
						/std::pow(min_semimajor, 2));
	return std::valarray<double>(
			(semimajor-min_semimajor)/min_semimajor, 1);
}

void reached(short deriv_sign, unsigned index=0)
{
	__system.secondary_died();
}
