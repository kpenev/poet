#include "WindSaturationCondition.h"

std::valarray<double> WindSaturationCondition::operator()(
		EvolModeType evol_mode,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode!=LOCKED_SURFACE_SPIN);
	if(evole_mode!=BINARY) assert(__primary);
#endif 
	unsigned num_zones=__body.number_zones();
	if(evol_mode==BINARY) num_zones+=__other_body.number_zones();
	double wsurf=__body.spin_frequency(),
		   surf_angmom_deriv=
			   stop_deriv[1 + 2*num_zones
			   			  +
						  (__primary ? 0
						   			 : __other_body.number_zones()
									   -__other_body.number_locked_zones())
						 ];
	stop_deriv.resize(1,
		   	(surf_angmom_deriv - __body.zone(0).moment_of_inertia(1)*wsurf)
			/(__body.zone(0).moment_of_inertia()*wsat));
	if(std::isinf(wsat)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>((wsurf-wsat)/wsat, 1);
}
