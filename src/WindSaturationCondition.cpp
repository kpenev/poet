#include "WindSaturationCondition.h"

std::valarray<double> WindSaturationCondition::operator()(
		EvolModeType evol_mode,
		const std::valarray<double> &,
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode!=LOCKED_SURFACE_SPIN);
	if(evol_mode!=BINARY) assert(__primary);
#endif 
	unsigned angmom_index=1 + 2*__body.number_zones();
	if(evol_mode==BINARY) 
		angmom_index+=2*__other_body.number_zones()
					  +(__primary ? 0
							      : __other_body.number_zones()
								    -__other_body.number_locked_zones());
	else angmom_index-=3;
	assert(angmom_index>=0);
	assert(angmom_index<=derivatives.size());
	double wsurf=__body.spin_frequency(), surf_angmom_deriv,
		   result=(wsurf-__saturation_freq)/__saturation_freq;
	surf_angmom_deriv=derivatives[angmom_index];
	stop_deriv.resize(1,
		   	(surf_angmom_deriv - __body.zone(0).moment_of_inertia(1)*wsurf)
			/(__body.zone(0).moment_of_inertia()*__saturation_freq));
	if(std::isinf(__saturation_freq)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>(result, 1);
}
