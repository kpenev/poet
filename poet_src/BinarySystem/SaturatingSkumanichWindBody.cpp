#include "SaturatingSkumanichWindBody.h"
#include "WindSaturationCondition.h"

double SaturatingSkumanichWindBody::angular_momentum_loss(
		Dissipation::Derivative deriv) const
{
	double result=__wind_strength*std::sqrt(radius()/mass())
				  *
				  (__saturated
				   ? spin_frequency()*std::pow(__saturation_freq, 2)
				   : std::pow(spin_frequency(), 3));
	double freq_power=(__saturated ? 1.0 : 3.0);
	switch(deriv) {
		case Dissipation::NO_DERIV : return result;
		case Dissipation::SPIN_FREQUENCY : return freq_power*result
										   		  /spin_frequency();
		case Dissipation::RADIUS : return result/(2.0*radius());
		case Dissipation::MOMENT_OF_INERTIA :
								   return -freq_power*result/
									   	  zone(0).moment_of_inertia();
		case Dissipation::SPIN_ANGMOM : return freq_power*result
											   /zone(0).angular_momentum();
		default : return 0;
	}
}

CombinedStoppingCondition *SaturatingSkumanichWindBody::stopping_conditions(
		BinarySystem &system, bool primary)
{
#ifdef DEBUG
	if(primary) assert(this==&(system.primary()));
	else assert(this==&(system.secondary()));
#endif
	CombinedStoppingCondition *result=new CombinedStoppingCondition();
	if(system.evolution_mode()!=LOCKED_SURFACE_SPIN)
		(*result)|=new WindSaturationCondition(
				*this, (primary ? system.secondary() : system.primary()),
				primary);
	(*result)|=DissipatingBody::stopping_conditions(system, primary);
	return result;
}
