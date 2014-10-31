#include "SaturatingSkumanichWindBody.h"

double SaturatingSkumanichWindBody::angular_momentum_loss(
		Dissipation::Derivative deriv=Dissipation::NO_DERIV) const
{
	return __wind_strength*std::sqrt(radius()/mass())
		   *
		   (__saturated ? spin_frequency()*std::pow(__saturation_freq, 2)
						: std::pow(spin_frequency(), 3));
}

virtual CombinedStoppingCondition *
SaturatingSkumanichWindBody::stopping_conditions(BinarySystem &system, 
		bool primary)
{
	CombinedStoppingCondition *result=new CombinedStoppingCondition();
	(*result)|=WindSaturationCondition(
			*this, (primary ? system.primary() : system.secondary()),
			primary);
	(*result)|=DissipatingBody::stopping_conditions(system, primary);
}
