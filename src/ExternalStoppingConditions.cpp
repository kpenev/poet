#include "ExternalStoppingConditions.h"

std::valarray<double> RotFastCondition::operator()(EvolModeType,
		const std::valarray<double> &, const std::valarray<double> &,
		std::valarray<double> &stop_deriv) const
{
	if(!std::isfinite(__spin_thres)) return std::valarray<double>(-1, 1);
	double spin_freq=__zone.spin_frequency();
	stop_deriv.resize(1, NaN);
	return std::valarray<double>((spin_freq-__spin_thres)/__spin_thres, 1);
}
