#include "ExternalStoppingConditions.h"

std::valarray<double> RotFastCondition::operator()(
		const std::valarray<double> &orbit,
		//Accepts but never uses derivatives
		const std::valarray<double> &,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	if(!std::isfinite(spin_thres)) return std::valarray<double>(-1, 1);
	double spin_freq=__zone.spin_frequency();
	stop_deriv.resize(1, NaN);
	return std::valarray<double>((spin_freq-spin_thres)/spin_thres, 1);
}
