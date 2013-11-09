#include "ExternalStoppingConditions.h"

std::valarray<double> RotFastCondition::operator()(double age,
		const std::valarray<double> &orbit,
		//Accepts but never uses derivatives
		const std::valarray<double> &,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	if(!std::isfinite(spin_thres)) return std::valarray<double>(-1, 1);
	double wconv=convective_frequency(age, system, orbit, evol_mode);
	stop_deriv.resize(1, NaN);
	return std::valarray<double>((wconv-spin_thres)/spin_thres, 1);
}
