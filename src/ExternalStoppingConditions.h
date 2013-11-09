#ifndef __EXTERNAL_STOPPING_CONDITIONS_H
#define __EXTERNAL_STOPPING_CONDITIONS_H

#include "StoppingConditions.h"

/**\file
 *
 * \brief Users can define any stopping condition they wish the evolution to
 * search for in this file.
 *
 * \ingroup OrbitSolver_group
 */

///\brief Satisfied when the star is rotating faster than a threshold.
///
///\ingroup OrbitSolver_group
class RotFastCondition : public ExternalStoppingCondition {
private:
	///The spin threshold in rad/day.
	double spin_thres;
public:
	///Create a condition tied to the given threshold in rad/day.
	RotFastCondition(double spin_thres): spin_thres(spin_thres) {}

	///\brief Returns the difference between the convective zone spin and
	///the threshold divided by the latter.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;
};

///\brief How to construct the external stopping condition(s).
///
///Leave undefined if no external conditions need to be tracked.
//#define EXTERNAL_CONDITION RotFastCondition(4.0*M_PI)

#endif
