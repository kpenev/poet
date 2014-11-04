#ifndef __BREAK_LOCK_CONDITION_H
#define __BREAK_LOCK_CONDITION_H

#include "StoppingCondition.h"

class BinarySystem;

///\brief Satisfied when the maximum tidal torque that the planet can exert
///on the star is no longer sufficient to keep the lock.
///
///\ingroup OrbitSolver_group
class BreakLockCondition : public StoppingCondition {
private:
	///The binary system this condition is attached to.
	BinarySystem &__system;

	///The index within the list of locked zones of the checked zone.
	unsigned __locked_zone_index;
public: 
	///Create a condition monitoring for a lock breaking.
	BreakLockCondition(
			///The binary system this locking condition is attached to
			BinarySystem &system,

			///Index within the list of locked zones of the checked zone.
			unsigned locked_zone_index) :
		__system(system), __locked_zone_index(locked_zone_index) {}

	///\brief How far away from breaking the lock the system is.
	///
	///Two values are calculated: the fraction of the above the
	///lock terms which must be included in order to maintain the lock and
	///that minus 1.
	///
	///This way, if the lock is broken with a positive/negative value, the
	///future evolution is expected to proceed with positive/negative forcing
	///frequency.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	///
	///The evolution mode must be LOCKED_TO_PLANET.
	std::valarray<double> operator()(EvolModeType evol_mode, 
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			std::valarray<double> &stop_deriv) const;

	///The number of subconditions in the current condition.
	virtual size_t num_subconditions() const {return 2;}

	///Identify this as a BREAK_LOCK condition.
	StoppingConditionType type(unsigned =0) const {return BREAK_LOCK;}

	///See StoppingCondition::reached().
	void reached(short deriv_sign, unsigned index=0);

	///\brief See StoppingCondition::expected_crossing_deriv_sign().
	virtual short expected_crossing_deriv_sign(
			///Which sub-condition.
			unsigned index=0)
	{
#ifdef DEBUG
		assert(index<=1);
#endif
		return (index==0 ? -1 : 1);
	}
};

#endif
