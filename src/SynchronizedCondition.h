#ifndef __SYNCHRONIZED_CONDITION_H
#define __SYNCHRONIZED_CONDITION_H

#include "StoppingCondition.h"

class DissipatingZone;
class BinarySystem;

///\brief Satisfied when some multiples of the orbit and stellar rotation are
///synchronized.
///
///\ingroup OrbitSolver_group
class SynchronizedCondition : public StoppingCondition{
private:
	///The mutiplier in front of the orbital frequency in the lock.
	int __orbital_freq_mult,

		///The multiplier in front of the spin frequency in the lock.
		__spin_freq_mult;

	///Which body's spin is checked for locking.
	bool __primary;

	///Which zone is checked for locking.
	unsigned __zone_index;

	///The sign of the derivative expected at zerocrossing.
	short __expected_crossing_deriv_sign;

	///The zone whose spin is monitored.
	const DissipatingZone &__zone;

	///The binary system this locking condition is attached to.
	BinarySystem &__system;
public:
	///Create the synchronization condition for the given planet.
	SynchronizedCondition(
			///The mutiplier in front of the orbital frequency in the lock.
			int orbital_freq_mult,

			///The multiplier in front of the spin frequency in the lock.
			int spin_freq_mult,

			///The sign the first derivative should have if a crossing
			///occurs.
			short deriv_sign,

			///Which body's spin is checked for locking.
			bool primary,
			
			///Which zone's spin is checked for locking.
			unsigned zone_index,
			
			///The binary system this locking condition is attached to
			BinarySystem &system);

	///\brief Returns the difference between the orbital and multiplier
	///scaled stellar spin angular velocities divided by the orbital angular
	///velocity.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	///
	///The evolution mode must be FAST_PLANET, SLOW_PLANET or LOCKED_TO_DISK
	///or NO_PLANET.
	std::valarray<double> operator()(EvolModeType evol_mode,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			std::valarray<double> &stop_deriv) const;

	///Identify this as a SYNCHRONIZED condition.
	StoppingConditionType type(unsigned =0) const {return SYNCHRONIZED;}

	///The multiplier in front of the orbital frequency in the lock.
	int orbital_frequency_multiplier() const {return __orbital_freq_mult;}

	///The multiplier in front of the spin frequency in the lock.
	int spin_frequency_multiplier() const {return __spin_freq_mult;}

	///Which body's spin is checked for locking.
//	short body_index() const {return __body_index;}

	///See StoppingCondition::reached().
	void reached(short deriv_sign, unsigned index=0);

	///\brief The expected sign of the derivative at the next zero-crossing.
	///
	///Zero if unknown.
	virtual short expected_crossing_deriv_sign(
			///Which sub-condition.
			unsigned
#ifdef DEBUG
		   	index
#endif
			=0)
	{
#ifdef DEBUG
		assert(index==0);
#endif
		return __expected_crossing_deriv_sign;
	}
};


#endif
