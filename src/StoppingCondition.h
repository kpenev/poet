/**\file
 *
 * \brief Defines the various stopping conditions needed by OrbitSolver.
 *
 * \ingroup OrbitSolver_group
 */

#ifndef __STOPPING_CONDITION_H
#define __STOPPING_CONDITION_H

#include "Common.h"

#include <cassert>

///The reasons for stopping the evolution currently supported.
enum StoppingConditionType {
	///No reason to stop.
	NO_STOP,
	
	///Synchroneity between the stellar spin and the orbit was reached
	///(from either direction).
	SYNCHRONIZED,
	
	///The spin-orbit lock can no longer be maintaned.
	BREAK_LOCK,
	
	///The planet died.
	PLANET_DEATH,
	
	///The stellar spin frequency reached the wind saturation frequency
	///(from either direction).
	WIND_SATURATION,
	
	///The stellar rotation reached some critical value (from either
	///direction).
	EXTERNAL};

///More civilized output for StoppingConditionType variables.
std::ostream &operator<<(std::ostream &os,
		const StoppingConditionType &stop_cond_type);

///\brief A base class for all stopping conditions.
///
///A more detailed description of the stopping condition mechanism and why it
///is necessary are given in the \ref stop_cond section.
///
///\ingroup OrbitSolver_group
class StoppingCondition {
private:
	///\brief The sign of the first derivative at the next zero-crossing of
	///this condition.
	short __expected_crossing_deriv_sign;
public:
	///Create a generic stopping condition.
	StoppingCondition() : __expected_crossing_deriv_sign(0) {}

	///\brief The values of quantities which should cross zero when the
	///condition(s) is(are) satisfied.
	///
	///The input stellar system must already have its age set.
	virtual std::valarray<double> operator()(
			///The evolution mode for which the orbit and derivatives are
			///given. For some conditions some EvolModeType values will not
			///make sense and will result in an exception.
			EvolModeType evol_mode,

			///The variables which are currently being evolved. The content
			///depends on the evol_mode argument.
			const std::valarray<double> &orbit, 

			///The rate of change of the entries in orbit per the relevant
			///system of differential equations.
			const std::valarray<double> &derivatives,

			///On output contains the rate of change of the stopping
			///sub-conditions if known, or NaN if not.
			std::valarray<double> &stop_deriv) const=0;

	///The number of subconditions in the current condition.
	virtual size_t num_subconditions() const {return 1;}

	///What event is the index-th stopping sub-condition associated with.
	virtual StoppingConditionType type(unsigned index=0) const=0;

	///\brief Called when a stopping condition has been reached by the
	///evolution.
	///
	///Should perform any necessary changes to the further evolution (e.g.
	///switch wind saturatino state or check if a spin-orbit lock can be held
	///and lock it if so etc.).
	virtual void reached(
			///The sign of the first derivative when the condition was
			///reached.
			short deriv_sign,

			///The sub-condition reached for composite conditions.
			unsigned 
#ifdef DEBUG
			index
#endif
			=0)
	{
#ifdef DEBUG
		assert(index==0);
#endif
		__expected_crossing_deriv_sign=-deriv_sign;
	}

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

    virtual ~StoppingCondition() {}
};

///\brief A stopping condition that is never satisfied.
///
///\ingroup OrbitSolver_group
class NoStopCondition : public StoppingCondition {
public:
	///A single value is reternud that is always 1.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	std::valarray<double> operator()(
			EvolModeType, const std::valarray<double> &, 
			const std::valarray<double> &, 
			std::valarray<double> &stop_deriv) const
	{stop_deriv.resize(1, NaN); return std::valarray<double>(1, 1);}

	///Identifies the condition as NO_STOP.
	StoppingConditionType type(unsigned =0) const {return NO_STOP;}

	///See StoppingCondition::reached().
	void reached(short, unsigned=0)
	{
#ifdef DEBUG
		assert(false);
#endif
	}
};

///\brief A base class for all external stopping conditions.
///
///\ingroup OrbitSolver_group
class ExternalStoppingCondition : public StoppingCondition {
public:
	///Identify this as an EXTERNAL condition.
	StoppingConditionType type(unsigned =0) const {return EXTERNAL;}
};

#endif
