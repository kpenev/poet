/**\file
 *
 * \brief Defines the various stopping conditions needed by OrbitSolver.
 *
 * \ingroup OrbitSolver_group
 */

#ifndef __STOPPING_CONDITIONS_H
#define __STOPPING_CONDITIONS_H

#include "Common.h"
#include "StellarSystem.h"

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
	///Not used for now.
	double __interpolation_range;
public:
	///Create a generic stopping condition.
	StoppingCondition(
			///Not used for now.
			double interp_range=Inf) :
		__interpolation_range(interp_range) {}

	///\brief The values of quantities which should cross zero when the
	///condition(s) is(are) satisfied.
	virtual std::valarray<double> operator()(
			///System age in Gyr.
			double age,
			
			///The variables which are currently being evolved. The content
			///depends on the evol_mode argument.
			const std::valarray<double> &orbit, 

			///The rate of change of the entries in orbit per the relevant
			///system of differential equations.
			const std::valarray<double> &derivatives,

			///The planet-star system being evolved.
			const StellarSystem &system,

			///On output contains the rate of change of the stopping
			///sub-conditions if known, or NaN if not.
			std::valarray<double> &stop_deriv,

			///The evolution mode for which the orbit and derivatives are
			///given. For some conditions some EvolModeType values will not
			///make sense and will result in an exception.
			EvolModeType evol_mode) const=0;

	///The number of subconditions in the current condition.
	virtual size_t num_subconditions() const {return 1;}

	///What event is the index-th stopping sub-condition associated with.
	virtual StoppingConditionType type(unsigned index=0) const=0;

	///\brief Not used for now.
	///
	///The idea was that if a zero crossing is found only points within
	///__interpolation_range*(t_after-t_before) should be considered, where
	///t_after is tha age of the first point found after the zero crossing
	///and t_before is the age of the last point before the zero crossing.
	///In addition the interpolation should only include monotonic points, in
	///addition if derivatives are provided, no points other than the ones at
	///t_before and t_after should be considered.
	virtual double interpolation_range(unsigned =0) const
	{return __interpolation_range;}

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
	std::valarray<double> operator()(double,
			const std::valarray<double>&,
			const std::valarray<double>&,
			const StellarSystem&,
			std::valarray<double> &stop_deriv,
			EvolModeType) const
	{stop_deriv.resize(1, NaN); return std::valarray<double>(1, 1);}

	///Identifies the condition as NO_STOP.
	StoppingConditionType type(unsigned =0) const {return NO_STOP;}
};

///\brief Satisfied when the orbit and stellar rotation are synchronized.
///
///\ingroup OrbitSolver_group
class SynchronizedCondition : public StoppingCondition{
private:
	///The planet whose orbit is checked for synchronization with the star.
	const Planet *planet;
	
	///The semimajor axis at which the planet first appears.
	double __initial_semimajor;
public:
	///Create the synchronization condition for the given planet.
	SynchronizedCondition(const Planet *p, double initial_semimajor) :
		planet(p), __initial_semimajor(initial_semimajor) {}

	///\brief Returns the difference between the orbital and stellar spin
	///angular velocities divided by the orbital angular velocity.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	///
	///The evolution mode must be FAST_PLANET or SLOW_PLANET
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Identify this as a SYNCHRONIZED condition.
	StoppingConditionType type(unsigned =0) const {return SYNCHRONIZED;}
};

///\brief Satisfied when the maximum tidal torque that the planet can exert
///on the star is no longer sufficient to keep the lock.
///
///\ingroup OrbitSolver_group
class BreakLockCondition : public StoppingCondition {
public: 
	///\brief The difference between the maximum evolution of the semimajor
	///axis and the evolution required to keep the star synchronized divided
	///by the latter.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	///
	///The evolution mode must be LOCKED_TO_PLANET.
	std::valarray<double> operator()(
			double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Identify this as a BREAK_LOCK condition.
	StoppingConditionType type(unsigned =0) const {return BREAK_LOCK;}
};

///\brief Satisfied when the planet enters below either the roche sphere or
///the stellar photosphere.
///
///\ingroup OrbitSolver_group
class PlanetDeathCondition : public StoppingCondition {
public:
	///\brief The difference between the semimajor axis and the larger of
	///the roche radius and the stellar radius divided by the latter.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	///
	///The evolution mode must be FAST_PLANET, LOCKED_TO_PLANET or
	///SLOW_PLANET.
	std::valarray<double> operator()(
			double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Identify this as a PLANET_DEATH condition.
	StoppingConditionType type(unsigned =0) const {return PLANET_DEATH;}
};

///\brief Satisfied when the stellar convective zone is spinning at exactly
///the wind saturation frequency.
///
///\ingroup OrbitSolver_group
class WindSaturationCondition : public StoppingCondition {
public:
	///\brief The difference between the convective and wind saturation
	///angular velocities divided by the latter.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Identify this as a WIND_SATURATION condition.
	StoppingConditionType type(unsigned =0) const
	{return WIND_SATURATION;}
};

///\brief A base class for all external stopping conditions.
///
///\ingroup OrbitSolver_group
class ExternalStoppingCondition : public StoppingCondition {
public:
	///Identify this as an EXTERNAL condition.
	StoppingConditionType type(unsigned =0) const {return EXTERNAL;}
};

///\brief A class combining the the outputs of multiple stopping conditions.
///
///\ingroup OrbitSolver_group
class CombinedStoppingCondition : public StoppingCondition {
private:
	///The conditions that are to be combined
	std::vector<const StoppingCondition *> __sub_conditions;

	///\brief The number of subconditinos included, allowing for
	///subconditions with multiple entries.
	unsigned __num_subconditions;

	///Whether to delete the sub-conditions then *this is destroyed.
	bool __delete_subcond;

	///\brief The types of the subconditinos, including subconditions of
	///subconditions.
	std::vector<StoppingConditionType> __types;

	///\brief The interpolation ranges of the subconditinos, including
	///subconditions of subconditions.
	std::vector<double> __interpolation_ranges;
public:
	///Create an empty stopping condition (identical to NoStopCondition).
	CombinedStoppingCondition() :
		__sub_conditions(), __num_subconditions(0), __delete_subcond(true) {}

	///Adds the conditions in RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(
			const CombinedStoppingCondition &rhs);

	///Adds  RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(const StoppingCondition *rhs);

	///\brief Returns the values of all stopping sub_conditions.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Disables the destruction of the subconditions when *this is destroyed.
	void no_delete_subcond() {__delete_subcond=false;}

	virtual size_t num_subconditions() const {return __num_subconditions;}

	StoppingConditionType type(unsigned index=0) const
	{return __types[index];}

	double interpolation_range(unsigned index=0) const
	{return __interpolation_ranges[index];}

	///\brief Deletes all subconditions, unless no_delete_subcond has been
	///previously called.
	~CombinedStoppingCondition();
};

///Returns the spin frequency of the stellar convective zone in rad/day.
double convective_frequency(
		///The system age in Gyr.
		double age,
		
		///The stellar system.
		const StellarSystem &system, 

		///The present orbit as appropriate for the evol_mode parameter.
		const std::valarray<double> &orbit,
		
		///The present evolution mode.
		EvolModeType evol_mode);

#endif
