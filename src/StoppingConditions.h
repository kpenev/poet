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
#include "TidalDissipation.h"
#include "OrbitalExpressions.h"

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
	///
	///The input stellar system must already have its age set.
	virtual std::valarray<double> operator()(
			///The variables which are currently being evolved. The content
			///depends on the evol_mode argument.
			const std::valarray<double> &orbit, 

			///The rate of change of the entries in orbit per the relevant
			///system of differential equations.
			const std::valarray<double> &derivatives,

			///The rates of change of various quantities at the last step,
			///split into locked and not locked components and per body.
			///The convention to follow is that the planet should always be
			///body2. This argument bust be ignored by all stopping
			///conditions evaluated at evolution modes where the planet is
			///not present.
			const TidalDissipation &tidal_dissipation,

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
	std::valarray<double> operator()(
			const std::valarray<double> &,
			const std::valarray<double> &,
			const TidalDissipation &,
			const StellarSystem&,
			std::valarray<double> &stop_deriv,
			EvolModeType) const
	{stop_deriv.resize(1, NaN); return std::valarray<double>(1, 1);}

	///Identifies the condition as NO_STOP.
	StoppingConditionType type(unsigned =0) const {return NO_STOP;}
};

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
	short __body_index;

public:
	///Create the synchronization condition for the given planet.
	SynchronizedCondition(
			///The mutiplier in front of the orbital frequency in the lock.
			int orbital_freq_mult,

			///The multiplier in front of the spin frequency in the lock.
			int spin_freq_mult,

			///Which body's spin is checked for locking.
			short body_index) :
		__orbital_freq_mult(orbital_freq_mult),
		__spin_freq_mult(spin_freq_mult), __body_index(body_index) {}

	///\brief Returns the difference between the orbital and multiplier
	///scaled stellar spin angular velocities divided by the orbital angular
	///velocity.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	///
	///The evolution mode must be FAST_PLANET, SLOW_PLANET or LOCKED_TO_DISK
	///or NO_PLANET.
	std::valarray<double> operator()(
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const TidalDissipation &tidal_dissipation,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Identify this as a SYNCHRONIZED condition.
	StoppingConditionType type(unsigned =0) const {return SYNCHRONIZED;}

	///The multiplier in front of the orbital frequency in the lock.
	int orbital_frequency_multiplier() const {return __orbital_freq_mult;}

	///The multiplier in front of the spin frequency in the lock.
	int spin_frequency_multiplier() const {return __spin_freq_mult;}

	///Which body's spin is checked for locking.
	short body_index() const {return __body_index;}

};

///\brief Satisfied when the maximum tidal torque that the planet can exert
///on the star is no longer sufficient to keep the lock.
///
///\ingroup OrbitSolver_group
class BreakLockCondition : public StoppingCondition {
public: 
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
	std::valarray<double> operator()(
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const TidalDissipation &tidal_dissipation,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///The number of subconditions in the current condition.
	virtual size_t num_subconditions() const {return 2;}

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
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const TidalDissipation &tidal_dissipation,
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
	std::valarray<double> operator()(
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const TidalDissipation &tidal_dissipation,
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
	///The non synchronized conditions that are to be combined.
	std::vector<const StoppingCondition *> __generic_sub_conditions;

	///The synchronized conditions that are to be combined.
	std::vector<const SynchronizedCondition *> __sync_sub_conditions;

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

	///Updates all private methods except the sub_condition lists for adding
	///rhs.
	void update_meta_information(const StoppingCondition *rhs);

	///Adds the values of the given sub-condition to the given array.
	void add_subcondition_values(
			///The stopping condition to process.
			const StoppingCondition *cond,

			///See operator()
			const std::valarray<double> &orbit,
			
			///See operator()
			const std::valarray<double> &derivatives,

			///See operator()
			const TidalDissipation &tidal_dissipation,

			///See operator()
			const StellarSystem &system,

			///See operator()
			EvolModeType evol_mode,

			///The index in the values and derivs arrays where to start
			///inserting. On output it is updated to the position immediately
			///after the last insert.
			size_t &first_index,

			///An array to overwrite with the sub-condition values.
			std::valarray<double> &values,

			///An array to overwrite with the sub-condition derivatives.
			std::valarray<double> &derivs) const;
public:
	///Create an empty stopping condition (identical to NoStopCondition).
	CombinedStoppingCondition() :
		__generic_sub_conditions(), __sync_sub_conditions(),
		__num_subconditions(0), __delete_subcond(true) {}

	///Adds the conditions in RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(
			const CombinedStoppingCondition &rhs);

	///Adds  RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(const StoppingCondition *rhs);

	///Adds  RHS to the conditions of *this special case for synchronization.
	CombinedStoppingCondition &operator|=(
			const SynchronizedCondition *rhs);

	///\brief Returns the values of all stopping sub_conditions, not
	///necessarily in the order added.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	std::valarray<double> operator()(
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const TidalDissipation &tidal_dissipation,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Disables the destruction of the subconditions when *this is destroyed.
	void no_delete_subcond() {__delete_subcond=false;}

	virtual size_t num_subconditions() const {return __num_subconditions;}

	StoppingConditionType type(unsigned index=0) const
	{return __types[index];}

	///Returns the index-th condition if it is of SYNCHRONIZED type.
	const SynchronizedCondition &get_sync_condition(unsigned index)
	{
#ifdef DEBUG
		assert(index>=__generic_sub_conditions.size() &&
				index<__num_subconditions);
#endif
		return *__sync_sub_conditions[index-__generic_sub_conditions.size()];
	}

	double interpolation_range(unsigned index=0) const
	{return __interpolation_ranges[index];}

	///\brief Deletes all subconditions, unless no_delete_subcond has been
	///previously called.
	~CombinedStoppingCondition();
};

///Returns the spin frequency of the stellar convective zone in rad/day.
///
///The system age must already be set.
double convective_frequency(
		///The stellar system.
		const StellarSystem &system, 

		///The present orbit as appropriate for the evol_mode parameter.
		const std::valarray<double> &orbit,
		
		///The present evolution mode.
		EvolModeType evol_mode);

#endif
