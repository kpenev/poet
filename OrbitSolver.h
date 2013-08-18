#ifndef __ORBIT_SOLVER_H
#define __ORBIT_SOLVER_H

#include "AstronomicalConstants.h"
#include "Common.h"
#include "StellarSystem.h"
#include <math.h>
#include <list>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <sstream>
#include <iostream>
#include <limits>

///A class that solves the system of ODEs describing the orbital evolution of
///a single planet around a single star.

typedef int (*GSL_ODE_TYPE)(double, const double*, double*, void*);
typedef int (*GSL_JAC_TYPE)(double, const double*, double*, double*, void*);
typedef bool (*STOP_EVOL_TYPE)(double, const double*, void*);

enum EvolModeType {FAST_PLANET=-1, LOCKED_TO_PLANET, SLOW_PLANET, NO_PLANET,
	LOCKED_TO_DISK, TABULATION};

enum StoppingConditionType {NO_STOP, SYNCHRONIZED, BREAK_LOCK, PLANET_DEATH,
	WIND_SATURATION, ROT_FAST};

enum EvolVarType {AGE=-1, SEMIMAJOR, LCONV, LRAD};

///More civilized output for EvolModeType variables.
std::ostream &operator<<(std::ostream &os, const EvolModeType &evol_mode);

///More civilized output for StoppingConditionType variables.
std::ostream &operator<<(std::ostream &os,
		const StoppingConditionType &stop_cond_type);

///More civilized output for EvolVarType variables.
std::ostream &operator<<(std::ostream &os, const EvolVarType &evol_var);

///Define the minimum requirements of all stopping conditions.
class StoppingCondition {
private:
	double __interpolation_range;
public:
	StoppingCondition(double interp_range=Inf) :
		__interpolation_range(interp_range) {}

	virtual std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit, 
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const=0;

	virtual size_t num_subconditions() const {return 1;}

	virtual StoppingConditionType type(unsigned index=0) const=0;

	///A if a zero crossing is found only points within
	///__interpolation_range*(t_after-t_before) should be considered, where
	///t_after is tha age of the first point found after the zero crossing
	///and t_before is the age of the last point before the zero crossing.
	///In addition the interpolation should only include monotonic points, in
	///addition if derivatives are provided, no points other than the ones at
	///t_before and t_after should be considered.
	virtual double interpolation_range(unsigned index=0) const
	{return __interpolation_range;}
};

///A stopping condition that is never satisfied.
class NoStopCondition : public StoppingCondition {
public:
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const
	{stop_deriv.resize(1, NaN); return std::valarray<double>(1, 1);}
	StoppingConditionType type(unsigned index=0) const {return NO_STOP;}
};

///A stopping condition that is satisfied when the orbit and stellar
///rotation are synchronized. The orbit should be describe in terms of
///(a^6.5, Lconv, Lrad).
class SynchronizedCondition : public StoppingCondition{
private:
	///The planet whose orbit is checked for synchronization with the star.
	const Planet *planet;
	
	double __initial_semimajor;
public:
	///Create the synchronization condition for the given planet.
	SynchronizedCondition(const Planet *p, double initial_semimajor) :
		planet(p), __initial_semimajor(initial_semimajor) {}

	///Returns the difference between the orbital and stellar spin angular
	///velocities divided by the orbital angular velocity.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;
	StoppingConditionType type(unsigned index=0) const {return SYNCHRONIZED;}
};

///A stopping condition that is satisfied when the maximum tidal torque that
///the planet can exert on the star is no longer sufficient to keep the lock.
class BreakLockCondition : public StoppingCondition {
public: 
	///Returns the difference between the maximum tidal evolution of the
	///semimajor axis and the evolution required to keep the star spinning
	///synchronously with the orbit divided by the evolution required to keep
	///the lock.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;
	StoppingConditionType type(unsigned index=0) const {return BREAK_LOCK;}
};

///A stopping condition that is satisfied when the planet enters below either
///the roche sphere or the stellar photosphere.
class PlanetDeathCondition : public StoppingCondition {
public:
	///Returns the difference between the semimajor axis and the larger of
	///the roche radius and the stellar radius divided by the latter.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;
	StoppingConditionType type(unsigned index=0) const {return PLANET_DEATH;}
};

///A stopping condition that is satisfied when the stellar convective zone is
///spinning at exactly the wind saturation frequency.
class WindSaturationCondition : public StoppingCondition {
public:
	///Returns the difference between the convective angular velocity and the
	///wind saturation angular velocity divided by the wind saturation
	///angular velocity.
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;
	StoppingConditionType type(unsigned index=0) const
	{return WIND_SATURATION;}
};

//A stopping condition satisfied when the star is rotating faster than the
//threshold
class RotFastCondition : public StoppingCondition {
	///Returns the difference between the maximum tidal evolution of the
	///semimajor axis and the evolution required to keep the star spinning
	///synchronously with the orbit divided by the evolution required to keep
	///the lock.
private:
	double spin_thres;
public:
	RotFastCondition(double spin_thres): spin_thres(spin_thres) {}
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;
	StoppingConditionType type(unsigned index=0) const {return ROT_FAST;}
};

///A class combining the the outputs of multiple stopping conditions.
class CombinedStoppingCondition : public StoppingCondition {
private:
	///The conditions that are to be combined
	std::vector<const StoppingCondition *> sub_conditions;

	///Whether to delete the sub-conditions then *this is destroyed.
	bool delete_subcond;
public:
	///Create an empty stopping condition (identical to NoStopCondition).
	CombinedStoppingCondition() : sub_conditions(), delete_subcond(true) {}

	///Adds the conditions in RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(
			const CombinedStoppingCondition &rhs);

	///Adds  RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(const StoppingCondition *rhs);

	///Returns the values of all stopping sub_conditions
	std::valarray<double> operator()(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system,
			std::valarray<double> &stop_deriv,
			EvolModeType evol_mode) const;

	///Disables the destruction of the subconditions when *this is destroyed.
	void no_delete_subcond() {delete_subcond=false;}

	virtual size_t num_subconditions() const {return sub_conditions.size();}

	StoppingConditionType type(unsigned index=0) const
	{return sub_conditions[index]->type();}

	double interpolation_range(unsigned index=0) const
	{return sub_conditions[index]->interpolation_range();}

	///Deletes all subconditions, unless no_delete_subcond has been
	///previously called.
	~CombinedStoppingCondition();
};

///A wrapper tha allows the stellar system differential equation to be passed
///to the GSL ODE solver.
int stellar_system_diff_eq(double age, const double *orbital_parameters,
		double *orbital_derivatives, void *system_mode);

///A wrapper tha allows the stellar system jacobian to be passed
///to the GSL ODE solver.
int stellar_system_jacobian(double age, const double *orbital_parameters,
		double *param_derivs, double *age_derivs, void *system_mode);

///A structure to hold the information about why and where the evolution
///should stop.
class StopInformation {
private:
	double __stop_age, __stop_condition_precision;
	StoppingConditionType __stop_reason;
	bool __is_crossing;
	size_t __stop_condition_index;
public:
	StopInformation(double stop_age=Inf, double stop_precision=NaN,
			StoppingConditionType stop_reason=NO_STOP,
			bool is_crossing=false, size_t stop_condition_index=0) :
		__stop_age(stop_age),
		__stop_condition_precision(stop_precision),
		__stop_reason(stop_reason),
		__is_crossing(is_crossing),
		__stop_condition_index(stop_condition_index) {}

	double stop_age() const {return __stop_age;}
	double &stop_age() {return __stop_age;}

	double stop_condition_precision() const
	{return __stop_condition_precision;}

	double &stop_condition_precision()
	{return __stop_condition_precision;}

	StoppingConditionType stop_reason() const {return __stop_reason;}
	StoppingConditionType &stop_reason() {return __stop_reason;}

	bool is_crossing() const {return __is_crossing;}
	bool &is_crossing() {return __is_crossing;}

	size_t stop_condition_index() const {return __stop_condition_index;}
	size_t &stop_condition_index() {return __stop_condition_index;}
};

std::ostream &operator<<(std::ostream &os, const StopInformation &stop);

class ExtremumInformation {
private:
	double __x, __y;
public:
	ExtremumInformation(double x=Inf, double y=NaN) : __x(x), __y(y) {}

	double x() const {return __x;}
	double &x() {return __x;}

	double y() const {return __y;}
	double &y() {return __y;}
};

class StopHistoryInterval {
private:
	size_t __num_points,///< Number of points in the interval
		   __point_i;///< The index of the current point

	///The first age in the interval
	std::list<double>::const_iterator __first_age,

		///The last age in the interval.
		__last_age,

		///The one past last element of the history of stoppnig condition
		///ages.
		__history_age_end,

		///The first age in the discarded stopping conditions.
		__discarded_age_begin,
		
		///The age of the current point.
		__age_i;

	///The first stopping condition value in the interval
	std::list< std::valarray<double> >::const_iterator __first_stop_cond,

		///The last stopping condition value in the interval
		__last_stop_cond,

		///The one past last element of the history of stoppnig conditions
		__stop_cond_history_end,

		///The first of the discarded stopping condition values.
		__stop_cond_discarded_begin,
		
		///The current stopping condition value.
		__stop_cond_i,

		///The first stopping condition derivative in the interval
		__first_stop_deriv,

		///The last stopping condition derivative in the interval
		__last_stop_deriv,
		
		///The one past last element of the history of stoppnig conditions
		__stop_deriv_history_end,

		///The first of the discarded stopping condition values.
		__stop_deriv_discarded_begin,
			
		///The current stop condition derivative
		__stop_deriv_i;

	///Increments all the iterators passed arguments, taking care of the
	///switch from history to discarded if necessary.
	void advance_iterator_set(std::list<double>::const_iterator &age_i,
			std::list< std::valarray<double> >::const_iterator &cond_i,
			std::list< std::valarray<double> >::const_iterator &deriv_i);

	///Decrements all the iterators passed arguments, taking care of the
	///switch from discarded to history if necessary.
	void retreat_iterator_set(std::list<double>::const_iterator &age_i,
			std::list< std::valarray<double> >::const_iterator &cond_i,
			std::list< std::valarray<double> >::const_iterator &deriv_i);

public:
	StopHistoryInterval(
			///The number of points in the interval.
			size_t num_points=0,

			///An iterator pointing to the first age in the interval.
			std::list<double>::const_iterator
			first_age=std::list<double>::const_iterator(),

			///An iterator pointing to one past the last age in the history.
			std::list<double>::const_iterator
			history_age_end=std::list<double>::const_iterator(),

			///An iterator pointing to the first age in the discarded list.
			std::list<double>::const_iterator
			discarded_age_begin=std::list<double>::const_iterator(),	

			///An iterator pointing to the first stopping condition in the
			///interval.
			std::list< std::valarray<double> >::const_iterator
			first_stop_cond=std::list< std::valarray<double> >::const_iterator(),

			///An iterator pointing to one past the last stopping condition
			///in the history.
			std::list< std::valarray<double> >::const_iterator
			stop_cond_history_end=
			std::list< std::valarray<double> >::const_iterator(),

			///An iterator pointing to the first stopping condition in the
			///discarded list.
			std::list< std::valarray<double> >::const_iterator
			stop_cond_discarded_begin=
			std::list< std::valarray<double> >::const_iterator(),

			///An iterator pointing to the first stopping derivative in the
			///interval.
			std::list< std::valarray<double> >::const_iterator
			first_stop_deriv=
			std::list< std::valarray<double> >::const_iterator(),

			///An iterator pointing to one past the last stopping derivative
			///in the history.
			std::list< std::valarray<double> >::const_iterator
			stop_deriv_history_end=
			std::list< std::valarray<double> >::const_iterator(),

			///An iterator pointing to the first stopping derivative in the
			///discarded list.
			std::list< std::valarray<double> >::const_iterator
			stop_deriv_discarded_begin=
			std::list< std::valarray<double> >::const_iterator());

	StopHistoryInterval(const StopHistoryInterval &orig);

	///Makes the current point the first point in the interval
	void reset();

	///Advances to the next point in the interval.
	StopHistoryInterval &operator++();

	///Advances to the next point in the interval.
	StopHistoryInterval operator++(int);

	///Advances to the next point in the interval.
	StopHistoryInterval &operator--();

	///Advances to the next point in the interval.
	StopHistoryInterval operator--(int);

	///Moves the entire interval, along with the current point left n points,
	///gaining n new points at the front and losing n at the back. If there
	///are not enough points in the history undefined behavior results.
	StopHistoryInterval &operator<<(size_t n);

	///Moves the entire interval, along with the current point right n
	///points, gaining n new points at the back and losing n at the front. If
	///there are not enough points in the discarded list undefined behavior
	///results.
	StopHistoryInterval &operator>>(size_t n);

	///Copies rhs to this.
	StopHistoryInterval &operator=(const StopHistoryInterval &rhs);

	///Checks if the RHS is the same interval and is at the same point in it.
	bool operator==(const StopHistoryInterval &rhs);

	///Adds the n points before the first point to the interval.
	void grow_left(size_t n=1);

	///Adds the n points before the first point to the interval.
	void grow_right(size_t n=1);

	///Returns the number of points in the interval.
	size_t num_points() {return __num_points;}

	///Returns the index of the current point within the interval.
	size_t current_point_index() {return __point_i;}

	///Returns the number of conditions at the first point.
	size_t number_conditions() {return __first_stop_cond->size();}

	///Returns true iff this is the invalid point marking the end of the
	///interval.
	bool end() {return __point_i==__num_points;}

	///Returns the age of the first point in the interval.
	double first_age() const {return *__first_age;}

	///Returns the age of the last point in the interval.
	double last_age() const {return *__last_age;}

	///Returns the age of the current point.
	double age() const {return *__age_i;}

	///Returns the value of the stop condition with the given index for the
	///first point in the interval.
	double first_stop_condition_value(size_t condition_index) const
	{return (*__first_stop_cond)[condition_index];}

	///Returns the value of the stop condition with the given index for the
	///last point in the interval.
	double last_stop_condition_value(size_t condition_index) const
	{return (*__last_stop_cond)[condition_index];}

	///Returns the value of the stop condition with the given index for the
	///current point.
	double stop_condition_value(size_t condition_index) const
	{return (*__stop_cond_i)[condition_index];}

	///Returns the derivative of the stop condition with the given index for
	///the first point in the interval.
	double first_stop_condition_deriv(size_t condition_index) const
	{return (*__first_stop_deriv)[condition_index];}

	///Returns the derivative of the stop condition with the given index for
	///the last point in the interval.
	double last_stop_condition_deriv(size_t condition_index) const
	{return (*__last_stop_deriv)[condition_index];}

	///Returns the derivative of the stop condition with the given index for
	///the current point.
	double stop_condition_deriv(size_t condition_index) const
	{return (*__stop_deriv_i)[condition_index];}
};

std::ostream &operator<<(std::ostream &os, StopHistoryInterval interval);

class OrbitSolver {
private:
	static const double MAX_END_AGE = 10;
	///The first and last age for which evolution is required
	double start_age, end_age, 
	       precision,///< The precision required of the solution

		   ///A threshold for the stellar spin to note while calculating the evolution
	       spin_thres, 

		   ///A definition of the ZAMS age (a point at this age is forced to
		   ///be included in the tabulated evolution.
		   main_seq_start;

	///The ages at which solution is tabulated
	std::list<double> tabulated_ages;

	///The evolution mode corresponding to the matching tabulated age.
	std::list<EvolModeType> tabulated_evolution_mode;

	///The orbital ODE variables at the tabulated ages
	std::vector< std::list<double> > tabulated_orbit, 

	///The derivatives of the orbital ODE variables
	tabulated_deriv;

	///How many points after the start of the history have to be skipped when
	///lookng for a zero crossing for each condition
	std::valarray<size_t> skip_history_zerocrossing;

	///The age after which to look for extrema for each condition
	std::valarray<double> skip_history_extremum;


	///The ages at which the stop condition history is kept
	std::list<double> stop_history_ages,
		discarded_stop_ages;

	std::list< std::valarray<double> >
		orbit_history,///< Past orbits
		orbit_deriv_history, ///< Past orbital derivatives
		stop_cond_history,///< Past values of the stop conditions
		stop_deriv_history,///< Past values of the stop condition derivatives

		///Discarded values of the stop conditions (useful for
		///interpolating to zeroes and extrema).
		stop_cond_discarded,

		///Discarded derivatives of the stop conditions (useful for
		///interpolating to zeroes and extrema).
		stop_deriv_discarded;

	///Generates a nicely formatted table of the contents of the discarded
	///and history stopping condition information.
	void output_history_and_discarded(std::ostream &os);

	///Removes all stored discarded stop condition information.
	void clear_discarded();

	///Inserts a new entry in the discarded ages, stop conditions and
	///derivatives, making sure the ages remain ordered.
	void insert_discarded(double age,
			const std::valarray<double> &current_stop_cond,
			const std::valarray<double> &current_stop_deriv);

	///Appends the given orbit and derivatives to tabulated_orbit and
	///tabulated_deriv respectively assuming the orbit contains the
	///quantities evolved for the given evolution mode.
	void append_to_orbit(const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			EvolModeType evolution_mode, double age,
			const StellarSystem &system,
			double planet_formation_semimajor=NaN);

	///Clears the current stopping condition history.
	void clear_history();

	///Rewinds the evlution to the last step before the given age, setting 
	///the orbit and derivatives to what they were at that step and removing
	///any items from the histories and tabulations that are later than 
	///max_age. Returns the age of the last non-erased step.
	double go_back(double max_age, std::valarray<double> &orbit,
			std::valarray<double> &derivatives);

	///Returns the dimension of the ODEs governing the evolution of the
	///given type.
	size_t ode_dimension(EvolModeType evolution_mode);

	///Finds the smallest possible interval that contains a zero crossing/or
	///an extremum straddling the history and discarded stop conditions, 
	///containing at most the specified number of points (could be less if 
	///there are not enough points). The interval is also guaranteed to
	///contain at least one point in the history and one point in the
	///discarded list.
	StopHistoryInterval select_stop_condition_interval(bool crossing, 
			size_t cond_ind, size_t max_points) const;

	///Estimates the value and age of an extremum if it potentially can cross
	///a zero by using the last and past tabulated points. If no extremum is
	///indicated by the points, or if it is in the wrong direction, returns
	///the result of the default constructor of ExtremumInformation.
	ExtremumInformation extremum_from_history_no_deriv(
			size_t condition_index) const;

	///Estimates the value and age of an extremum if it potentially can cross
	///a zero by using the last and past tabulated points. If no extremum is
	///indicated by the points, or if it is in the wrong direction, returns
	///the result of the default constructor of ExtremumInformation.
	ExtremumInformation extremum_from_history(size_t condition_index) const;

	///Estimates the age at which the stopping condition with the given index
	///crossed zero. If no zero-crossing is indicated, Inf is returned.
	double crossing_from_history_no_deriv(size_t condition_index) const;

	///Estimates the age at which the stopping condition with the given index
	///crossed zero. If no zero-crossing is indicated, Inf is returned.
	double crossing_from_history(size_t condition_index) const;

	///Initializes the skip_history_zerocrossing and skip_history_extremum
	///arrays appropriately after a mode change.
	void initialize_skip_history(const StoppingCondition &stop_cond,
			StoppingConditionType stopp_reason);

	///Updates the skip_history_zerocrossing and skip_history_extremum
	///arrays appropriately after an acceptable step.
	void update_skip_history(
			const std::valarray<double> &current_stop_cond,
			const StoppingCondition &stop_cond,
			const StopInformation &stop_info);

	///Return true iff the step with the corresponding stop information is
	///acceptable.
	bool acceptable_step(double age, const StopInformation &stop_info);

	///This function is called after every GSL step. It updates the
	///stop_cond_history and stop_deriv_history variables appropriately.
	///Returns the full information about the closest estimated age where a
	///condition is zero or an extremum exists which might have crossed zero.
	StopInformation update_stop_condition_history(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system, EvolModeType evolution_mode,
			const StoppingCondition &stop_cond,
			StoppingConditionType stop_reason=NO_STOP);

	///Calculates the evolution of a planetary system starting at an initial
	///age of start_age, with orbital parameters start_orbit, where the star
	///and the planet should already be members of the given stellar system,
	///stopping the evolution if an age of max_age is reached or a quantity
	///(stop cond) that evolves with age crosses zero.
	///
	///The orbital parameters should be:
	/// * semimajor axis^6.5, convective and radiative angular momenta if
	///   evolution_mode is FAST_PLANET or SLOW_PLANET
	/// * semimajor axis and radiative angular momentum if evolution_mode
	///   is LOCKED_TO_PLANET
	/// * convective and radiative angular momenta if evolution_mode is
	///   NO_PLANET
	/// * radiative angular momentum if evolution mode is LOCKED_TO_DISK
	///
	///The stopping condition quantity should be a smooth function of age
	///and should support evalution with
	///arguments:age, orbital parameters, orbital_derivatives and
	///stellar system as well as a last return argument that should be set
	///to the total derivative of the quantity with respect to age if it is
	///available or should be left unchanged if it is not.
	///
	///Does not reset the current orbit, but appends to the currently
	///tabulated values.
	///
	///On input step_size should be the initial step size required an on exit
	///it is the size of the next suggested step.
	///
	///An exit max_age is overwritten with the age at which the evolution
	///was stopped, and orbit is overwritten with the orbital state at which
	///the evolution stopped.
	///
	///The stop_reason argument on input should be the reason why the 
	///last evolution stopped. It should be NO_STOP if this is the first
	///piece of evolution being calculated. On exit it is overwritten
	///with the value appropriate for the next run.
	///The return value is true if the last step finished after the stopping
	///condition crossed zero and false if it ended before that.
	bool evolve_until(StellarSystem *system, double start_age,
			double &max_age, std::valarray<double> &orbit,
			double &stop_condition_value, StoppingConditionType &stop_reason,
			double max_step=Inf, EvolModeType evolution_mode=FAST_PLANET,
			WindSaturationState wind_state=UNKNOWN,
			const StoppingCondition &stop_cond=NoStopCondition(),
			double planet_formation_semimajor=NaN);

	///Returns the stopping condition which ends the given evolution mode.
	CombinedStoppingCondition 
		*get_stopping_condition(EvolModeType evolution_mode,
				double initial_semimajor,
				const Planet *planet) const;

	///Returns the evolution mode that the system is entering, assuming that
	///some critical age is reached (e.g. the disk dissipated). The last 
	///orbital state should be orbit (in the old evolution mode), the 
	///semimajor at which the planet starts after the disk dissipates is
	///initial_semimajor (in AU) and the previous mode was evolution_mode.
	EvolModeType critical_age_evol_mode(double age, 
			const std::valarray<double> &orbit,
			double initial_semimajor, const StellarSystem &system, 
			EvolModeType evolution_mode, double planet_formation_age) const;

	///Returns the evolution mode that the system is entering, assuming that
	///the last orbital state is orbit (in the old evolution mode), the
	///semimajor at which the planet starts after the disk dissipates is
	///initial_semimajor (in AU) and the previous mode was evolution_mode.
	EvolModeType next_evol_mode(double age, 
			const std::valarray<double> &orbit,
			double initial_semimajor, const StellarSystem &system, 
			EvolModeType evolution_mode,
			StoppingConditionType condition_type,
			double condition_value, bool stopped_before,
			double planet_formation_age) const;

	///Returns what age the evolution with the given mode should stop if
	///no other stopping condition occurs.
	double stopping_age(double age, EvolModeType evolution_mode,
			const StellarSystem &system, double planet_formation_age);

	///Transforms orbital parameters from one evolution mode (from_mode) to
	///another (to_mode).
	std::valarray<double> transform_orbit(EvolModeType from_mode,
			EvolModeType to_mode, double age,
			const std::valarray<double> &from_orbit,
			double initial_semimajor, const StellarSystem &system) const;

	///Transforms the deriatives of the orbital parameters from one evolution
	///mode (from_mode) to another (to_mode).
	std::valarray<double> transform_derivatives(EvolModeType from_mode,
			EvolModeType to_mode, double age,
			const std::valarray<double> &from_orbit, 
			const std::valarray<double> &from_deriv,
			double initial_semimajor, const StellarSystem &system);

	///Clears any previously calculated evolution.
	void reset();
public:
	///Prepare to solve for the orbital evolution over the given range
	///and to the required precision.
	OrbitSolver(double min_age, double max_age, double required_precision,
			double spin_thres=4*M_PI, double main_seq_start=0.1);

	///Actually solves the given differential equation with the given
	///boundary conditions.
	void operator()(
			///The stellar system to calculate the evolution for
			StellarSystem &system,

			///The maximum size of the time steps allowed (useful if finer
			///sampling of the output than default is necessary).
			double max_step=Inf,

			///The age at which a planet magically appears in a perfectly
			///circularized orbit and starts affecting the system. By
			///default, the planet is assumed to always be there and
			///planet_formation_semimajor has no effect.
			double planet_formation_age=0,

			///The semimajor axis at which the planet first appears. This
			///argument is ignored if planet_formation_age<=start_age.
			double planet_formation_semimajor=NaN,

			///The age at which to start the evolution. Use NaN (default) to
			///start when the radiative core first starts to appear, in
			///which case, start_orbit and initial_evol_mode should be left
			///at their default values as well, but planet_formation_age and
			///planet_formation_semimajor must be specified.
			double start_age=NaN, 

			///The initial evolution mode of the system
			EvolModeType initial_evol_mode=LOCKED_TO_DISK,

			///The initial state to start the system in. The contents
			///depends on initial_evol_mode.
			const std::valarray<double>
				&start_orbit=std::valarray<double>(0.0, 1),
			
			///If given, the evolution mode is kept constant regardless
			///of what actually occurs with the evolution
			bool no_evol_mode_change=false);


	/* Evolves the given system from start to end, with its initial orbit
	 * given in y.  The final orbit is put into y.*/
	void evolve_to(double start, double end, double* y,
			StellarSystem* system) {};

	///Helper function that returns the time, in Gyr, that the convective
	///rotation rate is above spin_thres.  The other arguments are identical
	///to those of operator().
	double fast_time(StellarSystem &system,
			double max_step=Inf, double planet_formation_age=0,
			double planet_formation_semimajor=NaN, double start_age=NaN,
			EvolModeType initial_evol_mode=LOCKED_TO_DISK,
			std::valarray<double>
				start_orbit=std::valarray<double>(0.0, 1));

	const std::list<double>
		*get_tabulated_var(EvolVarType var_type) const;

	///Returns the derivative of the independent variable at the points
	///specified by get_tabulated_indep_var.
	const std::list<double> *get_tabulated_var_deriv(
			EvolVarType var_type) const;

	///Returns a list of the evolution modes.
	const std::list<EvolModeType> *get_tabulated_evolution_mode() const
	{return &tabulated_evolution_mode;}

	/*returns the L2 norm of the last fractional error between the last
	simulated orbit, and the arguments (a, Lc). If the last simulated age
	is not age, returns NaN. */
	double last_error(double age, double a, double Lc);
};

///Converts a semimajor axis given in AU to the variable used in solving for
///the orbital evolution.
inline double AU_to_Rsun6p5(double au) {
	return std::pow(au*AstroConst::AU/AstroConst::solar_radius, 6.5);
};

///Converts the variable used in solving for the orbital evolution to a
///semimajor axis given in AU.
inline double Rsun6p5_to_AU(double rsun) {
  return std::pow(rsun, 1.0/6.5)*AstroConst::solar_radius/AstroConst::AU;
};

#endif
