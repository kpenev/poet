/**\file
 *
 * \brief Defines the OrbitSolver class, the various stopping conditions and
 * a number of other classes used while calculating the orbital evolution.
 * 
 * \ingroup OrbitSolver_group
 */

#ifndef __ORBIT_SOLVER_H
#define __ORBIT_SOLVER_H

#include "AstronomicalConstants.h"
#include "Common.h"
#include "StellarSystem.h"
#include "StoppingConditions.h"
#include "ExternalStoppingConditions.h"
#include "TidalDissipation.h"
#include "OrbitalExpressions.h"
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

typedef int (*GSL_ODE_TYPE)(double, const double*, double*, void*);
typedef int (*GSL_JAC_TYPE)(double, const double*, double*, double*, void*);
typedef bool (*STOP_EVOL_TYPE)(double, const double*, void*);


///Tags for the variables being evolved.
enum EvolVarType {
	AGE=-1,		///< The age.
	SEMIMAJOR, 	///< The semimajor axis.
	INCLINATION,///< Surface spin to orbital angular momentum angle
	LCONV, 		///< Angular momentum of the convective envelope.
	LRAD_PAR,	///< Angular momentum of the core along LCONV
	LRAD_PERP,	///< Angular momentum of the core perpendicular to LCONV
	NUM_EVOL_VAR///< The number of evolution variables.
};

///More civilized output for EvolVarType variables.
std::ostream &operator<<(std::ostream &os, const EvolVarType &evol_var);

///\brief (Re-)Initializes the given tidal dissipation.
///
///Must only be called when the planet is present and after the stellar
///system has been assigned an age.
void get_tidal_dissipation(
		///The planet-star system evolving.
		StellarSystem &system, 

		///The age when tidal dissipation must be computed.
		double age,

		///The current values of the parameters being evolved.
		const double *parameters,

		///The tidal dissipation variable to (re-)initialize.
		TidalDissipation &dissipation);

///\brief A wrapper tha allows the stellar system differential equation to be
///passed to the GSL ODE solver.
int stellar_system_diff_eq(
		///System age in Gyr.
		double age,
		
		///A C-style array of the variables evolved under the current set of
		///equations.
		const double *parameters,

		///A C-style array of the rates of change of the orbital_parameters.
		double *derivatives,
		
		///A C-style array of pointers giving the stellar system, the
		///current evolution mode, the wind saturation state to assume,
		///the spin-orbit lock of the star and the dissipation rates.
		void *system_mode_windstate_lock_dissipation);

///\brief A wrapper tha allows the stellar system jacobian to be passed
///to the GSL ODE solver.
int stellar_system_jacobian(
		///System age in Gyr.
		double age,
		
		///A C-style array of the variables evolved under the current set of
		///equations.
		const double *parameters,

		///A C-style array of the Jacobian matrix entries.
		double *param_derivs,
		
		///A C-style array of the partial age derivatives of the differential
		///equations. Age is in Gyr.
		double *age_derivs,
		
		///A C-style array of pointers giving the stellar system and the
		///current evolution mode.
		void *system_mode);

///\brief The information about why and where the evolution should stop.
///
///\ingroup OrbitSolver_group
class StopInformation {
private:
	///The target stopping age in Gyr.
	double __stop_age,

		   ///\brief The precision up to which the reason to stop is satisfied.
		   ///
		   ///If we are stopping because some stopping condition is actually
		   ///satisfied, this is the value of the condition. If we are
		   ///stopping because of an extremum this is the difference between
		   ///the estimated value of the condition at the extremum and the
		   ///value at the closest step.
		   __stop_condition_precision;

	///The reason for stopping.
	StoppingConditionType __stop_reason;

	///Is the reason for stopping that a condition actually crossed zero?
	bool __is_crossing, __crossed_zero;

	///The index of the condition which caused us to stop.
	size_t __stop_condition_index;
public:
	///\brief Create an object with the information about why we evolution
	///should be stopped.
	StopInformation(
			///The target stopping age in Gyr.
			double stop_age=Inf,
			
			///The precision up to which the reason to stop is satisfied.
			double stop_precision=NaN,

			///The reason for stopping.
			StoppingConditionType stop_reason=NO_STOP,

			///Is the reason for stopping that a condition actually crossed
			///zero?
			bool is_crossing=false,

			///Did we stop after a zero-crossing.
			bool crossed_zero=false,
			
			///The index of the condition which caused us to stop.
			size_t stop_condition_index=0) :
		__stop_age(stop_age),
		__stop_condition_precision(stop_precision),
		__stop_reason(stop_reason),
		__is_crossing(is_crossing),
		__crossed_zero(crossed_zero),
		__stop_condition_index(stop_condition_index) {}

	///Copy orig to *this.
	StopInformation(const StopInformation &orig) :
		__stop_age(orig.__stop_age),
		__stop_condition_precision(orig.__stop_condition_precision),
		__stop_reason(orig.__stop_reason),
		__is_crossing(orig.__is_crossing),
		__stop_condition_index(orig.__stop_condition_index) {}

	///The target stopping age in Gyr.
	double stop_age() const {return __stop_age;}

	///The target stopping age in Gyr.
	double &stop_age() {return __stop_age;}

	///\brief The precision up to which the reason to stop is satisfied.
	double stop_condition_precision() const
	{return __stop_condition_precision;}

	///\brief The precision up to which the reason to stop is satisfied.
	double &stop_condition_precision()
	{return __stop_condition_precision;}

	///The reason for stopping.
	StoppingConditionType stop_reason() const {return __stop_reason;}
	///The reason for stopping.
	StoppingConditionType &stop_reason() {return __stop_reason;}

	///Is the reason for stopping that a condition actually crossed zero?
	bool is_crossing() const {return __is_crossing;}
	///Is the reason for stopping that a condition actually crossed zero?
	bool &is_crossing() {return __is_crossing;}

	///The index of the condition which caused us to stop.
	size_t stop_condition_index() const {return __stop_condition_index;}
	///The index of the condition which caused us to stop.
	size_t &stop_condition_index() {return __stop_condition_index;}

	///\brief Did we stop after the stopping condition crossed zero (always
	///false for extrema).
	bool crossed_zero() const {return __crossed_zero;}
	///\brief Did we stop after the stopping condition crossed zero (always
	///false for extrema).
	bool &crossed_zero() {return __crossed_zero;}

	///Copy rhs to *this.
	StopInformation &operator=(const StopInformation &rhs)
	{
		__stop_age=rhs.__stop_age;
		__stop_condition_precision=rhs.__stop_condition_precision;
		__stop_reason=rhs.__stop_reason;
		__is_crossing=rhs.__is_crossing;
		__stop_condition_index=rhs.__stop_condition_index;
		return *this;
	}
};

///Civilized output of a StopInformation object.
std::ostream &operator<<(std::ostream &os, const StopInformation &stop);

///\brief Infomation about an extremum of a function.
///
///\ingroup OrbitSolver_group
class ExtremumInformation {
private:
	///The value of the argument where the extremum occurs.
	double __x,
		   
		   ///The value of the function at the extremum.
		   __y;
public:
	///Create an object to hold information about a function extremum.
	ExtremumInformation(
			///The value of the argument where the extremum occurs.
			double x=Inf,

			///The value of the function at the extremum.
			double y=NaN) : __x(x), __y(y) {}

	///The value of the argument where the extremum occurs.
	double x() const {return __x;}
	///The value of the argument where the extremum occurs.
	double &x() {return __x;}

	///The value of the function at the extremum.
	double y() const {return __y;}
	///The value of the function at the extremum.
	double &y() {return __y;}
};

///\brief A collection of accepted and discarded evolution steps which
///contain some reason to stop.
///
///\ingroup OrbitSolver_group
class StopHistoryInterval {
private:
	size_t __num_points,///< Number of points in the interval
		   __point_i;///< The index of the current point

	///The first age in the interval.
	std::list<double>::const_iterator __first_age,

		///The last age in the interval.
		__last_age,

		///\brief The one past last element of the history of stoppnig
		///condition ages.
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

	///\brief Increments all the iterators passed as arguments.
	///
	///Takes care of the switch from history to discarded if necessary.
	void advance_iterator_set(std::list<double>::const_iterator &age_i,
			std::list< std::valarray<double> >::const_iterator &cond_i,
			std::list< std::valarray<double> >::const_iterator &deriv_i);

	///\brief Decrements all the iterators passed as arguments.
	///
	///Takes care of the switch from discarded to history if necessary.
	void retreat_iterator_set(std::list<double>::const_iterator &age_i,
			std::list< std::valarray<double> >::const_iterator &cond_i,
			std::list< std::valarray<double> >::const_iterator &deriv_i);

public:
	///Construct an interval of steps with some reason to stop the evolution.
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

	///Copy orig to *this.
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

	///\brief Moves the entire interval, along with the current point left n
	///points.
	///
	///Gains n new points at the front and loses n at the back.
	///
	///If there are not enough points in the history undefined behavior
	///results.
	StopHistoryInterval &operator<<(size_t n);

	///\brief Moves the entire interval, along with the current point right n
	///points.
	///
	///Gains n new points at the back and loses n at the front.
	///
	///If there are not enough points in the discarded list undefined
	///behavior results.
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

	///\brief Returns true iff this is the invalid point marking the end of
	///the interval.
	bool end() {return __point_i==__num_points;}

	///Returns the age of the first point in the interval.
	double first_age() const {return *__first_age;}

	///Returns the age of the last point in the interval.
	double last_age() const {return *__last_age;}

	///Returns the age of the current point.
	double age() const {return *__age_i;}

	///\brief Returns the value of the stop condition with the given index
	///for the first point in the interval.
	double first_stop_condition_value(size_t condition_index) const
	{return (*__first_stop_cond)[condition_index];}

	///\brief Returns the value of the stop condition with the given index
	///for the last point in the interval.
	double last_stop_condition_value(size_t condition_index) const
	{return (*__last_stop_cond)[condition_index];}

	///\brief Returns the value of the stop condition with the given index
	///for the current point.
	double stop_condition_value(size_t condition_index) const
	{return (*__stop_cond_i)[condition_index];}

	///\brief Returns the derivative of the stop condition with the given
	///index for the first point in the interval.
	double first_stop_condition_deriv(size_t condition_index) const
	{return (*__first_stop_deriv)[condition_index];}

	///\brief Returns the derivative of the stop condition with the given
	///index for the last point in the interval.
	double last_stop_condition_deriv(size_t condition_index) const
	{return (*__last_stop_deriv)[condition_index];}

	///\brief Returns the derivative of the stop condition with the given
	///index for the current point.
	double stop_condition_deriv(size_t condition_index) const
	{return (*__stop_deriv_i)[condition_index];}
};

///Civilized output of a StopHistoryInterval.
std::ostream &operator<<(std::ostream &os, StopHistoryInterval interval);

///\brief Solves the system of ODEs describing the evolution of a 
///single planet around a single star.
///
///\ingroup OrbitSolver_group
class OrbitSolver {
private:
	///A hard limit for the maximum age imposed on all solvers.
	static const double MAX_END_AGE;

	
	double end_age, ///< The last  age for which evolution is required.
	       precision,///< The precision required of the solution

		   ///\brief A threshold for the stellar spin to note while
		   ///calculating the evolution
	       spin_thres;

	bool adjust_end_age;

	///The ages at which solution is tabulated
	std::list<double> __tabulated_ages;

	///The evolution mode corresponding to the matching tabulated age.
	std::list<EvolModeType> __tabulated_evolution_mode;

    ///The wind saturation state corresponding to the matching tabulated age.
    std::list<WindSaturationState> __tabulated_wind_saturation;

	///The orbital ODE variables at the tabulated ages
	std::vector< std::list<double> > __tabulated_orbit, 

		///The derivatives of the orbital ODE variables
		__tabulated_deriv;

	std::list< SpinOrbitLockInfo > __tabulated_lock;

	///\brief The number of points at the start of the history to skip when
	///lookng for a zero crossing for each condition
	std::valarray<size_t> skip_history_zerocrossing;

	///The age after which to look for extrema for each condition
	std::valarray<double> skip_history_extremum;

	///The ages at which the stop condition history is kept
	std::list<double> stop_history_ages,

		///\brief The ages of steps which were discarded becauset they are
		///past a zero or an extremum of a stopping condition .
		discarded_stop_ages;

	std::list< std::valarray<double> >
		orbit_history,///< Past orbits
		orbit_deriv_history, ///< Past orbital derivatives
		stop_cond_history,///< Past values of the stop conditions
		stop_deriv_history,///< Past values of the stop condition derivatives

		///\brief Discarded values of the stop conditions. 
		///
		///Useful for interpolating to zeroes and extrema.
		stop_cond_discarded,

		///\brief Discarded derivatives of the stop conditions.
		///
		///Useful for interpolating to zeroes and extrema.
		stop_deriv_discarded;

	///The tidal dissipation with suitable locks for the current evolution.
	TidalDissipation __dissipation;

#ifdef DEBUG
	///\brief Generates a nicely formatted table of the contents of the
	///discarded and history stopping condition information.
	///
	///Only exists if DEBUG is \#defined.
	void output_history_and_discarded(std::ostream &os);
#endif

	///Removes all stored discarded stop condition information.
	void clear_discarded();

	///\brief Adds an entry in the discarded ages, stop conditions and
	///derivatives.
	///
	///Makes sure the ages remain ordered.
	void insert_discarded(double age,
			const std::valarray<double> &current_stop_cond,
			const std::valarray<double> &current_stop_deriv);

	///\brief Appends the given orbit and derivatives to tabulated_orbit and
	///tabulated_deriv respectively.
	void append_to_orbit(
			///The orbital state to append. The expected content depends on
			///the evolution_mode and star_lock arguments. See
			//#StellarSystem::differential_equations() for details.
			const std::valarray<double> &orbit,

			///The rates of change of the orbit entries per Gyr.
			const std::valarray<double> &derivatives,

			///The evolution mode represented in orbit.
			EvolModeType evolution_mode,

			///Whether the star is in some spin-orbit locked state, and which
			///one.
			const SpinOrbitLockInfo &star_lock,

			///The saturation state of the stellar wind.
			WindSaturationState wind_state,

			///The age of the evolving stellar system.
            double age,
			
			///The stellar system being evolved.
			const StellarSystem &system,

			///The semimajor axis at which the planet froms in AU.
			double planet_formation_semimajor=NaN,
			
			///The inclination with which the planet froms.
			double planet_formation_inclination=NaN);

	///Clears the current stopping condition history.
	void clear_history();

	///\brief Rewinds the evlution to the last step before the given age and
	///returns the age of that step.
	///
	///Sets the orbit and derivatives to what they were at that step and
	///removes any items from the histories and tabulations that are later
	///than max_age.
	double go_back(double max_age, std::valarray<double> &orbit,
			std::valarray<double> &derivatives);

	///\brief Returns the dimension of the ODEs governing the evolution of
	///the given type.
	size_t ode_dimension(EvolModeType evolution_mode,
			const SpinOrbitLockInfo &star_lock);

	///\brief Selects a history interval for interpolating to a reason to
	///stop the evolution.
	///
	///Finds the smallest possible interval that contains a zero crossing/or
	///an extremum straddling the history and discarded stop conditions, 
	///containing at most the specified number of points (could be less if 
	///there are not enough points). The interval is also guaranteed to
	///contain at least one point in the history and one point in the
	///discarded list.
	StopHistoryInterval select_stop_condition_interval(bool crossing, 
			size_t cond_ind, size_t max_points) const;

	///\brief Estimates the value and age of an extremum of a stopping
	///condition for which no derivative information is available.
	///
	///If no extremum is indicated by the points, or if it is in the wrong
	///direction to cause a zero-crossing, returns the result of the default
	///constructor of ExtremumInformation.
	ExtremumInformation extremum_from_history_no_deriv(
			size_t condition_index) const;

	///\brief Estimates the value and age of an extremum of a stopping
	///condition for which derivative information is available.
	///
	///If no extremum is indicated by the points, or if it is in the wrong
	///direction to cause a zero-crossing, returns the result of the default
	///constructor of ExtremumInformation.
	ExtremumInformation extremum_from_history(size_t condition_index) const;

	///\brief Estimates the age at which a stopping condition with no
	///derivative information crossed zero.
	///
	///If no zero-crossing is indicated, Inf is returned.
	double crossing_from_history_no_deriv(
			///The index of the condition to search.
			size_t condition_index) const;

	///\brief Estimates the age at which a stopping condition with
	///derivative information crossed zero.
	///
	///If no zero-crossing is indicated, Inf is returned.
	double crossing_from_history(
			///The index of the condition to search.
			size_t condition_index) const;

	///\brief Initializes the skip_history_zerocrossing and
	///skip_history_extremum arrays appropriately after a mode change.
	void initialize_skip_history(const StoppingCondition &stop_cond,
			StoppingConditionType stopp_reason);

	///\brief Updates the skip_history_zerocrossing and
	///skip_history_extremum arrays appropriately after an acceptable step.
	void update_skip_history(
			const std::valarray<double> &current_stop_cond,
			const StoppingCondition &stop_cond,
			const StopInformation &stop_info);

	///\brief Return true iff the step with the given stop information is
	///acceptable.
	bool acceptable_step(double age, const StopInformation &stop_info);

	///\brief Updates stop_cond_history and stop_deriv_history after a GSL
	///step, returning if/where the evolution needs to stop.
	///
	///Returns the full information about the closest estimated age where a
	///condition is zero or has an extremum exists which might have crossed
	///zero.
	StopInformation update_stop_condition_history(
			///System age in Gyr.
			double age,
			
			///The values of the current evolution variables.
			const std::valarray<double> &orbit,

			///The rates of change of the evolution variables per Gyr.
			const std::valarray<double> &derivatives,

			///The stellar system being evolved.
			const StellarSystem &system,
			
			///The current evolution mode.
			EvolModeType evolution_mode,

			///The conditions indicating where evolution should stop.
			const StoppingCondition &stop_cond,

			///For the first call of this function for an evolution stretch,
			///this should indicate the reason why the previous stretch was
			///stopped. For subsequent calls during the same evolution
			///stretch it should be NO_STOP.
			StoppingConditionType stop_reason=NO_STOP);

	///\brief Evolves a system until either some age cut-off is reached or
	///some stopping condition crosses zero.
	///
	///
	///Appends each accepted step to tabulated_ages,
	///tabulated_evolution_mode, tabulated_orbit and tabulated_deriv.
	///
	///The return value is true if the last step finished after the stopping
	///condition crossed zero and false if it ended before that.
	StopInformation evolve_until(
			///The planet-star system to evolve.
			StellarSystem &system,
			
			///The age to start this part of the evolution from in Gyr.
			double start_age,

			///The age at which to stop this part of the evolution. On
			///exit, it is overwritten with the age of the last accepted
			///step.
			double &max_age,
			
			///The initial conditions. The contents depends on the value of
			///evolution_mode. See #StellarSystem.differential equations for
			///details.			
			///
			///On exit, it is overwritten with the orbit of the last
			///accepted step.
			std::valarray<double> &orbit,
		
			///On input should be the reason why the last evolution stopped.
			///It should be NO_STOP if this is the first piece of evolution
			///being calculated. On exit it is overwritten with the value
			///appropriate for the next run.
			StoppingConditionType &stop_reason,

			///The maximum step that GSL is allowed to take.
			double max_step,

			///The evolution mode for this part of the evolution.
			EvolModeType evolution_mode,

			///Whether some harmonic of the stellar spin is locked to some
			///harmonic of the orbital frequency and which harmonics.
			const SpinOrbitLockInfo &star_lock,

			///The wind saturation state for this part of the evolution.
			WindSaturationState wind_state,

			///The conditions indicating where evolution should stop.
			///Each sub-condition should be a smooth function of age.
			const StoppingCondition &stop_cond=NoStopCondition(),

			///The semimajor axis where the planet should first appear.
			double planet_formation_semimajor=NaN,

			///The inclination with which the planet should first appear.
			double planet_formation_inclination=NaN);

	///Returns the stopping conditions which end the given evolution mode.
	CombinedStoppingCondition *get_stopping_condition(
				///The evolution mode we need the stopping conditions for.
				EvolModeType evolution_mode,

				///Whether a harmonic of the stellar spin is locked to a
				///harmonic of the orbital frequency and which harmonics.
				const SpinOrbitLockInfo &star_lock) const;

	///\brief Returns the evolution mode that the system is entering,
	///assuming that some critical age is reached (e.g. the disk dissipated).
	EvolModeType critical_age_evol_mode(
			///The critical age in Gyr when something happens.
			double age, 

			///The last orbital state in the old evolution mode.
			const std::valarray<double> &orbit,
			
			///The semimajor at which the planet starts when
			///planet_formation_age is reached.
			double initial_semimajor,
			
			///The planet-star system being evolved.
			const StellarSystem &system, 

			///The old evolution mode.
			EvolModeType evolution_mode,

			///Whether the planet is in an orbit which matches the rotation
			///of the star and which harmonic to which.
			SpinOrbitLockInfo &star_lock,
			
			///The age at which the planet forms.
			double planet_formation_age);

	///\brief Returns the evolution mode that the system is entering and
	///updates the lock.
	///
	///assuming that
	///the last orbital state is orbit (in the old evolution mode), the
	///semimajor at which the planet starts after the disk dissipates is
	///initial_semimajor (in AU) and the previous mode was evolution_mode.
	EvolModeType next_evol_mode(
			///The age of the last step in the old mode in Gyr.
			double age, 

			///The values of the old evolution mode variables being evolved.
			const std::valarray<double> &parameters,

			///The semimajor axis at which the planet first appears.
			double initial_semimajor,
			
			///The planet-star system being evolved.
			StellarSystem &system, 

			///The old evolution mode.
			EvolModeType evolution_mode,

			///Whether the planet is in an orbit which matches the rotation
			///of the star and which harmonic to which.
			SpinOrbitLockInfo &star_lock,

			///The current wind saturation state.
			WindSaturationState wind_state,
			
			///The reason for stopping the evolution.
			StoppingConditionType condition_type,

			///The value of the condition which caused the stop at the last
			///step in the old mode.
			double condition_value,
			
			///Was the last step before the condition actually crossed zero.
			bool stopped_before,
			
			///The age at which the planet suddenly appears in Gyr.
			double planet_formation_age);

	///\brief The age at which the evolution should stop next if no other
	///stopping condition occurs.
	double stopping_age(
			///The age from which the next part of the evolution starts.
			double age, 
			
			///The current evolution mode.
			EvolModeType evolution_mode,

			///The stellar system being evolved.
			const StellarSystem &system,
			
			///The age at which the planet forms.
			double planet_formation_age,

			///A sorted list of ages which must be stopped at.
			const std::list<double> &required_ages);

	///\brief Converts the entries in an input orbit or its rate of evolution
	///to a standard set of parameters.
	void parse_orbit_or_derivatives(
			///The mode the input orbit is in. If the mode is such that the
			///planet is not present a and inclination are set to NaN.
			EvolModeType evolution_mode,

			///Does the input orbit represent spin-orbit locked star.
			const SpinOrbitLockInfo &star_lock,

			///The age at which both the input and output orbits apply
			double age,

			///The current orbit or derivatives
			const std::valarray<double> &orbit_deriv,
			
			///The stellar system being evolved.
			const StellarSystem &system,

			///Are we transforming an orbit or its evolution rates.
			bool evolution,
			
			///If transforming derivatives on input should contain the
			///semimajor axis to assume.
			///On output contains the semimajor axis or its derivative.
			double &a,

			///On output contains the inclination or its derivative.
			double &inclination,

			///On output contains the stellar convective zone angular
			///momentum or its derivative.
			double &Lconv,

			///On output contains the stellar radiative zone angular
			///momentum component along Lconv or its derivative.
			double &Lrad_parallel,

			///On output contains the stellar radiative zone angular
			///momentum component perpendicular to Lconv or its derivative.
			double &Lrad_perpendicular) const;

	///\brief Generates an orbit array with the appropriate content for
	///running evolution calculations.
	void collect_orbit_or_derivatives(
			///The evolution mode of the desired orbit.
			EvolModeType evolution_mode,

			///Should the output orbit represent spin-orbit locked star.
			const SpinOrbitLockInfo &star_lock,

			///The semimajor axis or its derivative.
			double a,

			///The inclination or its derivative.
			double inclination,

			///The stellar convective zone angular momentum or its
			///derivative.
			double Lconv,

			///The stellar radiative zone angular momentum component along
			///Lconv or its derivative.
			double Lrad_parallel,

			///The stellar radiative zone angular momentum component
			///perpendicular to Lconv or its derivative.
			double Lrad_perpendicular,

			///The array to fill with the orbit. Resized as necessary.
			std::valarray<double> &result,

			///If we are collecitng derivatives, this should be the value of
			///the semimajor axis. Othrewise it should be NaN.
			double semimajor=NaN) const;

	///Transforms orbital parameters from one evolution mode to another.
	std::valarray<double> transform_orbit(
			///The mode the input orbit is in.
			EvolModeType from_mode,

			///The mode to transform to.
			EvolModeType to_mode, 
			
			///Does the input orbit represent spin-orbit locked star.
			const SpinOrbitLockInfo &from_star_lock,

			///Should the output orbit represent spin-orbit locked star.
			const SpinOrbitLockInfo &to_star_lock,
			
			///The age at which both the input and output orbits apply
			double age,

			///The current orbit
			const std::valarray<double> &from_orbit,

			///The semimajor axis at which the planet first appears in AU.
			double initial_semimajor,

			///The inclination with which the planet first appears in AU.
			double initial_inclination,
			
			///The stellar system being evolved.
			const StellarSystem &system) const;

	///\brief Transforms the deriatives of the orbital parameters from one
	///evolution mode to another.
	std::valarray<double> transform_derivatives(
			///The mode in which the input orbit and derivatives are.
			EvolModeType from_mode,

			///The mode which we need the output derivatives to be in.
			EvolModeType to_mode, 
			
			///Does the input orbit represent spin-orbit locked star.
			const SpinOrbitLockInfo &from_star_lock,

			///Should the output orbit represent spin-orbit locked star.
			const SpinOrbitLockInfo &to_star_lock,
			
			double age,
			const std::valarray<double> &from_orbit, 
			const std::valarray<double> &from_deriv,
			double initial_semimajor, const StellarSystem &system);

	///Clears any previously calculated evolution.
	void reset();
public:
	///\brief Prepare to solve for the orbital evolution.
	OrbitSolver(
			///The end age for the evolution.
			double max_age,
			
			///The precision which to require of the solution.
			double required_precision);

	///\brief Actually solves the given differential equation with the given
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

			///The inclination with which the planet first appears. This
			///argument is ignored if planet_formation_age<=start_age.
			double planet_formation_inclination=NaN,

			///The age at which to start the evolution. Use NaN (default) to
			///start when the radiative core first starts to appear, in
			///which case, start_orbit and initial_evol_mode should be left
			///at their default values as well, but planet_formation_age and
			///planet_formation_semimajor must be specified.
			double start_age=NaN, 

			///The initial evolution mode of the system
			EvolModeType initial_evol_mode=LOCKED_TO_DISK,

			///Whether we are starting with some harmonic of the star's
			///surface rotation locked to some harmonic of the orbit, and
			///if yes, which harmonics.
			const SpinOrbitLockInfo &initial_lock=SpinOrbitLockInfo(1, 1, 0),

			///The initial state to start the system in. The contents
			///depends on initial_evol_mode.
			const std::valarray<double>
				&start_orbit=std::valarray<double>(0.0, 1),

			///A sorted list of ages to include in the tabulated evolution.
			const std::list<double> &required_ages=std::list<double>(),
			
			///If given, the evolution mode is kept constant regardless
			///of what actually occurs with the evolution
			bool no_evol_mode_change=false);

	///Returns the values of a variable at the tabulated ages.
	const std::list<double>
		*get_tabulated_var(EvolVarType var_type) const;

	///\brief Returns the derivative of an independent variable at the
	///points specified by get_tabulated_indep_var.
	const std::list<double> *get_tabulated_var_deriv(
			EvolVarType var_type) const;

	///Returns a list of the evolution modes.
	const std::list<EvolModeType> *get_tabulated_evolution_mode() const
	{return &__tabulated_evolution_mode;}

	///Returns a list of the wind saturation states.
	const std::list<WindSaturationState> *get_tabulated_wind_state() const
	{return &__tabulated_wind_saturation;}

	///Returns a list of the lock states the system evolved through.
	const std::list<SpinOrbitLockInfo> *get_locks() const
	{return &__tabulated_lock;}
	
};

///\brief Converts a semimajor axis given in AU to the variable used in
///solving for the orbital evolution.
inline double AU_to_Rsun6p5(double au) {
	return std::pow(au*AstroConst::AU/AstroConst::solar_radius, 6.5);
};

///\brief Converts the variable used in solving for the orbital evolution to
///a semimajor axis given in AU.
inline double Rsun6p5_to_AU(double rsun) {
  return std::pow(rsun, 1.0/6.5)*AstroConst::solar_radius/AstroConst::AU;
};

#endif
