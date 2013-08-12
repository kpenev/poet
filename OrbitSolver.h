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
public:
	StopInformation(double stop_age=Inf, double stop_precision=NaN,
			StoppingConditionType stop_reason=NO_STOP,
			bool is_crossing=false) :
		__stop_age(stop_age),
		__stop_condition_precision(stop_precision),
		__stop_reason(stop_reason),
		__is_crossing(is_crossing) {}

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
};

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

	///The ages at which the stop condition history is kept
	std::list<double> stop_history_ages;

	std::list< std::valarray<double> >
		stop_cond_history,///< Past values of the stop conditions
		stop_deriv_history;///< Past values of the stop condition derivatives

	///Appends the given orbit and derivatives to tabulated_orbit and
	///tabulated_deriv respectively assuming the orbit contains the
	///quantities evolved for the given evolution mode.
	void append_to_orbit(const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			EvolModeType evolution_mode, double age,
			const StellarSystem &system,
			double planet_formation_semimajor=NaN);

	///Clears the current stopping condition history.
	void clear_stop_history();

	///Returns the dimension of the ODEs governing the evolution of the
	///given type.
	size_t ode_dimension(EvolModeType evolution_mode);

	///Estimates the value and age of an extremum if it potentially can cross
	///a zero by using the last and past tabulated points. If no extremum is
	///indicated by the points, or if it is in the wrong direction, returns
	///the result of the default constructor of ExtremumInformation.
	ExtremumInformation extremum_from_history(double age,
			double stop_condition_value, size_t condition_index)
		const;

	///Estimates the value and age of an extremum if it potentially can cross
	///a zero by using the last and past tabulated points. If no extremum is
	///indicated by the points, or if it is in the wrong direction, returns
	///the result of the default constructor of ExtremumInformation.
	ExtremumInformation extremum_from_history(double age,
			double stop_condition_value, double stop_condition_derivative,
			size_t condition_index) const;

	///Estimates the age at which the stopping condition with the given index
	///crossed zero. If no zero-crossing is indicated, Inf is returned.
	double crossing_from_history(double age, double stop_condition_value,
			size_t condition_index) const;

	///Estimates the age at which the stopping condition with the given index
	///crossed zero. If no zero-crossing is indicated, Inf is returned.
	double crossing_from_history(double age, double stop_condition_value,
			double stop_condition_derivative, size_t condition_index) const;

	///This function is called after every GSL step. It updates the
	///stop_cond_history and stop_deriv_history variables appropriately.
	///Returns the full information about the closest estimated age where a
	///condition is zero or an extremum exists which might have crossed zero.
	StopInformation update_stop_condition_history(double age,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			const StellarSystem &system, EvolModeType evolution_mode,
			const StoppingCondition &stop_cond,
			StoppingConditionType stopped_before=NO_STOP);

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
	///The stopped before argument on input should be true if the previous 
	///evolve_until stopped before the stopping condition changed sign, and
	///on exit it is overwritten with the value appropriate for the next run
	void evolve_until(StellarSystem *system, double start_age,
			double &max_age, std::valarray<double> &orbit,
			StoppingConditionType &stop_reason, double &stop_condition_value,
			StoppingConditionType &stopped_before, double max_step=Inf,
			EvolModeType evolution_mode=FAST_PLANET,
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
				&start_orbit=std::valarray<double>(0.0, 1));


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
