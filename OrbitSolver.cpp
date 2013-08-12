#include "OrbitSolver.h"
#include <iostream>
#include <iomanip>

const double Rsun_AU=AstroConst::solar_radius/AstroConst::AU;

///More civilized output for EvolModeType variables.
std::ostream &operator<<(std::ostream &os, const EvolModeType &evol_mode)
{
	switch(evol_mode) {
		case FAST_PLANET: os << "FAST_PLANET"; break;
		case LOCKED_TO_PLANET: os << "LOCKED_TO_PLANET"; break;
		case SLOW_PLANET: os << "SLOW_PLANET"; break;
		case NO_PLANET: os << "NO_PLANET"; break;
		case LOCKED_TO_DISK: os << "LOCKED_TO_DISK"; break;
		case TABULATION : os << "TABULATION";
	}
	return os;
}

///More civilized output for StoppingConditionType variables.
std::ostream &operator<<(std::ostream &os,
		const StoppingConditionType &stop_cond_type)
{
	switch(stop_cond_type) {
		case NO_STOP: os << "NO_STOP"; break;
		case SYNCHRONIZED: os << "SYNCHRONIZED"; break;
		case BREAK_LOCK: os << "BREAK_LOCK"; break;
		case PLANET_DEATH: os << "PLANET_DEATH"; break;
		case WIND_SATURATION: os << "WIND_SATURATION"; break;
		case ROT_FAST: os << "ROT_FAST";
	}
	return os;
}

///More civilized output for EvolVarType variables.
std::ostream &operator<<(std::ostream &os, const EvolVarType &evol_var)
{
	switch(evol_var) {
		case AGE : os << "AGE"; break;
		case SEMIMAJOR : os << "SEMIMAJOR"; break;
		case LCONV : os << "LCONV"; break;
		case LRAD : os << "LRAD";
	};
	return os;
}

///Returns the difference between the orbital and stellar spin angular
///velocities divided by the orbital angular velocity.
std::valarray<double> SynchronizedCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double a=(evol_mode==LOCKED_TO_DISK ? __initial_semimajor :
			std::pow(orbit[0], 1.0/6.5)*Rsun_AU),
		   worb=system.get_planet().orbital_angular_velocity_semimajor(a),
		   Lconv=(evol_mode==LOCKED_TO_DISK ?
				   system.get_star().get_disk_lock_frequency()*
				   system.get_star().moment_of_inertia(age, convective) :
				   orbit[1]),
		   wconv=system.get_star().spin_frequency(age, convective, Lconv);

	if(derivatives.size()>0) {
		double da_dt=(evol_mode==LOCKED_TO_DISK ? 0 :
				derivatives[0]/(6.5*std::pow(orbit[0],5.5/6.5))*Rsun_AU),
			   dworb_da=system.get_planet().
				   orbital_angular_velocity_semimajor_deriv(a),
			   dworb_dt=dworb_da*da_dt,
			   dwconv_dLconv=system.get_star().spin_frequency_angmom_deriv(
					   age, convective, Lconv),
			   dLconv_dt=(evol_mode==LOCKED_TO_DISK ? 
					   worb*system.get_star().moment_of_inertia_deriv(age,
						   convective) : derivatives[1]);
		stop_deriv.resize(1, (dworb_dt - dwconv_dLconv*dLconv_dt)/worb -
				(worb-wconv)/(worb*worb)*dworb_dt);
	} else stop_deriv.resize(1, NaN);
	return std::valarray<double>((worb-wconv)/worb, 1);
}

double no_planet_dwconv_dt(double age, const std::valarray<double> &orbit,
		const StellarSystem &system)
{
	double no_planet_diff_eq[2], no_planet_orbit[2];
	double semimajor=orbit[0]*Rsun_AU,
		   worb=system.get_planet().orbital_angular_velocity_semimajor(
				   semimajor),
		   Iconv=system.get_star().moment_of_inertia(age, convective),
		   dIconv_dt=system.get_star().moment_of_inertia_deriv(age,
				   convective);
	no_planet_orbit[0]=worb*Iconv;
	no_planet_orbit[1]=orbit[1];
	system.no_planet_differential_equation(age, no_planet_orbit,
			no_planet_diff_eq);
	return (no_planet_diff_eq[0] - dIconv_dt*worb)/Iconv;
}

///Returns the difference between the maximum tidal evolution of the
///semimajor axis and the evolution required to keep the star spinning
///synchronously with the orbit divided by the evolution required to keep
///the lock.
std::valarray<double> BreakLockCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system, 
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double semimajor=orbit[0]*Rsun_AU;
	double dwconv_dt=no_planet_dwconv_dt(age, orbit, system),
		   max_semi_evol = system.get_planet().tidal_decay(age, semimajor,
				   (dwconv_dt>0 ? Inf : -Inf), false)/Rsun_AU;
	stop_deriv.resize(1, NaN);
	if(orbit[0]==0) *static_cast<int*>(NULL)=0;
	if(derivatives[0]==0) return std::valarray<double>(Inf, 1);
	return std::valarray<double>(
			(max_semi_evol - derivatives[0])/derivatives[0], 1);
}

///Returns the difference between the semimajor axis and the larger of
///the roche radius and the stellar radius divided by the latter.
std::valarray<double> PlanetDeathCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double semimajor=orbit[0],
		   min_semimajor=system.get_planet().minimum_semimajor(age);

	if(evol_mode==SLOW_PLANET || evol_mode==FAST_PLANET) {
		semimajor*=std::pow(Rsun_AU, 6.5);
		min_semimajor=std::pow(min_semimajor, 6.5);
	} else {
		assert(evol_mode==LOCKED_TO_PLANET);
		semimajor*=Rsun_AU;
	}
	stop_deriv.resize(1, NaN);
	if(std::isinf(min_semimajor)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>((semimajor-min_semimajor)/min_semimajor, 1);
}

///Returns the difference between the convective angular velocity and the
///wind saturation angular velocity divided by the wind saturation
///angular velocity.
std::valarray<double> WindSaturationCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double wconv, wsat=system.get_star().get_wind_saturation_frequency();
	switch(evol_mode) {
		case FAST_PLANET: case SLOW_PLANET: case TABULATION:
			wconv=system.get_star().spin_frequency(age, convective,
					orbit[1]);
			break;
		case LOCKED_TO_PLANET:
			wconv=system.get_planet().orbital_angular_velocity_semimajor(
					orbit[0]*Rsun_AU);
			break;
		case NO_PLANET:
			wconv=system.get_star().spin_frequency(age, convective,
					orbit[0]);
			break;
		case LOCKED_TO_DISK:
			std::cout << "BAD FUNCTION ARG17S" << std::endl;
			throw Error::BadFunctionArguments("Wind saturation condition is "
					"not defined for disk locked evolution mode.");
		default:
			throw Error::BadFunctionArguments(
					"Unrecognized wind saturation argument");
	}
	stop_deriv.resize(1, NaN);
	return std::valarray<double>((wconv-wsat)/wsat, 1);
}

///Returns the difference between the convective angular velocity and the
//threshold velocity
std::valarray<double> RotFastCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double wconv = system.get_star().get_wind_saturation_frequency();
	switch(evol_mode) {
		case FAST_PLANET: case SLOW_PLANET: case TABULATION:
			wconv=system.get_star().spin_frequency(age, convective,
					orbit[1]);
			break;
		case LOCKED_TO_PLANET:
			wconv=system.get_planet().orbital_angular_velocity_semimajor(
					orbit[0]*Rsun_AU);
			break;
		case NO_PLANET:
			wconv=system.get_star().spin_frequency(age, convective,
					orbit[0]);
			break;
		case LOCKED_TO_DISK:
			std::cout << "BAD FUNCTION ARGS" << std::endl;
			throw Error::BadFunctionArguments("RotFast condition is "
					"not defined for disk locked evolution mode.");
	}
	stop_deriv.resize(1, NaN);
	return std::valarray<double>((wconv-spin_thres)/spin_thres, 1);
}

///Adds the conditions in RHS to the conditions of *this.
CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const CombinedStoppingCondition &rhs)
{
	for(std::vector<const StoppingCondition *>::const_iterator
			i=rhs.sub_conditions.begin(); i!=rhs.sub_conditions.end(); i++)
		sub_conditions.push_back(*i);
	return *this;
}

///Adds RHS to the conditions of *this.
CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const StoppingCondition *rhs)
{
	sub_conditions.push_back(rhs);
	return *this;
}

///Deletes all subconditions, unless no_delete_subcond has been
///previously called.
CombinedStoppingCondition::~CombinedStoppingCondition()
{
	for(std::vector<const StoppingCondition *>::const_iterator
			i=sub_conditions.begin(); i!=sub_conditions.end(); i++)
		delete *i;
}

///Returns the values of all stopping sub_conditions
std::valarray<double> CombinedStoppingCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	std::valarray<double> result(sub_conditions.size());
	stop_deriv.resize(sub_conditions.size(), NaN);
	size_t i=0;
	for(std::vector<const StoppingCondition *>::const_iterator
			cond=sub_conditions.begin(); cond!=sub_conditions.end(); cond++) {
		std::valarray<double> sub_stop_deriv(1),
			temp_array=(**cond)(age, orbit, derivatives, system,
					sub_stop_deriv, evol_mode);
		result[i]=temp_array[0];
		stop_deriv[i]=sub_stop_deriv[0];
		i++;
	}
	return result;
}

///A wrapper that allows the stellar system differential equation to be passed
///to the GSL ODE solver.
int stellar_system_diff_eq(double age, const double *orbital_parameters,
		double *orbital_derivatives, void *system_mode_windstate)
{
	void **input_params=static_cast<void **>(system_mode_windstate);
	StellarSystem *system=static_cast< StellarSystem* >(input_params[0]);
	EvolModeType evol_mode=*static_cast< EvolModeType* >(input_params[1]);
	WindSaturationState wind_state=*static_cast<WindSaturationState*>(
			input_params[2]);

	int result;
	switch(evol_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			result=system->orbit_differential_equation(
					age, orbital_parameters, orbital_derivatives, evol_mode,
					wind_state);
			break;
		case LOCKED_TO_PLANET :
			result=system->locked_orbit_differential_equation(
					age, orbital_parameters, orbital_derivatives,wind_state);
			break;
		case NO_PLANET :
			result=system->no_planet_differential_equation(
					age, orbital_parameters, orbital_derivatives,wind_state);
			break;
		case LOCKED_TO_DISK :
			result=system->locked_conv_differential_equation(
					age, orbital_parameters, orbital_derivatives);
			break;
		default :
			throw Error::BadFunctionArguments("Unrecognized evolution mode "
					"encountered in stellar_system_diff_eq!");
	}
	return result;
}

///A wrapper tha allows the stellar system jacobian to be passed
///to the GSL ODE solver.
int stellar_system_jacobian(double age, const double *orbital_parameters,
		double *param_derivs, double *age_derivs,void *system_mode_windstate)
{
	void **input_params=static_cast<void **>(system_mode_windstate);
	StellarSystem *system=static_cast< StellarSystem* >(input_params[0]);
	EvolModeType evol_mode=*static_cast< EvolModeType* >(input_params[1]);
	WindSaturationState wind_state=*static_cast<WindSaturationState*>(
			input_params[2]);

	int result;
	switch(evol_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			result=system->orbit_jacobian(age, orbital_parameters,
					param_derivs, age_derivs, evol_mode, wind_state);
			break;
		case LOCKED_TO_PLANET :
			result=system->locked_orbit_jacobian(age, orbital_parameters,
					param_derivs, age_derivs, wind_state);
			break;
		case NO_PLANET :
			result=system->no_planet_jacobian(age, orbital_parameters,
					param_derivs, age_derivs, wind_state);
			break;
		case LOCKED_TO_DISK :
			result=system->locked_conv_jacobian(age, orbital_parameters,
					param_derivs, age_derivs);
			break;
		default :
			throw Error::BadFunctionArguments("Unrecognized evolution mode "
					"encountered in stellar_system_jacobian!");
	}
	return result;
}

///Appends the given orbit and derivatives to tabulated_orbit and
///tabulated_deriv respectively assuming the orbit contains the
///quantities evolved for the given evolution mode.
void OrbitSolver::append_to_orbit(const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		EvolModeType evolution_mode, double age, const StellarSystem &system,
		double planet_formation_semimajor)
{
	std::valarray<double> orbit_to_tabulate=transform_orbit(evolution_mode,
			TABULATION, age, orbit, planet_formation_semimajor, system),
		deriv_to_tabulate=transform_derivatives(evolution_mode, TABULATION, 
				age, orbit, derivatives, planet_formation_semimajor, system);
	size_t nargs=orbit_to_tabulate.size();
	tabulated_ages.push_back(age);
	for(size_t i=0; i<nargs; i++) {
		tabulated_orbit[i].push_back(orbit_to_tabulate[i]);
		tabulated_deriv[i].push_back(deriv_to_tabulate[i]);
	}
	tabulated_evolution_mode.push_back(evolution_mode);
}

///Clears the current stopping condition history.
void OrbitSolver::clear_stop_history()
{
	stop_history_ages.clear();
	stop_cond_history.clear();
	stop_deriv_history.clear();
}

///Returns the dimension of the ODEs governing the evolution of the
///given type.
size_t OrbitSolver::ode_dimension(EvolModeType evolution_mode)
{
	switch(evolution_mode) {
		case FAST_PLANET : case SLOW_PLANET : return 3;
		case LOCKED_TO_PLANET : case NO_PLANET : return 2;
		case LOCKED_TO_DISK : return 1;
		default :
			throw Error::BadFunctionArguments("Unrecognized evolution "
					"mode in OrbitSolver::ode_dimension");
	}
}

///Estimates the value and age of an extremum if it potentially can cross
///a zero by using the last and past tabulated points. If no extremum is
///indicated by the points, or if it is in the wrong direction, returns
///the result of the default constructor of ExtremumInformation.
ExtremumInformation OrbitSolver::extremum_from_history(double age,
		double stop_condition_value, size_t condition_index) const
{
	ExtremumInformation result;
	if(stop_history_ages.size()<2) return result;
	double prev1_stop_cond=stop_cond_history.back()[condition_index];
	if(std::abs(prev1_stop_cond)>std::abs(stop_condition_value)) return result;
	std::list< std::valarray<double> >::const_reverse_iterator
		stop_val_iter=stop_cond_history.rbegin();
	++stop_val_iter;
	double prev2_stop_cond=(*stop_val_iter)[condition_index];
	if(std::abs(prev2_stop_cond)<std::abs(prev1_stop_cond)) return result;
	std::list<double>::const_reverse_iterator
		stop_age_iter=stop_history_ages.rbegin();
	double prev1_age=*(++stop_age_iter), prev2_age=*(++stop_age_iter);
	if(stop_history_ages.size()==2)
		if((prev1_stop_cond-prev2_stop_cond)*
				(stop_condition_value-prev1_stop_cond)>0) return result;
		else return quadratic_extremum(prev2_age, prev2_stop_cond,
				prev1_age, prev1_stop_cond, age, stop_condition_value,
				&(result.y()));
	double prev3_age=*(++stop_age_iter),
		   prev3_stop_cond=(*(++stop_val_iter))[condition_index];
	if((prev2_stop_cond-prev3_stop_cond)*(prev1_stop_cond-prev2_stop_cond)>0
			&&
			(prev1_stop_cond-prev2_stop_cond)*
			(stop_condition_value-prev1_stop_cond)>0) return result;
	result.x()=cubic_extremum(prev3_age, prev3_stop_cond, prev2_age,
			prev2_stop_cond, prev1_age, prev1_stop_cond, age,
			stop_condition_value, &(result.y()));
	return result;
}

///Estimates the value and age of an extremum if it potentially can cross
///a zero by using the last and past tabulated points. If no extremum is
///indicated by the points, or if it is in the wrong direction, returns
///the result of the default constructor of ExtremumInformation.
ExtremumInformation OrbitSolver::extremum_from_history(double age,
		double stop_condition_value, double stop_condition_derivative,
		size_t condition_index) const
{
	ExtremumInformation result;
	if(stop_history_ages.size()==0) return result;
	double prev_stop_cond=stop_cond_history.back()[condition_index];
	if(stop_condition_value*prev_stop_cond<0) return result;
	if(std::isnan(stop_condition_derivative))
		return extremum_from_history(age, stop_condition_value,
				condition_index);
	double prev_stop_deriv=stop_deriv_history.back()[condition_index];
	if(stop_condition_value*stop_condition_derivative<0 ||
			stop_condition_derivative*prev_stop_deriv>=0) return result;
	result.x()=estimate_extremum(stop_history_ages.back(), prev_stop_cond,
			age, stop_condition_value, prev_stop_deriv,
			stop_condition_derivative, &(result.y()));
	return result;
}

///Estimates the age at which the stopping condition with the given index
///crossed zero. If no zero-crossing is indicated, Inf is returned.
double OrbitSolver::crossing_from_history(double age,
		double stop_condition_value, size_t condition_index) const
{
	std::list< std::valarray<double> >::const_reverse_iterator
		stop_val_iter=stop_cond_history.rbegin();
	double prev1_stop_cond=(*stop_val_iter)[condition_index];
	if(stop_condition_value*prev1_stop_cond>0) return Inf;
	std::list<double>::const_reverse_iterator
		stop_age_iter=stop_history_ages.rbegin();
	double prev1_age=*stop_age_iter;
	if(stop_history_ages.size()==1) return estimate_zerocrossing(
			prev1_age, prev1_stop_cond, age, stop_condition_value);
	double prev2_stop_cond=(*(++stop_val_iter))[condition_index],
		   prev2_age=*(++stop_age_iter);
	if(stop_history_ages.size()==2) return quadratic_zerocrossing(
			prev2_age, prev2_stop_cond, prev1_age, prev1_stop_cond,
			age, stop_condition_value);
	double prev3_stop_cond=(*(++stop_val_iter))[condition_index],
		   prev3_age=*(++stop_age_iter);
	return cubic_zerocrossing(prev3_age, prev3_stop_cond,
			prev2_age, prev2_stop_cond, prev1_age, prev1_stop_cond,
			age, stop_condition_value);
}

///Estimates the age at which the stopping condition with the given index
///crossed zero. If no zero-crossing is indicated, Inf is returned.
double OrbitSolver::crossing_from_history(double age,
		double stop_condition_value, double stop_condition_derivative,
		size_t condition_index) const
{
	if(stop_history_ages.size()==0) return Inf;
	if(std::isnan(stop_condition_derivative)) return crossing_from_history(
			age, stop_condition_value, condition_index);
	double prev_stop_cond=stop_cond_history.back()[condition_index];
	if(stop_condition_value*prev_stop_cond>0) return Inf;
	double prev_stop_deriv=stop_deriv_history.back()[condition_index],
		   prev_age=stop_history_ages.back();
	return estimate_zerocrossing(prev_age, prev_stop_cond, age,
			stop_condition_value, prev_stop_deriv,stop_condition_derivative);
}

///This function is called after every GSL step. It updates the
///stop_cond_history and stop_deriv_history variables appropriately.
///Returns the full information about the closest estimated age where a
///condition is zero or an extremum exists which might have crossed zero.
StopInformation OrbitSolver::update_stop_condition_history(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system, EvolModeType evolution_mode,
		const StoppingCondition &stop_cond,
		StoppingConditionType stopped_before)
{
	std::valarray<double> current_stop_cond(stop_cond.num_subconditions()),
		current_stop_deriv;
	current_stop_cond=stop_cond(age, orbit, derivatives, system,
			current_stop_deriv, evolution_mode);
	StopInformation result;
	for(size_t cond_ind=0; cond_ind<current_stop_cond.size(); cond_ind++) {
		StoppingConditionType
			stop_cond_type=stop_cond.type(cond_ind);
		if((stopped_before==BREAK_LOCK && stop_cond_type==SYNCHRONIZED) || 
				(stopped_before!=NO_STOP && stop_cond_type==stopped_before))
			current_stop_cond[cond_ind]=0;
/*		else if()
			current_stop_cond[cond_ind]*=-1;*/
		double stop_cond_value=current_stop_cond[cond_ind],
			   stop_cond_deriv=current_stop_deriv[cond_ind],
			   crossing_age=crossing_from_history(age, stop_cond_value,
					   stop_cond_deriv, cond_ind);
		ExtremumInformation extremum=extremum_from_history(age,
				stop_cond_value, stop_cond_deriv, cond_ind);
		if(std::min(crossing_age, extremum.x())<result.stop_age()) {
			result.stop_age()=std::min(crossing_age, extremum.x());
			result.stop_reason()=stop_cond_type; 
			if(crossing_age<=extremum.x()) {
				result.is_crossing()=true;
				result.stop_condition_precision()=stop_cond_value;
			} else {
				result.is_crossing()=false;
				result.stop_condition_precision()=
					extremum.y()-stop_cond_value;
			}
		}
	}
	if(result.stop_age()>age) {
		stop_history_ages.push_back(age);
		stop_cond_history.push_back(current_stop_cond);
		stop_deriv_history.push_back(current_stop_deriv);
	}
	return result;
}

void OrbitSolver::evolve_until(StellarSystem *system, double start_age,
		double &max_age, std::valarray<double> &orbit,
		StoppingConditionType &stop_reason, double &stop_condition_value,
		StoppingConditionType &stopped_before, double max_step,
		EvolModeType evolution_mode, WindSaturationState wind_state,
		const StoppingCondition &stop_cond,
		double planet_formation_semimajor)
{
	size_t nargs=orbit.size();

	const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_bsimp;
//	const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;

	gsl_odeiv2_step *step=gsl_odeiv2_step_alloc(step_type, nargs);
	gsl_odeiv2_control *step_control=gsl_odeiv2_control_standard_new(
			precision, precision, 1, 1);
	gsl_odeiv2_evolve *evolve=gsl_odeiv2_evolve_alloc(nargs);

	void *sys_mode_windstate[3]={system, &evolution_mode, &wind_state};
	gsl_odeiv2_system ode_system={stellar_system_diff_eq,
		stellar_system_jacobian, nargs, sys_mode_windstate};
	double t=start_age; 
	std::valarray<double> derivatives(nargs), param_derivatives(nargs),
		age_derivatives(nargs);
	stellar_system_diff_eq(t, &(orbit[0]), &(derivatives[0]),
			sys_mode_windstate);
	append_to_orbit(orbit, derivatives, evolution_mode, t, *system,
			planet_formation_semimajor);

	update_stop_condition_history(t, orbit, derivatives, *system,
			evolution_mode, stop_cond, stopped_before);

	double step_size=0.01*(max_age-start_age);

	while(t<max_age) {
		double max_next_t=std::min(t + max_step, max_age), old_t=t;
		std::valarray<double> old_orbit(orbit), old_derivatives(derivatives);
		int status=GSL_SUCCESS;
		StopInformation stop;
		do {
			status=gsl_odeiv2_evolve_apply(
					evolve, step_control, step, &ode_system,
					&t, max_next_t, &step_size, &(orbit[0]));
			stellar_system_diff_eq(t, &(orbit[0]), &(derivatives[0]),
					sys_mode_windstate);
			stop=update_stop_condition_history(t, orbit, derivatives, *system,
					evolution_mode, stop_cond, stopped_before);
			if(stop.stop_age()<t) {
				orbit=old_orbit;
				derivatives=old_derivatives;
				t=old_t;
				if(stop.is_crossing()) stopped_before=stop.stop_reason();
			} else stopped_before=
				(stop.stop_reason()==BREAK_LOCK ? BREAK_LOCK : NO_STOP);
		} while(stop.stop_reason()!=NO_STOP &&
				std::abs(stop.stop_condition_precision())>precision);
		if (status != GSL_SUCCESS) {
			std::ostringstream msg;
			msg << "GSL signaled failure while evolving (error code " <<
				status << ")";
			throw Error::Runtime(msg.str());
		}
		append_to_orbit(orbit, derivatives, evolution_mode, t, *system,
				planet_formation_semimajor);
		if(stop.is_crossing() && stop.stop_reason()!=NO_STOP) break;
	}
	max_age=t;
	clear_stop_history();

	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(step_control);
	gsl_odeiv2_step_free(step);
}

///Returns the stopping condition which ends the given evolution mode.
CombinedStoppingCondition *OrbitSolver::get_stopping_condition(
		EvolModeType evolution_mode, double initial_semimajor,
		const Planet *planet) const
{
	CombinedStoppingCondition *result=new CombinedStoppingCondition();
	if(evolution_mode==LOCKED_TO_DISK) return result;
	(*result)|=new WindSaturationCondition;
	(*result)|= new RotFastCondition(spin_thres);
	if(evolution_mode==NO_PLANET) return result;
	(*result)|=new PlanetDeathCondition();
	switch(evolution_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			*result|=new SynchronizedCondition(planet, initial_semimajor);
			break;
		case LOCKED_TO_PLANET :
			*result|=new BreakLockCondition();
			break;
		default :
			*result|=new NoStopCondition();
	}
	return result;
}

///Returns the evolution mode that the system is entering, assuming that
///some critical age is reached (e.g. the disk dissipated). The last 
///orbital state should be orbit (in the old evolution mode), the 
///semimajor at which the planet starts after the disk dissipates is
///initial_semimajor (in AU) and the previous mode was evolution_mode.
EvolModeType OrbitSolver::critical_age_evol_mode(double age, 
		const std::valarray<double> &orbit,
		double initial_semimajor, const StellarSystem &system, 
		EvolModeType evolution_mode, double planet_formation_age) const
{
	if(system.get_star().get_disk_dissipation_age()>age)
		return LOCKED_TO_DISK;
	else if(planet_formation_age>age) return NO_PLANET;
	SynchronizedCondition sync_condition(&(system.get_planet()),
			initial_semimajor);
	std::valarray<double> dummy_cond_deriv;
	double in_sync=sync_condition(age, orbit, std::valarray<double>(),system,
			dummy_cond_deriv, evolution_mode)[0];
	if(std::abs(in_sync)<precision) {
		BreakLockCondition locked_cond;
		std::valarray<double> locked_orbit=transform_orbit(evolution_mode,
				LOCKED_TO_PLANET, age, orbit, initial_semimajor, system);
		std::valarray<double> locked_deriv(2);
		WindSaturationState wind_sat=UNKNOWN;
		EvolModeType locked_mode=LOCKED_TO_PLANET;
		void *diff_eq_params[]={
			static_cast<void*>(const_cast<StellarSystem*>(&system)),
			static_cast<void*>(&locked_mode),
			static_cast<void*>(&wind_sat)};
		stellar_system_diff_eq(age, &locked_orbit[0], &locked_deriv[0],
				static_cast<void*>(diff_eq_params));

		if(locked_cond(age, locked_orbit, locked_deriv, system,
					dummy_cond_deriv, evolution_mode)[0]>=0)
			return LOCKED_TO_PLANET;
		else return (no_planet_dwconv_dt(age, locked_orbit, system) > 0 ?
				SLOW_PLANET : FAST_PLANET);
	}
	return (in_sync>0 ? FAST_PLANET : SLOW_PLANET);
}

///Returns the evolution mode that the system is entering, assuming that
///the last orbital state is orbit (in the old evolution mode), the
///semimajor at which the planet starts after the disk dissipates is
///initial_semimajor (in AU) and the previous mode was evolution_mode.
EvolModeType OrbitSolver::next_evol_mode(double age,
		const std::valarray<double> &orbit,
		double initial_semimajor, const StellarSystem &system,
		EvolModeType evolution_mode,
		StoppingConditionType condition_type,
		double condition_value, bool stopped_before,
		double planet_formation_age) const
{
	if(condition_type==NO_STOP) return critical_age_evol_mode(age, orbit,
			initial_semimajor, system, evolution_mode, planet_formation_age);
	if(condition_type==SYNCHRONIZED) {
		std::valarray<double> locked_orbit=transform_orbit(evolution_mode,
				LOCKED_TO_PLANET, age, orbit, initial_semimajor, system);
		std::valarray<double> locked_orbit_diff_eq(2), stop_deriv;
		WindSaturationState wind_sat=UNKNOWN;
		EvolModeType locked_mode=LOCKED_TO_PLANET;
		void *diff_eq_params[]={
			static_cast<void*>(const_cast<StellarSystem*>(&system)),
			static_cast<void*>(&locked_mode),
			static_cast<void*>(&wind_sat)};
		stellar_system_diff_eq(age, &locked_orbit[0], &locked_orbit_diff_eq[0],
				static_cast<void*>(diff_eq_params));
		if(BreakLockCondition()(age, locked_orbit, locked_orbit_diff_eq,
					system, stop_deriv, evolution_mode)[0]>=0)
			return LOCKED_TO_PLANET;
		else return (condition_value>0 ? FAST_PLANET : SLOW_PLANET);
	} else if(condition_type==BREAK_LOCK) {
		return (no_planet_dwconv_dt(age, orbit, system)>0 ? SLOW_PLANET :
				FAST_PLANET);
	} else if(condition_type==PLANET_DEATH) return NO_PLANET;
	else if(condition_type==WIND_SATURATION || condition_type==ROT_FAST)
		return evolution_mode;
	else {
		std::ostringstream msg;
		msg << "Invalid stopping condition type encountered in "
			"OrbitSolver::next_evol_mode: " << condition_type;
		throw Error::BadFunctionArguments(msg.str());
	}
}

///Returns what age the evolution with the given mode should stop if
///no other stopping condition occurs.
double OrbitSolver::stopping_age(double age, EvolModeType evolution_mode,
		const StellarSystem &system, double planet_formation_age)
{
	double result;
	switch(evolution_mode) {
		case FAST_PLANET : case LOCKED_TO_PLANET : case SLOW_PLANET :
			result=end_age; break;
		case NO_PLANET :
			result=(std::isnan(planet_formation_age) || 
					age>planet_formation_age ? end_age :
					std::min(planet_formation_age, end_age));
			break;
		case LOCKED_TO_DISK :
			result=std::min(system.get_star().get_disk_dissipation_age(),
					end_age); break;
		default :
			throw Error::BadFunctionArguments("Unrecognized evolution "
					"mode in OrbitSolver::stopping_age");
	}
	double core_formation_age=system.get_star().core_formation_age();
	if(age < core_formation_age)
		result = std::min(result, core_formation_age);
	//we want 1 data point right at main sequence start
	if (age < main_seq_start)
		result = std::min(result, main_seq_start);
	return result;
}

///Transforms orbital parameters from one evolution mode (from_mode) to
///another (to_mode).
std::valarray<double> OrbitSolver::transform_orbit(EvolModeType from_mode,
		EvolModeType to_mode, double age,
		const std::valarray<double> &from_orbit, 
		double initial_semimajor, const StellarSystem &system) const
{
	double a, Lconv, Lrad;
	switch(from_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			assert(from_orbit.size()==3);
			a=std::pow(from_orbit[0], 1.0/6.5);
			Lconv=from_orbit[1]; Lrad=from_orbit[2];
			break;
		case LOCKED_TO_PLANET :
			assert(from_orbit.size()==2);
			a=from_orbit[0]; Lrad=from_orbit[1];
			Lconv=system.get_planet().orbital_angular_velocity_semimajor(
					a*Rsun_AU)*
				system.get_star().moment_of_inertia(age, convective);
			break;
		case NO_PLANET :
			Lconv=from_orbit[0]; Lrad=from_orbit[1];
			if (to_mode != TABULATION)
				a = initial_semimajor/Rsun_AU;
			else
				a = NaN;
			break;
		case LOCKED_TO_DISK :
			assert(from_orbit.size()==1);
			a=(to_mode==TABULATION ? NaN : initial_semimajor/Rsun_AU);
			Lconv=system.get_star().get_disk_lock_frequency()*
				system.get_star().moment_of_inertia(age, convective);
			Lrad=from_orbit[0];
			break;
		default :
			throw Error::BadFunctionArguments("Urecognized evolution mode "
					"(from_mode) in transform_orbit.");
	}
	std::valarray<double> to_orbit;
	switch(to_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			to_orbit.resize(3);
			to_orbit[0]=std::pow(a, 6.5); to_orbit[1]=Lconv; to_orbit[2]=Lrad;
			break;
		case LOCKED_TO_PLANET :
			to_orbit.resize(2);
			to_orbit[0]=a; to_orbit[1]=Lrad;
			break;
		case NO_PLANET :
			to_orbit.resize(2);
			to_orbit[0]=Lconv+(from_mode==FAST_PLANET || 
					from_mode==LOCKED_TO_PLANET || 
					from_mode==SLOW_PLANET ?
					system.get_planet().orbital_angular_momentum(a*Rsun_AU) :
					0.0);
			to_orbit[1]=Lrad;
			break;
		case LOCKED_TO_DISK :
			to_orbit.resize(1, Lrad);
			break;
		case TABULATION:
			to_orbit.resize(3);
			to_orbit[0]=a; to_orbit[1]=Lconv; to_orbit[2]=Lrad;
	}
	return to_orbit;
}

///Transforms the deriatives of the orbital parameters from one evolution
///mode (from_mode) to another (to_mode).
std::valarray<double> OrbitSolver::transform_derivatives(
		EvolModeType from_mode, EvolModeType to_mode, double age,
		const std::valarray<double> &from_orbit, 
		const std::valarray<double> &from_deriv,
		double initial_semimajor, const StellarSystem &system)
{
	if(from_mode==NO_PLANET && to_mode!=TABULATION) 
		throw Error::BadFunctionArguments(
				"No evolution mode can possibly follow NO_PLANET in "
				"transform_derivatives.");

	double a, da_dt, dLconv_dt, dLrad_dt;
	assert(from_orbit.size()==from_deriv.size());
	switch(from_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			assert(from_orbit.size()==3);
			a=std::pow(from_orbit[0], 1.0/6.5);
			da_dt=1.0/6.5*std::pow(from_orbit[0], -5.5/6.5)*from_deriv[0];
			dLconv_dt=from_deriv[1]; dLrad_dt=from_deriv[2];
			break;
		case LOCKED_TO_PLANET :
			assert(from_orbit.size()==2);
			a=from_orbit[0];
			da_dt=from_deriv[0]; dLrad_dt=from_deriv[1];
			dLconv_dt=
				system.get_planet().orbital_angular_velocity_semimajor_deriv(
						a*Rsun_AU)*
				system.get_star().moment_of_inertia(age, convective)*
				da_dt*Rsun_AU
				+
				system.get_planet().orbital_angular_velocity_semimajor(
					a*Rsun_AU)*
				system.get_star().moment_of_inertia_deriv(age, convective);
			break;
		case NO_PLANET :
			a=NaN; da_dt=NaN;
			dLconv_dt=from_deriv[0]; dLrad_dt=from_deriv[1];
			break;
		case LOCKED_TO_DISK :
			assert(from_orbit.size()==1);
			a=(to_mode==TABULATION ? NaN : initial_semimajor/Rsun_AU);
			da_dt=NaN;
			dLconv_dt=system.get_star().get_disk_lock_frequency()*
				system.get_star().moment_of_inertia_deriv(age, convective);
			dLrad_dt=from_deriv[0];
			break;
		default :
			throw Error::BadFunctionArguments("Urecognized evolution mode "
					"(from_mode) in transform_orbit.");
	}
	std::valarray<double> to_deriv;
	switch(to_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			to_deriv.resize(3);
			to_deriv[0]=6.5*std::pow(a, 5.5)*da_dt;
			to_deriv[1]=dLconv_dt; to_deriv[2]=dLrad_dt;
			break;
		case LOCKED_TO_PLANET :
			to_deriv.resize(2);
			to_deriv[0]=da_dt; to_deriv[1]=dLrad_dt;
			break;
		case NO_PLANET :
			to_deriv.resize(2);
			to_deriv[0]=dLconv_dt; to_deriv[1]=dLrad_dt;
			break;
		case LOCKED_TO_DISK :
			to_deriv.resize(1, dLrad_dt);
			break;
		case TABULATION:
			to_deriv.resize(3);
			to_deriv[0]=da_dt; to_deriv[1]=dLconv_dt; to_deriv[2]=dLrad_dt;
	}
	return to_deriv;
}

///Clears any previously calculated evolution.
void OrbitSolver::reset()
{
	tabulated_ages.clear();
	tabulated_evolution_mode.clear();
	for (int i=0; i < 3; i++) {
		tabulated_orbit[i].clear();
		tabulated_deriv[i].clear();
	}
}

OrbitSolver::OrbitSolver(
		double min_age, double max_age, double required_precision,
		double spin_thres, double main_seq_start) :
	start_age(min_age), end_age(max_age), precision(required_precision),
	spin_thres(spin_thres), main_seq_start(main_seq_start),
	tabulated_orbit(3), tabulated_deriv(3)
{
	/* Prepare to solve for the orbital evolution over the given range
	and to the required precision. */
	if (max_age > MAX_END_AGE) end_age = MAX_END_AGE;
}

///Actually solves the given differential equation with the given
///boundary conditions.
void OrbitSolver::operator()(
		///The stellar system to calculate the evolution for
		StellarSystem &system,

		///The maximum size of the time steps allowed (useful if finer
		///sampling of the output than default is necessary).
		double max_step,

		///The age at which a planet magically appears in a perfectly
		///circularized orbit and starts affecting the system. By
		///default, the planet is assumed to always be there and
		///planet_formation_semimajor has no effect.
		double planet_formation_age,

		///The semimajor axis at which the planet first appears. This
		///argument is ignore if planet_formation_age<=start_age.
		double planet_formation_semimajor,

		///The age at which to start the evolution. Use NaN (default) to
		///start when the radiative core first starts to appear, in
		///which case, start_orbit and initial_evol_mode should be left
		///at their default values as well, but planet_formation_age and
		///planet_formation_semimajor must be specified.
		double start_age,

		///The initial evolution mode of the system
		EvolModeType initial_evol_mode,

		///The initial state to start the system in. The contents
		///depends on initial_evol_mode.
		const std::valarray<double> &start_orbit)
{
	if(std::isnan(start_age))
		start_age=system.get_star().core_formation_age();
	double wconv;
	switch(initial_evol_mode) {
		case FAST_PLANET : case SLOW_PLANET :
			assert(start_orbit.size()==3); 
			wconv=system.get_star().spin_frequency(start_age, convective,
					start_orbit[1]);
			break;
		case LOCKED_TO_PLANET :
			assert(start_orbit.size()==2);
			wconv=system.get_planet().orbital_angular_velocity_semimajor(
					start_orbit[0]*Rsun_AU);
			break;
		case NO_PLANET :
			assert(start_orbit.size()==2);
			wconv=system.get_star().spin_frequency(start_age, convective,
					start_orbit[0]);
			break;
		case LOCKED_TO_DISK :
			assert(start_orbit.size()==1);
			wconv=system.get_star().get_disk_lock_frequency();
			break;
		default:
			throw Error::BadFunctionArguments("Unsupported initial evolution"
					" mode encountered in OrbitSolver::operator().");
	}
	WindSaturationState wind_state=(
			wconv<system.get_star().get_wind_saturation_frequency() ?
			NOT_SATURATED : SATURATED);

	reset();
	StoppingConditionType stopped_before=NO_STOP;
	double last_age=start_age;
	std::valarray<double> orbit=start_orbit;
	EvolModeType evolution_mode=initial_evol_mode;
	while(last_age<end_age) {
		if (evolution_mode == FAST_PLANET || evolution_mode == SLOW_PLANET ||
				evolution_mode == LOCKED_TO_PLANET) {
			//if planet starts out below minimum semimajor axis,
			//PlanetDeathCondition never reaches 0, so we have to manually
			//check for it
			PlanetDeathCondition deathTest;
			std::valarray<double> fakeArr;
			double deathResult = deathTest(last_age, orbit, fakeArr,
					system, fakeArr, evolution_mode)[0];
			if (deathResult < 0) return;
		}
		//std::cout << deathResult << std::endl;
		double stop_evol_age=stopping_age(last_age, evolution_mode,
				system, planet_formation_age);
		CombinedStoppingCondition *stopping_condition=get_stopping_condition(
				evolution_mode, planet_formation_semimajor,
				&system.get_planet());
		StoppingConditionType stop_reason = NO_STOP;
		double stop_condition_value;
		evolve_until(&system, start_age, stop_evol_age, orbit, stop_reason,
				stop_condition_value, stopped_before, max_step,
				evolution_mode, wind_state, *stopping_condition,
				planet_formation_semimajor);
		last_age=stop_evol_age;
		if(last_age<end_age) {
			EvolModeType old_evolution_mode=evolution_mode;
			evolution_mode=next_evol_mode(last_age, orbit,
					planet_formation_semimajor, system, evolution_mode,
					stop_reason, stop_condition_value, stopped_before,
					planet_formation_age);
			if(old_evolution_mode!=evolution_mode) {
				std::valarray<double> new_orbit=transform_orbit(
						old_evolution_mode, evolution_mode, last_age, orbit,
						planet_formation_semimajor, system);
				orbit.resize(new_orbit.size());
				orbit=new_orbit;
			}
			if(stop_reason==WIND_SATURATION)
				wind_state=static_cast<WindSaturationState>(-wind_state);
			start_age=last_age;
		}
		delete stopping_condition;
	}
}


double OrbitSolver::fast_time(StellarSystem &system,
		double max_step, double planet_formation_age,
		double planet_formation_semimajor, double start_age,
		EvolModeType initial_evol_mode,
		std::valarray<double> start_orbit) {
  (*this)(system, max_step, planet_formation_age, 
	  planet_formation_semimajor, start_age, initial_evol_mode,
	  start_orbit);
	Star star = system.get_star();
	std::list<double>::iterator a = tabulated_orbit[0].begin();
	std::list<double>::iterator Lc = tabulated_orbit[1].begin();
	std::list<double>::iterator Lr = tabulated_orbit[2].begin();
	double tot_time = 0;
	double prev_age = -1;
	for(std::list<double>::iterator ages = tabulated_ages.begin();
			ages != tabulated_ages.end(); ages++, a++, Lc++,  Lr++) {
		double age = *ages;
		if (age < main_seq_start) continue;
		double spin = (*Lc)/star.moment_of_inertia(age, convective);
		//std::cout << age << " " << Rsun6p5_to_AU(*a) << " " << spin << std::endl;
		if (spin > spin_thres && prev_age != -1)
			tot_time += age - prev_age;
		prev_age = age;

	}
	return tot_time;

}

///Returns the values of the variables at the tabulated ages.
const std::list<double>
*OrbitSolver::get_tabulated_var(EvolVarType var_type) const
{
	if(var_type==AGE) return &tabulated_ages;
	return &(tabulated_orbit[var_type]);
}

///Returns the derivative of the independent variable at the points
//specified by get_tabulated_indep_var.
const std::list<double> *OrbitSolver::get_tabulated_var_deriv(
		EvolVarType var_type) const
{
	assert(var_type!=AGE);
	return &(tabulated_deriv[var_type]);
}

double OrbitSolver::last_error(double age, double a, double Lc)
{
	/*returns the L2 norm of the last fractional error between the last
	simulated orbit, and the arguments (a, Lc). If the last simulated age
	is not age, returns NaN. */
	double last_age = tabulated_ages.back();
	if (abs(age - last_age) > 1e-6) return NaN;
	double last_a6p5 = (tabulated_orbit[0]).back();
	double last_a = pow(last_a6p5, 1.0/6.5)*
			AstroConst::solar_radius/AstroConst::AU;
	double last_Lc = (tabulated_orbit[1]).back();
	double sqr_err = pow((last_a - a)/a, 2) + pow((last_Lc - Lc)/2, 2);
	return std::sqrt(sqr_err);
}


