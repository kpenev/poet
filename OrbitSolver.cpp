/**\file
 *
 * \brief Implements some of the members of the OrbitSolver class, the
 * various stopping conditions and a number of other classes used while
 * calculating the orbital evolution.
 * 
 * \ingroup OrbitSolver_group
 */

#include "OrbitSolver.h"
#include <iostream>
#include <iomanip>

const double Rsun_AU=AstroConst::solar_radius/AstroConst::AU;

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

std::ostream &operator<<(std::ostream &os, const StopInformation &stop)
{
	std::streamsize orig_precision=os.precision();
	os.precision(16);
	std::ios_base::fmtflags orig_flags=os.flags();
	os.setf(std::ios_base::scientific);
	os << "Stop at t=" << stop.stop_age()
		<< (stop.is_crossing() ? ", crossing" : ", extremum")
		<< " of " << stop.stop_reason() << ", precision="
		<< stop.stop_condition_precision(); 
	os.precision(orig_precision);
	os.flags(orig_flags);
	return os;
}

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
				   (dwconv_dt>0 ? Inf : -Inf), false)/Rsun_AU,
		   da_dt=derivatives[0];
	stop_deriv.resize(1, NaN);
	if(orbit[0]==0) *static_cast<int*>(NULL)=0;
	if(derivatives[0]==0) return std::valarray<double>(Inf, 1);
	return std::valarray<double>(
			(max_semi_evol - da_dt)/da_dt, 1);
}

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

double convective_frequency(double age, const StellarSystem &system, 
		const std::valarray<double> &orbit, EvolModeType evol_mode)
{
	switch(evol_mode) {
		case FAST_PLANET: case SLOW_PLANET: case TABULATION:
			return system.get_star().spin_frequency(age, convective,
					orbit[1]);
		case LOCKED_TO_PLANET:
			return system.get_planet().orbital_angular_velocity_semimajor(
					orbit[0]*Rsun_AU);
		case NO_PLANET:
			return system.get_star().spin_frequency(age, convective,
					orbit[0]);
		case LOCKED_TO_DISK:
			throw Error::BadFunctionArguments("The convective spin frequency"
					" is not defined for disk locked evolution mode.");
		default:
			throw Error::BadFunctionArguments(
					"Unrecognized wind saturation argument");
	}
}

std::valarray<double> WindSaturationCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double wsat=system.get_star().get_wind_saturation_frequency(),
		   wconv=convective_frequency(age, system, orbit, evol_mode);

	stop_deriv.resize(1, NaN);
	if(std::isinf(wsat)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>((wconv-wsat)/wsat, 1);
}

std::valarray<double> RotFastCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	if(!std::isfinite(spin_thres)) return std::valarray<double>(-1, 1);
	double wconv=convective_frequency(age, system, orbit, evol_mode);
	stop_deriv.resize(1, NaN);
	return std::valarray<double>((wconv-spin_thres)/spin_thres, 1);
}

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const CombinedStoppingCondition &rhs)
{
	for(std::vector<const StoppingCondition *>::const_iterator
			i=rhs.sub_conditions.begin(); i!=rhs.sub_conditions.end(); i++)
		sub_conditions.push_back(*i);
	return *this;
}

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const StoppingCondition *rhs)
{
	sub_conditions.push_back(rhs);
	return *this;
}

CombinedStoppingCondition::~CombinedStoppingCondition()
{
	for(std::vector<const StoppingCondition *>::const_iterator
			i=sub_conditions.begin(); i!=sub_conditions.end(); i++)
		delete *i;
}

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

void StopHistoryInterval::advance_iterator_set(
		std::list<double>::const_iterator &age_i,
		std::list< std::valarray<double> >::const_iterator &cond_i,
		std::list< std::valarray<double> >::const_iterator &deriv_i)
{
	++age_i; ++cond_i; ++deriv_i;
	if(age_i==__history_age_end) {
		assert(cond_i==__stop_cond_history_end &&
				deriv_i==__stop_deriv_history_end);
		age_i=__discarded_age_begin;
		cond_i=__stop_cond_discarded_begin;
		deriv_i=__stop_deriv_discarded_begin;
	}
}

void StopHistoryInterval::retreat_iterator_set(
		std::list<double>::const_iterator &age_i,
		std::list< std::valarray<double> >::const_iterator &cond_i,
		std::list< std::valarray<double> >::const_iterator &deriv_i)
{
	if(age_i==__discarded_age_begin) {
		assert(cond_i==__stop_cond_discarded_begin &&
				deriv_i==__stop_deriv_discarded_begin);
		age_i=__history_age_end;
		cond_i=__stop_cond_history_end;
		deriv_i=__stop_deriv_history_end;
	}
	--age_i; --cond_i; --deriv_i;
}

StopHistoryInterval::StopHistoryInterval(size_t num_points,
		std::list<double>::const_iterator first_age,
		std::list<double>::const_iterator history_age_end,
		std::list<double>::const_iterator discarded_age_begin,		
		std::list< std::valarray<double> >::const_iterator first_stop_cond,
		std::list< std::valarray<double> >::const_iterator
			stop_cond_history_end,
		std::list< std::valarray<double> >::const_iterator
			stop_cond_discarded_begin,
		std::list< std::valarray<double> >::const_iterator first_stop_deriv,
		std::list< std::valarray<double> >::const_iterator
			stop_deriv_history_end,
		std::list< std::valarray<double> >::const_iterator
		stop_deriv_discarded_begin) :
			__num_points(num_points), __point_i(0), __first_age(first_age),
			__last_age(first_age), __history_age_end(history_age_end),
			__discarded_age_begin(discarded_age_begin), __age_i(first_age),
			__first_stop_cond(first_stop_cond),
			__last_stop_cond(first_stop_cond),
			__stop_cond_history_end(stop_cond_history_end),
			__stop_cond_discarded_begin(stop_cond_discarded_begin),
			__stop_cond_i(first_stop_cond),
			__first_stop_deriv(first_stop_deriv),
			__last_stop_deriv(first_stop_deriv),
			__stop_deriv_history_end(stop_deriv_history_end),
			__stop_deriv_discarded_begin(stop_deriv_discarded_begin),
			__stop_deriv_i(first_stop_deriv)
{
	if(num_points) advance_iterator_set(__last_age, __last_stop_cond,
			__last_stop_deriv);
	reset();
}

StopHistoryInterval::StopHistoryInterval(const StopHistoryInterval &orig) :
	__num_points(orig.__num_points), __point_i(orig.__point_i),
	__first_age(orig.__first_age),
	__last_age(orig.__last_age),
	__history_age_end(orig.__history_age_end),
	__discarded_age_begin(orig.__discarded_age_begin), __age_i(orig.__age_i),
	__first_stop_cond(orig.__first_stop_cond),
	__last_stop_cond(orig.__last_stop_cond),
	__stop_cond_history_end(orig.__stop_cond_history_end),
	__stop_cond_discarded_begin(orig.__stop_cond_discarded_begin),
	__stop_cond_i(orig.__stop_cond_i),
	__first_stop_deriv(orig.__first_stop_deriv),
	__last_stop_deriv(orig.__last_stop_deriv),
	__stop_deriv_history_end(orig.__stop_deriv_history_end),
	__stop_deriv_discarded_begin(orig.__stop_deriv_discarded_begin),
	__stop_deriv_i(orig.__stop_deriv_i)
{}

void StopHistoryInterval::reset()
{
	__point_i=0;
	__age_i=__first_age;
	__stop_cond_i=__first_stop_cond;
	__stop_deriv_i=__first_stop_deriv;
}

StopHistoryInterval &StopHistoryInterval::operator++()
{
	++__point_i;
	if(__point_i>__num_points)
		throw Error::Runtime("Attempting to increment two points past the "
				"end of a StopHistoryInterval!");
	advance_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
	return *this;
}

StopHistoryInterval StopHistoryInterval::operator++(int)
{
	StopHistoryInterval result(*this);
	operator++();
	return result;
}

StopHistoryInterval &StopHistoryInterval::operator--()
{
	if(__point_i==0) throw Error::Runtime("Attempting to go before the "
			"beginning of a StopHistoryInterval.");
	--__point_i;
	retreat_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
	return *this;
}

StopHistoryInterval StopHistoryInterval::operator--(int)
{
	StopHistoryInterval result(*this);
	operator--();
	return result;
}

StopHistoryInterval &StopHistoryInterval::operator<<(size_t n)
{
	for(size_t i=0; i<n; i++) {
		retreat_iterator_set(__first_age, __first_stop_cond,
				__first_stop_deriv);
		retreat_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
		retreat_iterator_set(__last_age, __last_stop_cond,__last_stop_deriv);
	}
	return *this;
}

StopHistoryInterval &StopHistoryInterval::operator>>(size_t n)
{
	for(size_t i=0; i<n; i++) {
		advance_iterator_set(__last_age, __last_stop_cond, __last_stop_deriv);
		advance_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
		advance_iterator_set(__first_age, __first_stop_cond, __first_stop_deriv);
	}
	return *this;
}

StopHistoryInterval &StopHistoryInterval::operator=(
		const StopHistoryInterval &rhs)
{
	__num_points=rhs.__num_points;
	__point_i=rhs.__point_i;
	__first_age=rhs.__first_age;
	__last_age=rhs.__last_age;
	__history_age_end=rhs.__history_age_end;
	__discarded_age_begin=rhs.__discarded_age_begin;
	__age_i=rhs.__age_i;
	__first_stop_cond=rhs.__first_stop_cond;
	__last_stop_cond=rhs.__last_stop_cond;
	__stop_cond_history_end=rhs.__stop_cond_history_end;
	__stop_cond_discarded_begin=rhs.__stop_cond_discarded_begin;
	__stop_cond_i=rhs.__stop_cond_i;
	__first_stop_deriv=rhs.__first_stop_deriv;
	__last_stop_deriv=rhs.__last_stop_deriv;
	__stop_deriv_history_end=rhs.__stop_deriv_history_end;
	__stop_deriv_discarded_begin=rhs.__stop_deriv_discarded_begin;
	__stop_deriv_i=rhs.__stop_deriv_i;
	return *this;
}

bool StopHistoryInterval::operator==(const StopHistoryInterval &rhs)
{
	return __num_points==rhs.__num_points &&
		__point_i==rhs.__point_i &&
		__first_age==rhs.__first_age &&
		__last_age==rhs.__last_age &&
		__history_age_end==rhs.__history_age_end &&
		__discarded_age_begin==rhs.__discarded_age_begin &&
		__age_i==rhs.__age_i &&
		__first_stop_cond==rhs.__first_stop_cond &&
		__last_stop_cond==rhs.__last_stop_cond &&
		__stop_cond_history_end==rhs.__stop_cond_history_end &&
		__stop_cond_discarded_begin==rhs.__stop_cond_discarded_begin &&
		__stop_cond_i==rhs.__stop_cond_i &&
		__first_stop_deriv==rhs.__first_stop_deriv &&
		__last_stop_deriv==rhs.__last_stop_deriv &&
		__stop_deriv_history_end==rhs.__stop_deriv_history_end &&
		__stop_deriv_discarded_begin==rhs.__stop_deriv_discarded_begin &&
		__stop_deriv_i==rhs.__stop_deriv_i;
}

void StopHistoryInterval::grow_left(size_t n)
{
	__num_points+=n;
	for(size_t i=0; i<n; i++) {
		retreat_iterator_set(__first_age, __first_stop_cond,
				__first_stop_deriv);
		retreat_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
	}
}

void StopHistoryInterval::grow_right(size_t n)
{
	__num_points+=n;
	for(size_t i=0; i<n; i++)
		advance_iterator_set(__last_age, __last_stop_cond,__last_stop_deriv);
}

std::ostream &operator<<(std::ostream &os, StopHistoryInterval interval)
{
	std::streamsize orig_precision=os.precision();
	std::ios_base::fmtflags orig_flags=os.flags();
	os.setf(std::ios_base::scientific);
	os.precision(16);
	os << std::endl;
	size_t current_index=interval.current_point_index();
	os << std::setw(20) << "Age:";
	for(interval.reset(); !interval.end(); interval++) {
		if(current_index==interval.current_point_index())
			os << "|";
		os << std::setw(25) << interval.age();
		if(current_index==interval.current_point_index())
			os << "|";
	}
	os << std::endl;
	os << "number conditions: " << interval.number_conditions() << std::endl;
	for(size_t cond_ind=0; cond_ind<interval.number_conditions();
			cond_ind++) {
		os << std::setw(13) << "Condition[" << std::setw(5)
			<< cond_ind << "]:";
			for(interval.reset(); !interval.end(); interval++) {
				if(current_index==interval.current_point_index())
					os << "|";
				os << std::setw(25)
					<< interval.stop_condition_value(cond_ind);
				if(current_index==interval.current_point_index())
					os << "|";
			}
		os << std::endl;
		os << std::setw(13) << "Derivative[" << std::setw(5)
			<< cond_ind << "]:";
		for(interval.reset(); !interval.end(); interval++) {
			if(current_index==interval.current_point_index())
				os << "|";
			os << std::setw(25) << interval.stop_condition_deriv(cond_ind);
			if(current_index==interval.current_point_index())
				os << "|";
		}
		os << std::endl;
	}
	os.precision(orig_precision);
	os.flags(orig_flags);
	return os;
}

StopHistoryInterval OrbitSolver::select_stop_condition_interval(
		bool crossing, size_t cond_ind, size_t max_points) const
{
	if(crossing) assert(stop_cond_history.size()-
			skip_history_zerocrossing[cond_ind]>=1);
	else assert(stop_cond_history.size()-
			skip_history_zerocrossing[cond_ind]>=2);
	size_t num_points=std::min(max_points, stop_history_ages.size()+
			discarded_stop_ages.size()-
			(crossing ? skip_history_zerocrossing[cond_ind] : 0));
	std::list<double>::const_iterator first_age=stop_history_ages.end();
	std::list< std::valarray<double> >::const_iterator
		first_stop_cond=stop_cond_history.end(),
		first_stop_deriv=stop_deriv_history.end();
	int go_back=static_cast<int>(num_points)-
			static_cast<int>(discarded_stop_ages.size());
	go_back=std::max(go_back, (crossing ? 1 : 2));
	size_t failed_back=0;
	for(int i=0; i<go_back; i++) {
		--first_age;
		--first_stop_cond;
		--first_stop_deriv;
		if(!crossing && *first_age < skip_history_extremum[cond_ind]) {
			++failed_back;
			++first_age;
			++first_stop_cond;
			++first_stop_deriv;
			--num_points;
		}
	}
	if(go_back-failed_back<(crossing ? 1 : 2)) return StopHistoryInterval();
	StopHistoryInterval interval(num_points, first_age,
			stop_history_ages.end(), discarded_stop_ages.begin(),
			first_stop_cond, stop_cond_history.end(),
			stop_cond_discarded.begin(), first_stop_deriv,
			stop_deriv_history.end(), stop_deriv_discarded.begin()),
						result=interval;
	int max_left_shift=std::min(discarded_stop_ages.size()-1,
			num_points-(crossing ? 2 : 3));
	int history_limit=0;
	if(crossing)
		history_limit=stop_history_ages.size()-go_back+failed_back-
			   skip_history_zerocrossing[cond_ind];
	else while(first_age!=stop_history_ages.begin() &&
			*(--first_age)>skip_history_extremum[cond_ind])
		++history_limit;
	max_left_shift=std::min(history_limit, max_left_shift);
	for(int i=0; i<max_left_shift; i++) {
		interval << 1;
		if((interval.last_age()-interval.first_age())<
				result.last_age()-result.first_age()) result=interval;
	}
	return result;
}

void OrbitSolver::output_history_and_discarded(std::ostream &os)
{
	std::streamsize orig_precision=os.precision();
	os.precision(16);
	std::ios_base::fmtflags orig_flags=os.flags();
	os.setf(std::ios_base::scientific);

	os << "Stored stop condition information:" << std::endl
		<< std::setw(20) << "Age:";
	for(std::list<double>::const_iterator age_i=stop_history_ages.begin();
			age_i!=discarded_stop_ages.end(); age_i++) {
		if(age_i==stop_history_ages.end()) {
			os << "|";
			age_i=discarded_stop_ages.begin();
		}
		os << std::setw(28) << *age_i;
	}
	std::string hline;
	hline.assign(20+25*(stop_history_ages.size()+discarded_stop_ages.size()),
			'_');
	os << std::endl << hline << std::endl;

	for(size_t i=0; i<stop_cond_discarded.front().size(); i++) {
		os << std::setw(13) << "Condition[" << std::setw(5) << i
			<< "]:";
		size_t point_num=0;
		std::list<double>::const_iterator age_i=stop_history_ages.begin();
		bool marked_skip_extremum=false;
		for(std::list< std::valarray<double> >::const_iterator
				cond_i=stop_cond_history.begin();
				cond_i!=stop_cond_discarded.end(); cond_i++) {
			if(cond_i==stop_cond_history.end()) {
				os << "|";
				cond_i=stop_cond_discarded.begin();
			}
			bool marked=false;
			if(point_num==skip_history_zerocrossing[i]) {
				os << "z"; marked=true;
			} else os << " ";
			if((*age_i)>skip_history_extremum[i] && !marked_skip_extremum) {
				os << "e"; marked_skip_extremum=true; marked=true;
			} else os << (marked ? "-" : " ");
			os << (marked ? ">" : " ");
			os << std::setw(25) << (*cond_i)[i];
			point_num++;
			age_i++;
		}
		os << std::endl;
		os << std::setw(13) << "Derivative[" << std::setw(5) << i
			<< "]:";
		for(std::list< std::valarray<double> >::const_iterator
				deriv_i=stop_deriv_history.begin();
				deriv_i!=stop_deriv_discarded.end(); deriv_i++) {
			if(deriv_i==stop_deriv_history.end()) {
				os << "|";
				deriv_i=stop_deriv_discarded.begin();
			}
			os << std::setw(28) << (*deriv_i)[i];
		}
		os << std::endl;
	}

	os.precision(orig_precision);
	os.flags(orig_flags);
}

void OrbitSolver::clear_discarded()
{
	discarded_stop_ages.clear();
	stop_cond_discarded.clear();
	stop_deriv_discarded.clear();
}

void OrbitSolver::insert_discarded(double age,
		const std::valarray<double> &current_stop_cond,
		const std::valarray<double> &current_stop_deriv)
{
	std::list<double>::iterator age_i=discarded_stop_ages.begin();
	std::list< std::valarray<double> >::iterator
		cond_i=stop_cond_discarded.begin(),
		deriv_i=stop_deriv_discarded.begin();
	while(age_i!=discarded_stop_ages.end() && age>(*age_i)) {
		age_i++; cond_i++; deriv_i++;
	}
	discarded_stop_ages.insert(age_i, age);
	stop_cond_discarded.insert(cond_i, current_stop_cond);
	stop_deriv_discarded.insert(deriv_i, current_stop_deriv);
}

void OrbitSolver::append_to_orbit(const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		EvolModeType evolution_mode, double age, const StellarSystem &system,
		double planet_formation_semimajor)
{
	clear_discarded();
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

void OrbitSolver::clear_history()
{
	orbit_history.clear();
	orbit_deriv_history.clear();
	stop_history_ages.clear();
	stop_cond_history.clear();
	stop_deriv_history.clear();
	clear_discarded();
}

double OrbitSolver::go_back(double max_age, std::valarray<double> &orbit,
		std::valarray<double> &derivatives)
{
	while(max_age<tabulated_ages.back()) {
		tabulated_ages.pop_back();
		tabulated_evolution_mode.pop_back();
		for(size_t i=0; i<tabulated_orbit.size(); i++) {
			tabulated_orbit[i].pop_back();
			tabulated_deriv[i].pop_back();
		}
	}
	if(max_age<stop_history_ages.back()) clear_discarded();
	while(max_age<stop_history_ages.back()) {
		stop_history_ages.pop_back();
		orbit_history.pop_back();
		orbit_deriv_history.pop_back();
		stop_cond_history.pop_back();
		stop_deriv_history.pop_back();
	}
	for(size_t i=0; i<skip_history_zerocrossing.size(); i++) {
		skip_history_zerocrossing[i]=
			std::min(skip_history_zerocrossing[i], stop_history_ages.size());
		skip_history_extremum[i]=
			std::min(stop_history_ages.back(), skip_history_extremum[i]);
	}
	orbit=orbit_history.back();
	derivatives=orbit_deriv_history.back();
	return stop_history_ages.back();
}

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

std::ostream &operator<<(std::ostream &os, const std::list<double> &l)
{
	for(std::list<double>::const_iterator i=l.begin();
			i!=l.end(); i++) os << *i << ", ";
	return os;
}

ExtremumInformation OrbitSolver::extremum_from_history_no_deriv(
		size_t condition_index) const
{
	ExtremumInformation result;
	if(stop_history_ages.size()-
			skip_history_zerocrossing[condition_index]<2) return result;
	std::list< std::valarray<double> >::const_iterator
		stop_cond_i=stop_cond_history.end();
	double pre1_cond=(*(--stop_cond_i))[condition_index],
		   pre2_cond=(*(--stop_cond_i))[condition_index],
		   post_cond=stop_cond_discarded.front()[condition_index];
	if(std::abs(pre1_cond)>std::abs(pre2_cond) ||
			std::abs(pre1_cond)>std::abs(post_cond) ||
			(!std::isfinite(pre1_cond) && !std::isfinite(pre2_cond) &&
			 !std::isfinite(post_cond)))
		return result;

	StopHistoryInterval stop_interval=select_stop_condition_interval(
			false, condition_index, 4);
	if(stop_interval.num_points()<3) return result;
	double t0=stop_interval.age(),
		   c0=stop_interval.stop_condition_value(condition_index),
		   t1=(++stop_interval).age(),
		   c1=stop_interval.stop_condition_value(condition_index),
		   t2=(++stop_interval).age(),
		   c2=stop_interval.stop_condition_value(condition_index);
	if(stop_interval.num_points()==3) {
		if((c1-c0)*(c2-c1)>0) return result;
		result.x()=quadratic_extremum(t0, c0, t1, c1, t2, c2,
				&(result.y()));
	} else {
		double t3=(++stop_interval).age(),
			   c3=stop_interval.stop_condition_value(condition_index);
		if((c1-c0)*(c2-c1)>0 && (c2-c1)*(c3-c2)>0) return result;
		double range_low, range_high;
		if(std::abs(c1)<=std::abs(c0) && std::abs(c1)<=std::abs(c2)) {
			range_low=t0; range_high=t2;
		} else if(std::abs(c2)<=std::abs(c1) && std::abs(c2)<=std::abs(c3)) {
			range_low=t1; range_high=t3;
		}
		result.x()=cubic_extremum(t0, c0, t1, c1, t2, c2, t3, c3,
				&(result.y()), range_low, range_high);
	}
	return result;
}

ExtremumInformation OrbitSolver::extremum_from_history(
		size_t condition_index) const
{
	ExtremumInformation result;
	if(stop_history_ages.size()==0 ||
			stop_history_ages.back()<skip_history_extremum[condition_index])
		return result;
	double prev_stop_cond=stop_cond_history.back()[condition_index],
		   next_stop_cond=stop_cond_discarded.front()[condition_index];
	if(next_stop_cond*prev_stop_cond<=0)
		return result;
	if(std::isnan(stop_deriv_history.back()[condition_index]))
		return extremum_from_history_no_deriv(condition_index);
	double prev_stop_deriv=stop_deriv_history.back()[condition_index],
		   next_stop_deriv=stop_deriv_discarded.front()[condition_index];
	if(next_stop_cond*next_stop_deriv<0 ||next_stop_deriv*prev_stop_deriv>=0)
		return result;
	result.x()=estimate_extremum(stop_history_ages.back(), prev_stop_cond,
			discarded_stop_ages.front(), next_stop_cond, prev_stop_deriv,
			next_stop_deriv, &(result.y()));
	return result;
}

double OrbitSolver::crossing_from_history_no_deriv(size_t condition_index)
	const
{
	StopHistoryInterval stop_interval=select_stop_condition_interval(true,
			condition_index, 4);
	if(stop_interval.num_points()<2) return Inf;
	double t0=stop_interval.age(),
		   c0=stop_interval.stop_condition_value(condition_index),
		   t1=(++stop_interval).age(),
		   c1=stop_interval.stop_condition_value(condition_index);
	if(stop_interval.num_points()==2)
		return estimate_zerocrossing(t0, c0, t1, c1);
	double t2=(++stop_interval).age(),
		   c2=stop_interval.stop_condition_value(condition_index);
	double range_low, range_high;
	if(c0*c1<=0) {range_low=t0; range_high=t1;}
	else if(c1*c2<=0) {range_low=t1; range_high=t2;}
	if(stop_interval.num_points()==3)
		return quadratic_zerocrossing(t0, c0, t1, c1, t2, c2,
				range_low, range_high);
	double t3=(++stop_interval).age(),
		   c3=stop_interval.stop_condition_value(condition_index);
	if(c0*c1>0 && c1*c2>0) {range_low=t2; range_high=t3;}
	return cubic_zerocrossing(t0, c0, t1, c1, t2, c2, t3, c3,
			range_low, range_high);
}

double OrbitSolver::crossing_from_history(size_t condition_index) const
{
	if(stop_history_ages.size()==0 || 
			stop_history_ages.size()==
			skip_history_zerocrossing[condition_index]) return Inf;
	double prev_stop_cond=stop_cond_history.back()[condition_index];
	double next_stop_cond=stop_cond_discarded.front()[condition_index];
	if(next_stop_cond*prev_stop_cond>0) return Inf;
	if(std::isnan(stop_deriv_history.back()[condition_index]))
		return crossing_from_history_no_deriv(condition_index);
	double prev_stop_deriv=stop_deriv_history.back()[condition_index],
		   prev_age=stop_history_ages.back(),
		   next_stop_deriv=stop_deriv_discarded.front()[condition_index],
		   next_age=discarded_stop_ages.front();
	return estimate_zerocrossing(prev_age, prev_stop_cond, next_age,
			next_stop_cond, prev_stop_deriv, next_stop_deriv);
}

void OrbitSolver::initialize_skip_history(const StoppingCondition &stop_cond,
			StoppingConditionType stop_reason)
{
	skip_history_zerocrossing.resize(stop_cond.num_subconditions(), 0);
	skip_history_extremum.resize(stop_cond.num_subconditions(), 0);
	for(size_t cond_ind=0; cond_ind<stop_cond.num_subconditions();
			cond_ind++) {
		StoppingConditionType stop_cond_type=stop_cond.type(cond_ind);
		if((stop_reason==BREAK_LOCK && stop_cond_type==SYNCHRONIZED) ||
				(stop_reason!=NO_STOP && stop_cond_type==stop_reason)
				) {
			skip_history_zerocrossing[cond_ind]=1;
			skip_history_extremum[cond_ind]=stop_history_ages.front()*
				(1.0+std::numeric_limits<double>::epsilon());
		}
	}
}

void OrbitSolver::update_skip_history(
		const std::valarray<double> &current_stop_cond,
		const StoppingCondition &stop_cond,
		const StopInformation &stop_info)
{
	size_t history_size=stop_history_ages.size();
	for(size_t i=0; i<current_stop_cond.size(); i++) {
		if(skip_history_zerocrossing[i]==history_size &&
				std::abs(current_stop_cond[i])<precision)
			++skip_history_zerocrossing[i];
		if(stop_info.stop_reason()==stop_cond.type(i) &&
				!stop_info.is_crossing() &&
				std::abs(stop_info.stop_condition_precision())<=precision)
			skip_history_extremum[i]=stop_info.stop_age();
	}
}

bool OrbitSolver::acceptable_step(double age,
		const StopInformation &stop_info)
{
	return stop_info.stop_age()>=age ||
			std::abs(stop_info.stop_condition_precision())<=precision;
}

StopInformation OrbitSolver::update_stop_condition_history(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system, EvolModeType evolution_mode,
		const StoppingCondition &stop_cond,
		StoppingConditionType stop_reason)
{
	std::valarray<double> current_stop_cond(stop_cond.num_subconditions()),
		current_stop_deriv;
	current_stop_cond=stop_cond(age, orbit, derivatives, system,
			current_stop_deriv, evolution_mode);

	if(stop_history_ages.size()==0)
		initialize_skip_history(stop_cond, stop_reason);
	StopInformation result;
	insert_discarded(age, current_stop_cond, current_stop_deriv);
	for(size_t cond_ind=0; cond_ind<current_stop_cond.size(); cond_ind++) {
		StoppingConditionType
			stop_cond_type=stop_cond.type(cond_ind);
		double stop_cond_value=current_stop_cond[cond_ind],
			   crossing_age=crossing_from_history(cond_ind);
		ExtremumInformation extremum=extremum_from_history(cond_ind);
		double extremum_precision=std::min(
				std::abs(extremum.y()-stop_cond_value),
				std::abs(extremum.y()-
					stop_cond_history.back()[cond_ind]))/
			std::abs(extremum.y());
		bool is_crossing=crossing_age<=extremum.x();
		StopInformation stop_info(std::min(crossing_age, extremum.x()),
				(is_crossing ? stop_cond_value : extremum_precision),
				stop_cond_type, is_crossing, cond_ind);

		if((!acceptable_step(age, stop_info) || is_crossing) &&
				stop_info.stop_age()<result.stop_age()) result=stop_info;
	}
	if(acceptable_step(age, result)) {
		update_skip_history(current_stop_cond, stop_cond, result);
		stop_history_ages.push_back(age);
		stop_cond_history.push_back(current_stop_cond);
		stop_deriv_history.push_back(current_stop_deriv);
		orbit_history.push_back(orbit);
		orbit_deriv_history.push_back(derivatives);
	}
	return result;
}

bool OrbitSolver::evolve_until(StellarSystem *system, double start_age,
		double &max_age, std::valarray<double> &orbit,
		double &stop_condition_value, StoppingConditionType &stop_reason,
		double max_step, EvolModeType evolution_mode,
		WindSaturationState wind_state, const StoppingCondition &stop_cond,
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
			evolution_mode, stop_cond, stop_reason);
	clear_discarded();
	double step_size=0.01*(max_age-start_age);

	stop_reason=NO_STOP;
	bool result;
	while(t<max_age) {
		double max_next_t=std::min(t + max_step, max_age); 
		int status=GSL_SUCCESS;
		StopInformation stop;
		bool step_rejected=false;
		do {
			status=gsl_odeiv2_evolve_apply(
					evolve, step_control, step, &ode_system,
					&t, max_next_t, &step_size, &(orbit[0]));
			stellar_system_diff_eq(t, &(orbit[0]), &(derivatives[0]),
					sys_mode_windstate);
			stop=update_stop_condition_history(t, orbit, derivatives, *system,
					evolution_mode, stop_cond, stop_reason);
			if(!acceptable_step(t, stop)) {
				t=go_back(stop.stop_age(), orbit, derivatives);
				if(stop.is_crossing())
					stop.stop_condition_precision()=
						stop_cond_history.back()[
						stop.stop_condition_index()];
				max_next_t=stop.stop_age();
				step_rejected=true;
			} else step_rejected=false;
			if(stop.is_crossing())
				stop_condition_value=stop.stop_condition_precision();
		} while(step_rejected &&
				std::abs(stop.stop_condition_precision())>precision);
		if (status != GSL_SUCCESS) {
			std::ostringstream msg;
			msg << "GSL signaled failure while evolving (error code " <<
				status << ")";
			throw Error::Runtime(msg.str());
		}
		if(!step_rejected)
			append_to_orbit(orbit, derivatives,
				evolution_mode, t, *system, planet_formation_semimajor);
		if(stop.is_crossing() && stop.stop_reason()!=NO_STOP) {
			stop_reason=stop.stop_reason();
			result=(t>stop.stop_age());
			break;
		}
	}
	max_age=t;
	clear_history();

	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(step_control);
	gsl_odeiv2_step_free(step);
	return result;
}

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

EvolModeType OrbitSolver::next_evol_mode(double age,
		const std::valarray<double> &orbit,
		double initial_semimajor, const StellarSystem &system,
		EvolModeType evolution_mode,
		StoppingConditionType condition_type,
		double condition_value, bool stopped_before,
		double planet_formation_age) const
{
	if(stopped_before) condition_value*=-1;
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
	if (max_age > MAX_END_AGE) end_age = MAX_END_AGE;
}

void OrbitSolver::operator()(StellarSystem &system, double max_step,
		double planet_formation_age, double planet_formation_semimajor,
		double start_age, EvolModeType initial_evol_mode,
		const std::valarray<double> &start_orbit, bool no_evol_mode_change)
{
	if(std::isnan(start_age))
		start_age=system.get_star().core_formation_age();
	double wconv;
	if(initial_evol_mode==FAST_PLANET || initial_evol_mode==SLOW_PLANET
			|| initial_evol_mode==LOCKED_TO_PLANET) {
		if(start_orbit[0]<=
				system.get_planet().minimum_semimajor(start_age)/Rsun_AU)
			throw Error::BadFunctionArguments("Attempting to calculate the "
					"evolution of a system where the planet starts "
					"inside the Roche radius or its parent star!");
	} else if(std::isfinite(planet_formation_age)) {
		if(planet_formation_semimajor/Rsun_AU<=
				system.get_planet().minimum_semimajor(
					std::max(planet_formation_age,
						system.get_star().get_disk_dissipation_age()))/
				Rsun_AU)
			throw Error::BadFunctionArguments("Attempting to calculate the "
					"evolution of a system where the planet will be formed "
					"inside the Roche radius or its parent star!");
	}
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
	StoppingConditionType stop_reason=NO_STOP;
	double last_age=start_age;
	std::valarray<double> orbit=start_orbit;
	EvolModeType evolution_mode=initial_evol_mode;
	while(last_age<end_age) {
/*		if (evolution_mode == FAST_PLANET || evolution_mode == SLOW_PLANET ||
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
		//std::cout << deathResult << std::endl;*/
		double stop_evol_age=stopping_age(last_age, evolution_mode,
				system, planet_formation_age);
		CombinedStoppingCondition *stopping_condition=get_stopping_condition(
				evolution_mode, planet_formation_semimajor,
				&system.get_planet());
		double stop_condition_value;
		bool stopped_before=!evolve_until(&system, start_age, stop_evol_age,
				orbit, stop_condition_value, stop_reason, max_step,
				evolution_mode, wind_state, *stopping_condition,
				planet_formation_semimajor);
		last_age=stop_evol_age;
		if(last_age<end_age) {
			EvolModeType old_evolution_mode=evolution_mode;
			if(!no_evol_mode_change) evolution_mode=next_evol_mode(last_age,
					orbit, planet_formation_semimajor, system,
					evolution_mode, stop_reason, stop_condition_value,
					stopped_before, planet_formation_age);
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


const std::list<double>
*OrbitSolver::get_tabulated_var(EvolVarType var_type) const
{
	if(var_type==AGE) return &tabulated_ages;
	return &(tabulated_orbit[var_type]);
}

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


