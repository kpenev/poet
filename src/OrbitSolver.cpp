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

const double OrbitSolver::MAX_END_AGE = 10;

std::ostream &operator<<(std::ostream &os, const EvolVarType &evol_var)
{
	switch(evol_var) {
		case AGE : os << "AGE"; break;
		case SEMIMAJOR : os << "SEMIMAJOR"; break;
		case INCLINATION : os << "INCLINATION"; break;
		case LCONV : os << "LCONV"; break;
		case LRAD_PAR : os << "LRAD_PAR";
		case LRAD_PERP : os << "LRAD_PERP";
		default :
#ifdef DEBUG
						 assert(false)
#endif
							 ;
	};
	return os;
}

std::ostream &operator<<(std::ostream &os, const StopInformation &stop)
{
	return os;
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

void get_tidal_dissipation(StellarSystem &system, 
		double age, const double *parameters, SpinOrbitLockInfo &star_lock,
		TidalDissipation &dissipation)
{
	double semimajor=(star_lock ? parameters[0] :
			std::pow(parameters[0], 1.0/6.5))*Rsun_AU;
	SpinOrbitLockInfo planet_lock(1, 1, 0);
	dissipation.init(
			system.get_star(), system.get_planet(), semimajor, 0, star_lock,
			planet_lock);
}

int stellar_system_diff_eq(double age, const double *parameters,
		double *derivatives, void *system_mode_windstate_lock_dissipation)
{
	void **input_params=static_cast<void **>(
			system_mode_windstate_lock_dissipation);
	StellarSystem &system=*static_cast< StellarSystem* >(
			input_params[0]);
	EvolModeType evol_mode=*static_cast< EvolModeType* >(input_params[1]);
	WindSaturationState wind_state=*static_cast<WindSaturationState*>(
			input_params[2]);
	SpinOrbitLockInfo &star_lock=*static_cast<SpinOrbitLockInfo*>(
			input_params[3]);
	TidalDissipation &dissipation=*static_cast<TidalDissipation*>(
			input_params[4]);
	const Star &star=system.get_star();
	const Planet &planet=system.get_planet();
	double Lconv, inclination; 
	if(evol_mode==BINARY) {
		inclination=parameters[1];
		Lconv=(star_lock 
				? orbital_angular_velocity(star.mass(), planet.mass(),
					parameters[0]*Rsun_AU)*
				  star.moment_of_inertia(age, convective)
				: parameters[2]);
		if(std::isnan(inclination) || std::isnan(Lconv) || 
				(!star_lock && parameters[0]<0)) {
			for(int i=0; i<(star_lock ? 4 : 5); ++i) derivatives[i]=NaN;
			return GSL_EDOM;
		}
	} else {
		inclination=NaN;
		Lconv=(evol_mode==LOCKED_TO_DISK
				? star.disk_lock_frequency()*
				  star.moment_of_inertia(age, convective)
				: parameters[0]);
	}
	system(age, Lconv, inclination);
	if(evol_mode==BINARY) get_tidal_dissipation(system, age, &parameters[0],
											    star_lock, dissipation);
	int result=system.differential_equations(parameters, evol_mode,
			wind_state, star_lock, dissipation, derivatives);
	return result;
}

int stellar_system_jacobian(double age, const double *orbital_parameters,
		double *param_derivs, double *age_derivs,void *system_mode_windstate)
{
	throw Error::NotImplemented("Jacobian");
/*	void **input_params=static_cast<void **>(system_mode_windstate);
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
	return result;*/
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
	if(num_points==0) throw Error::BadFunctionArguments(
			"Attempt to contsruct a StopHistoryInterval of size 0.");
	for(size_t i=0; i<num_points-1; i++)
		advance_iterator_set(__last_age, __last_stop_cond,
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
	return os;
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

#ifdef DEBUG
void OrbitSolver::output_history_and_discarded(std::ostream &os)
{
	return;
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
#endif

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
		EvolModeType evolution_mode, const SpinOrbitLockInfo &star_lock,
		WindSaturationState wind_state, double age,
		const StellarSystem &system, double planet_formation_semimajor,
		double planet_formation_inclination)
{
	clear_discarded();
	std::valarray<double> orbit_to_tabulate=transform_orbit(evolution_mode,
			TABULATION, star_lock, star_lock, age, orbit,
			planet_formation_semimajor, planet_formation_inclination,
			system),
		deriv_to_tabulate=transform_derivatives(evolution_mode, TABULATION,
				star_lock, star_lock, age, orbit, derivatives,
				planet_formation_semimajor, system);
	size_t nargs=orbit_to_tabulate.size();
	static double last_age=0;
	bool print=age>last_age+1e-3;
	if(print) {
		std::cerr << std::setw(25) << age;
		last_age=age;
	}
	__tabulated_ages.push_back(age);
	for(size_t i=0; i<nargs; i++) {
		if(print) std::cerr << std::setw(25) << orbit_to_tabulate[i];
		__tabulated_orbit[i].push_back(orbit_to_tabulate[i]);
		__tabulated_deriv[i].push_back(deriv_to_tabulate[i]);
	}
	if(print) std::cerr << std::setw(25) << evolution_mode
		<< std::setw(25) << wind_state << std::endl;
	__tabulated_evolution_mode.push_back(evolution_mode);
    __tabulated_wind_saturation.push_back(wind_state);
	__tabulated_lock.push_back(star_lock);
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
	while(max_age<__tabulated_ages.back()) {
		__tabulated_ages.pop_back();
		__tabulated_evolution_mode.pop_back();
        __tabulated_wind_saturation.pop_back();
		__tabulated_lock.pop_back();
		for(size_t i=0; i<__tabulated_orbit.size(); i++) {
			__tabulated_orbit[i].pop_back();
			__tabulated_deriv[i].pop_back();
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

size_t OrbitSolver::ode_dimension(EvolModeType evolution_mode,
		const SpinOrbitLockInfo &star_lock)
{
	switch(evolution_mode) {
		case BINARY : return (star_lock ? 4 : 5);
		case NO_PLANET : return 3;
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
		} else throw Error::BadFunctionArguments("Searching for extremum "
				"among monotonic stopping condition values in "
				"OrbitSolver::extremum_from_history_no_deriv.");
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
	double range_low=NaN, range_high=NaN;
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
	current_stop_cond=stop_cond(orbit, derivatives, __dissipation, system,
			current_stop_deriv, evolution_mode);

	if(stop_history_ages.size()==0)
		initialize_skip_history(stop_cond, stop_reason);
	StopInformation result;
	insert_discarded(age, current_stop_cond, current_stop_deriv);
#ifdef DEBUG
//	std::cerr << std::string(77, '@') << std::endl;
	output_history_and_discarded(std::cerr);
#endif
	for(size_t cond_ind=0; cond_ind<current_stop_cond.size(); cond_ind++) {
		StoppingConditionType
			stop_cond_type=stop_cond.type(cond_ind);
		double stop_cond_value=current_stop_cond[cond_ind],
			   crossing_age=crossing_from_history(cond_ind);
		ExtremumInformation extremum=extremum_from_history(cond_ind);
		double extremum_precision;
		if(std::isnan(extremum.y())) extremum_precision=NaN;
		else extremum_precision=std::min(
				std::abs(extremum.y()-stop_cond_value),
				std::abs(extremum.y()-
					stop_cond_history.back()[cond_ind]))/
			std::abs(extremum.y());
		bool is_crossing=crossing_age<=extremum.x();
		StopInformation stop_info(std::min(crossing_age, extremum.x()),
				(is_crossing ? stop_cond_value : extremum_precision),
				stop_cond_type, is_crossing, true, cond_ind);

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

StopInformation OrbitSolver::evolve_until(StellarSystem &system,
		double start_age, double &max_age, std::valarray<double> &orbit,
		StoppingConditionType &stop_reason, double max_step,
		EvolModeType evolution_mode, const SpinOrbitLockInfo &star_lock,
		WindSaturationState wind_state, const StoppingCondition &stop_cond,
		double planet_formation_semimajor,
		double planet_formation_inclination)
{
	size_t nargs=orbit.size();
#ifdef DEBUG
	std::cerr << "Starting evolution leg in "
		<< (star_lock ? "locked " : "not locked ") << evolution_mode
		<< " mode with " << wind_state << " wind from t=" << start_age
		<< " initial orbit=";
	for(size_t i=0; i<nargs; ++i) {
		if(i) std::cerr << ", ";
		std::cerr << orbit[i];
	}
	std::cerr << std::endl;
#endif

//	const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_bsimp;
	const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;

	gsl_odeiv2_step *step=gsl_odeiv2_step_alloc(step_type, nargs);
	gsl_odeiv2_control *step_control=gsl_odeiv2_control_standard_new(
			precision, precision, 1, 0);
	gsl_odeiv2_evolve *evolve=gsl_odeiv2_evolve_alloc(nargs);

	void *sys_mode_windstate_lock_dissipation[5]={&system, &evolution_mode,
		&wind_state, const_cast<SpinOrbitLockInfo*>(&star_lock),
		&__dissipation};
	gsl_odeiv2_system ode_system={stellar_system_diff_eq,
		stellar_system_jacobian, nargs, sys_mode_windstate_lock_dissipation};
	double t=start_age; 
	std::valarray<double> derivatives(nargs), param_derivatives(nargs),
		age_derivatives(nargs);
	stellar_system_diff_eq(t, &(orbit[0]), &(derivatives[0]),
			sys_mode_windstate_lock_dissipation);
	append_to_orbit(orbit, derivatives, evolution_mode, star_lock, 
			wind_state, t, system, planet_formation_semimajor,
			planet_formation_inclination);

	update_stop_condition_history(t, orbit, derivatives, system,
			evolution_mode, stop_cond, stop_reason);
	clear_discarded();
	double step_size=0.01*(max_age-start_age);

	stop_reason=NO_STOP;
	StopInformation stop;
	while(t<max_age) {
		double max_next_t=std::min(t + max_step, max_age); 
		int status=GSL_SUCCESS;
		bool step_rejected=false;
		do {
			status=gsl_odeiv2_evolve_apply(
					evolve, step_control, step, &ode_system,
					&t, max_next_t, &step_size, &(orbit[0]));
			if (status != GSL_SUCCESS) {
				std::ostringstream msg;
				msg << "GSL signaled failure while evolving (error code " <<
					status << ")";
				throw Error::Runtime(msg.str());
			}
			stellar_system_diff_eq(t, &(orbit[0]), &(derivatives[0]),
					sys_mode_windstate_lock_dissipation);
			stop=update_stop_condition_history(t, orbit, derivatives, system,
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
		} while(step_rejected &&
				std::abs(stop.stop_condition_precision())>precision);
		if(!step_rejected)
			append_to_orbit(orbit, derivatives,
				evolution_mode, star_lock, wind_state, t, system,
                planet_formation_semimajor, planet_formation_inclination);
		if(stop.is_crossing() && stop.stop_reason()!=NO_STOP) {
			stop_reason=stop.stop_reason();
			stop.crossed_zero()=(t>stop.stop_age());
			break;
		}
	}
	max_age=t;
	clear_history();

	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(step_control);
	gsl_odeiv2_step_free(step);
	return stop;
}

CombinedStoppingCondition *OrbitSolver::get_stopping_condition(
		EvolModeType evolution_mode, const SpinOrbitLockInfo &star_lock)
	const
{
	CombinedStoppingCondition *result=new CombinedStoppingCondition();
	if(evolution_mode==LOCKED_TO_DISK) {
		*result|=new NoStopCondition();
		return result;
	}
	(*result)|=new WindSaturationCondition;
#ifdef EXTERNAL_CONDITION
	(*result)|= new EXTERNAL_CONDITION;
#endif
	if(evolution_mode==NO_PLANET) return result;
	(*result)|=new PlanetDeathCondition();
#ifdef DEBUG
	assert(evolution_mode==BINARY);
#endif
	if(star_lock)
		*result|=new BreakLockCondition();
	else {
/*		for(unsigned i=0; i<__dissipation.num_harmonics(); ++i)
			*result|=new SynchronizedCondition(
				__dissipation.harmonic(0, i).orbital_frequency_multiplier(),
				__dissipation.harmonic(0, i).spin_frequency_multiplier(), 0);*/
	}
	return result;
}

EvolModeType OrbitSolver::critical_age_evol_mode(double age, 
		const std::valarray<double> &orbit,
		double initial_semimajor, const StellarSystem &system, 
		EvolModeType evolution_mode, SpinOrbitLockInfo &star_lock,
		double planet_formation_age)
{
	if(system.get_star().disk_dissipation_age()>age)
		return LOCKED_TO_DISK;
	else if(planet_formation_age>age || 
			(planet_formation_age<age && evolution_mode==NO_PLANET))
		return NO_PLANET;
	const Star &star=system.get_star();
	double in_sync=Inf;
	std::valarray<double> dummy_cond_deriv;
	if(star_lock) in_sync=0;
	else {
		double wconv, worb=orbital_angular_velocity(
				system.get_planet().mass(), star.mass(), initial_semimajor);
		if(evolution_mode==LOCKED_TO_DISK) {
			wconv=star.disk_lock_frequency();
		} else if(evolution_mode==NO_PLANET) {
			wconv=orbit[0]/star.moment_of_inertia(age, convective);
		} else {
#ifdef DEBUG
					assert(evolution_mode==BINARY);
#endif
					wconv=orbit[2]/star.moment_of_inertia(age, convective);
		}
		__dissipation.init_harmonics(0, wconv/worb);
		star_lock=SpinOrbitLockInfo();
	}
	if(std::abs(in_sync)<precision) {
		throw Error::NotImplemented("Starting planet synchronized");
/*		BreakLockCondition locked_cond;
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
		double locked_cond_value=locked_cond(age, locked_orbit, locked_deriv,
				system, dummy_cond_deriv, evolution_mode)[0];
		if(locked_cond_value>=0)
			return LOCKED_TO_PLANET;
		else return (no_planet_dwconv_dt(age, locked_orbit, system) > 0 ?
				SLOW_PLANET : FAST_PLANET);*/
	}
	return BINARY;
}

EvolModeType OrbitSolver::next_evol_mode(double age,
		const std::valarray<double> &parameters,
		double initial_semimajor, StellarSystem &system,
		EvolModeType evolution_mode,
		SpinOrbitLockInfo &star_lock, WindSaturationState wind_state,
		StoppingConditionType condition_type,
		double condition_value, bool stopped_before,
		double planet_formation_age)
{
	if(stopped_before) condition_value*=-1;
	if(condition_type==NO_STOP) return critical_age_evol_mode(age,
			parameters, initial_semimajor, system, evolution_mode, star_lock,
			planet_formation_age);
	if(condition_type==SYNCHRONIZED) {
#ifdef DEBUG
		assert(star_lock);
		assert(evolution_mode==BINARY);
#endif
		system(age, parameters[2], parameters[1]);
		double semimajor=std::pow(parameters[0], 1.0/6.5);
		get_tidal_dissipation(system, age, &semimajor, star_lock,
				__dissipation);
		std::valarray<double> no_planet_diff_eq(3);
		system.differential_equations(&parameters[2], NO_PLANET, wind_state,
				star_lock, __dissipation, &no_planet_diff_eq[0]);
		double above_lock_fraction=system.above_lock_fraction(__dissipation, 
				no_planet_diff_eq[0]);
#ifdef DEBUG
		//#BEGIN_DEBUG_CODE#
		std::cerr << "Checking for " << star_lock
			<< ", at t=" << age << ", above_lock_fraction="
			<< above_lock_fraction << std::endl;
		//#END_DEBUG_CODE#
#endif
		if(above_lock_fraction<0) star_lock.lock_direction(-1);
		else if(above_lock_fraction>1) star_lock.lock_direction(1);
		else star_lock.lock_direction(0);
		__dissipation.init_harmonics(0, star_lock);
		return BINARY;
	} else if(condition_type==BREAK_LOCK) {
		star_lock.lock_direction(condition_value>0 ? 1 : -1);
		__dissipation.init_harmonics(0, star_lock);
		return BINARY;
	} else if(condition_type==PLANET_DEATH) return NO_PLANET;
	else if(condition_type==WIND_SATURATION || condition_type==EXTERNAL)
		return evolution_mode;
	else {
		std::ostringstream msg;
		msg << "Invalid stopping condition type encountered in "
			"OrbitSolver::next_evol_mode: " << condition_type;
		throw Error::BadFunctionArguments(msg.str());
	}
}

double OrbitSolver::stopping_age(double age, EvolModeType evolution_mode,
		const StellarSystem &system, double planet_formation_age,
		const std::list<double> &required_ages)
{
	double result;
	switch(evolution_mode) {
		case BINARY :
			result=end_age; break;
		case NO_PLANET :
			result=(std::isnan(planet_formation_age) || 
					age>planet_formation_age ? end_age :
					std::min(planet_formation_age, end_age));
			break;
		case LOCKED_TO_DISK :
			result=std::min(system.get_star().disk_dissipation_age(),
					end_age); break;
		default :
			throw Error::BadFunctionArguments("Unrecognized evolution "
					"mode in OrbitSolver::stopping_age");
	}
	double core_formation_age=system.get_star().core_formation_age();
	if(age < core_formation_age)
		result = std::min(result, core_formation_age);
	static std::list<double>::const_iterator
		next_required_age=required_ages.begin();
	if(age<=required_ages.front())
		next_required_age=required_ages.begin();
	if(next_required_age!=required_ages.end() && 
			result>*next_required_age) {
		if(age==*next_required_age) ++next_required_age;
		if(next_required_age!=required_ages.end()) 
			result=*next_required_age;
	}
	return result;
}

void OrbitSolver::parse_orbit_or_derivatives(EvolModeType evolution_mode,
		const SpinOrbitLockInfo &star_lock, double age,
		const std::valarray<double> &orbit_deriv, 
		const StellarSystem &system, bool evolution, double &a,
		double &inclination, double &Lconv, double &Lrad_parallel,
		double &Lrad_perpendicular) const
{
	const Star &star=system.get_star();
#ifdef DEBUG
	assert(evolution_mode!=TABULATION);
#endif
	if(evolution_mode==BINARY) {
		inclination=orbit_deriv[1];
		if(star_lock) {
#ifdef DEBUG
			assert(orbit_deriv.size()==4);
#endif
			double mplanet=system.get_planet().mass(),
				   mstar=star.mass(),
				   wconv=star_lock.spin(orbital_angular_velocity(mstar,
							   mplanet,
							   (evolution ? a : orbit_deriv[0])*Rsun_AU)),
				   Iconv=star.moment_of_inertia(age, convective);
			if(evolution) {
				double dIconv_dt=star.moment_of_inertia_deriv(age,
															  convective),
					   dwconv_da=star_lock.spin(orbital_angular_velocity(
								   mstar, mplanet, a*Rsun_AU, true));
				Lconv=dwconv_da*Iconv*orbit_deriv[0]*Rsun_AU+wconv*dIconv_dt;
			} else {
				Lconv=wconv*Iconv;
			}
			a=orbit_deriv[0];
		} else {
#ifdef DEBUG
			assert(orbit_deriv.size()==5);
#endif
			a=(evolution ? 1.0/6.5*std::pow(a, -5.5)*orbit_deriv[0] :
					std::pow(orbit_deriv[0], 1.0/6.5));
			Lconv=orbit_deriv[2];
		}
	} else if (evolution_mode==NO_PLANET) {
#ifdef DEBUG
		assert(orbit_deriv.size()==3);
#endif
		a = inclination = NaN;
		Lconv=orbit_deriv[0];
	} 
	if(evolution_mode==LOCKED_TO_DISK) {
#ifdef DEBUG
		assert(orbit_deriv.size()==1);
#endif
		a=inclination=NaN;
		Lconv=star.disk_lock_frequency()*
			(evolution ? star.moment_of_inertia_deriv(age, convective) :
			 star.moment_of_inertia(age, convective));
		Lrad_parallel=orbit_deriv[0];
		Lrad_perpendicular=0;
	} else {
		Lrad_parallel=orbit_deriv[orbit_deriv.size()-2];
		Lrad_perpendicular=orbit_deriv[orbit_deriv.size()-1];
	}
}

void OrbitSolver::collect_orbit_or_derivatives(EvolModeType evolution_mode,
		const SpinOrbitLockInfo &star_lock,
		double a, double inclination, double Lconv, double Lrad_parallel,
		double Lrad_perpendicular, std::valarray<double> &result,
		double semimajor) const
{
	if(evolution_mode==BINARY) {
		if(star_lock) {
			result.resize(4);
			result[0]=a;
		} else {
			result.resize(5);
			result[0]=(std::isnan(semimajor) ? std::pow(a, 6.5) :
					6.5*std::pow(semimajor, 5.5)*a);
			result[2]=Lconv;
		}
		result[1]=inclination;
	} else if(evolution_mode==NO_PLANET) {
		result.resize(3);
		result[0]=Lconv;
	} else if(evolution_mode==TABULATION) {
		result.resize(5);
		result[0]=a;
		result[1]=inclination;
		result[2]=Lconv;
	}

	if(evolution_mode==LOCKED_TO_DISK) 
		result.resize(Lrad_parallel, 1);
	else {
		result[result.size()-2]=Lrad_parallel;
		result[result.size()-1]=Lrad_perpendicular;
	}

}

std::valarray<double> OrbitSolver::transform_orbit(EvolModeType from_mode,
		EvolModeType to_mode, const SpinOrbitLockInfo &from_star_lock,
		const SpinOrbitLockInfo &to_star_lock, double age,
		const std::valarray<double> &from_orbit, 
		double initial_semimajor, double initial_inclination,
		const StellarSystem &system) const
{
	double a, inclination, Lconv, Lrad_parallel, Lrad_perpendicular;
	parse_orbit_or_derivatives(from_mode, from_star_lock, age, from_orbit,
			system, false, a, inclination, Lconv, Lrad_parallel,
			Lrad_perpendicular);
	if(to_mode!=TABULATION && std::isnan(a)) a=initial_semimajor/Rsun_AU;
	if(to_mode!=TABULATION && std::isnan(inclination)) 
		inclination=initial_inclination;

	std::valarray<double> to_orbit;
	collect_orbit_or_derivatives(to_mode, to_star_lock, a, inclination,
			Lconv, Lrad_parallel, Lrad_perpendicular, to_orbit);
	if(from_mode==BINARY && to_mode==NO_PLANET) {
		to_orbit[0]+=orbital_angular_momentum(
				system.get_star().mass(), 
				system.get_planet().mass(), a*Rsun_AU, 0);
	}

	return to_orbit;
}

std::valarray<double> OrbitSolver::transform_derivatives(
		EvolModeType from_mode, EvolModeType to_mode, 
		const SpinOrbitLockInfo &from_star_lock,
		const SpinOrbitLockInfo &to_star_lock, double age,
		const std::valarray<double> &from_orbit, 
		const std::valarray<double> &from_deriv,
		double initial_semimajor, const StellarSystem &system)
{
#ifdef DEBUG
	if(from_mode==NO_PLANET && to_mode!=TABULATION) 
		throw Error::BadFunctionArguments(
				"No evolution mode can possibly follow NO_PLANET in "
				"transform_derivatives.");
	assert(from_orbit.size()==from_deriv.size());
#endif

	double a, da_dt, dinclination_dt, dLconv_dt, dLrad_parallel_dt,
		   dLrad_perpendicular_dt;
	if(from_mode==BINARY) {
		a=da_dt=(from_star_lock ? from_orbit[0] :
				std::pow(from_orbit[0], 1.0/6.5));
	} else a=(to_mode==TABULATION ? NaN : initial_semimajor/Rsun_AU);
	parse_orbit_or_derivatives(from_mode, from_star_lock, age, from_deriv,
			system, true, da_dt, dinclination_dt, dLconv_dt,
			dLrad_parallel_dt, dLrad_perpendicular_dt);


	std::valarray<double> to_deriv;
	collect_orbit_or_derivatives(to_mode, to_star_lock, da_dt,
			dinclination_dt, dLconv_dt, dLrad_parallel_dt,
			dLrad_perpendicular_dt, to_deriv, a);
	return to_deriv;
}

void OrbitSolver::reset()
{
	__tabulated_ages.clear();
	__tabulated_evolution_mode.clear();
    __tabulated_wind_saturation.clear();
	__tabulated_lock.clear();
	for (int i=0; i < NUM_EVOL_VAR; i++) {
		__tabulated_orbit[i].clear();
		__tabulated_deriv[i].clear();
	}
}

OrbitSolver::OrbitSolver(double max_age, double required_precision) :
	end_age(std::abs(max_age)), precision(required_precision), 
	adjust_end_age(max_age<0), __tabulated_orbit(NUM_EVOL_VAR),
	__tabulated_deriv(NUM_EVOL_VAR)
{
	if (end_age > MAX_END_AGE) end_age = MAX_END_AGE;
}

void OrbitSolver::operator()(StellarSystem &system, double max_step,
		double planet_formation_age, double planet_formation_semimajor,
		double planet_formation_inclination, double start_age,
		EvolModeType initial_evol_mode,
		const SpinOrbitLockInfo &initial_lock,
		const std::valarray<double> &start_orbit,
		const std::list<double> &required_ages, bool no_evol_mode_change)
{
	const Planet &planet=system.get_planet();
	const Star &star=system.get_star();
	if(std::isnan(start_age)) start_age=star.core_formation_age();
	double stop_evol_age=(adjust_end_age ?
			std::min(end_age, star.lifetime()) : end_age);
	if(initial_evol_mode==BINARY) {
		if(start_orbit[0]<=planet.minimum_semimajor(start_age)/Rsun_AU)
			throw Error::BadFunctionArguments("Attempting to calculate the "
					"evolution of a system where the planet starts "
					"inside the Roche radius or its parent star!");
	} else if(std::isfinite(planet_formation_age) &&
			  planet_formation_age<stop_evol_age) {
		if(planet_formation_semimajor/Rsun_AU<=planet.minimum_semimajor(
					std::max(planet_formation_age,
						star.disk_dissipation_age()))/Rsun_AU)
			throw Error::BadFunctionArguments("Attempting to calculate the "
					"evolution of a system where the planet will be formed "
					"inside the Roche radius or its parent star!");
	}
	double wconv;
	switch(initial_evol_mode) {
		case BINARY :
#ifdef DEBUG
			assert(start_orbit.size()==(initial_lock ? 4 : 5)); 
#endif
			wconv=(initial_lock 
					? 
					initial_lock.spin(
						orbital_angular_velocity(planet.mass(), star.mass(),
							start_orbit[0]*Rsun_AU))
					:
					star.spin_frequency(start_age, convective,
						start_orbit[1])
					);
			break;
		case NO_PLANET :
			assert(start_orbit.size()==3);
			wconv=star.spin_frequency(start_age, convective, start_orbit[0]);
			break;
		case LOCKED_TO_DISK :
			assert(start_orbit.size()==1);
			wconv=star.disk_lock_frequency();
			break;
		default:
			throw Error::BadFunctionArguments("Unsupported initial evolution"
					" mode encountered in OrbitSolver::operator().");
	}
	WindSaturationState wind_state=(
			wconv<star.wind_saturation_frequency() ? NOT_SATURATED : 
			                                             SATURATED);

	reset();
	StoppingConditionType stop_reason=NO_STOP;
	double last_age=start_age;
	std::valarray<double> orbit=start_orbit;
	EvolModeType evolution_mode=initial_evol_mode;
	SpinOrbitLockInfo star_lock=initial_lock, old_star_lock=star_lock;
	while(last_age<stop_evol_age) {
		double next_stop_age=std::min(stopping_age(last_age, evolution_mode,
				system, planet_formation_age, required_ages), stop_evol_age);
		CombinedStoppingCondition *stopping_condition=get_stopping_condition(
				evolution_mode, star_lock);
		StopInformation stop_information=evolve_until(system, start_age,
				next_stop_age, orbit, stop_reason, max_step, evolution_mode,
				star_lock, wind_state, *stopping_condition,
				planet_formation_semimajor, planet_formation_inclination);
		last_age=next_stop_age;
		if(last_age<stop_evol_age) {
			EvolModeType old_evolution_mode=evolution_mode;
			old_star_lock=star_lock;
			if(stop_reason==SYNCHRONIZED) {
				const SynchronizedCondition &sync_cond=
					stopping_condition->get_sync_condition(
							stop_information.stop_condition_index());
				star_lock.set_lock(sync_cond.orbital_frequency_multiplier(),
						sync_cond.spin_frequency_multiplier(), 0);
			}
			if(!no_evol_mode_change) evolution_mode=next_evol_mode(last_age,
					orbit, planet_formation_semimajor, system,
					evolution_mode, star_lock, wind_state, stop_reason,
					stop_information.stop_condition_precision(),
					!stop_information.crossed_zero(), planet_formation_age);
#ifdef DEBUG
			std::cerr << "Changing evolution mode from "
				<< (old_star_lock ? "locked " : "not locked ") 
				<< old_evolution_mode << " to " 
				<< (star_lock ? "locked " : "not locked ") << evolution_mode
				<< std::endl;
#endif
			if(old_evolution_mode!=evolution_mode ||
					old_star_lock!=star_lock) {
				std::valarray<double> new_orbit=transform_orbit(
						old_evolution_mode, evolution_mode, old_star_lock,
						star_lock, last_age, orbit,
						planet_formation_semimajor,
						planet_formation_inclination, system);
#ifdef DEBUG
				std::cerr << "Transforming orbit from: " << orbit
					<< " to " << new_orbit << std::endl;
#endif
				orbit.resize(new_orbit.size());
				orbit=new_orbit;
			}
            if(evolution_mode==NO_PLANET && stop_reason!=WIND_SATURATION) {
                wconv=star.spin_frequency(last_age, convective,
                        orbit[0]); 
                wind_state=(
                    wconv<star.wind_saturation_frequency() ?
                    NOT_SATURATED : SATURATED);
            } else if(stop_reason==WIND_SATURATION)
				wind_state=static_cast<WindSaturationState>(-wind_state);
			start_age=last_age;
		}
		delete stopping_condition;
	}
}

const std::list<double> *OrbitSolver::get_tabulated_var(EvolVarType var_type)
	const
{
	if(var_type==AGE) return &__tabulated_ages;
	return &(__tabulated_orbit[var_type]);
}

const std::list<double> *OrbitSolver::get_tabulated_var_deriv(
		EvolVarType var_type) const
{
	assert(var_type!=AGE);
	return &(__tabulated_deriv[var_type]);
}
