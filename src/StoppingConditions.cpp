/**\file
 * 
 * \brief The implementations of the various stopping condition methods.
 *
 * \ingroup OrbitSolver_group
 */

#include "StoppingConditions.h"

std::ostream &operator<<(std::ostream &os,
		const StoppingConditionType &stop_cond_type)
{
	switch(stop_cond_type) {
		case NO_STOP: os << "NO_STOP"; break;
		case SYNCHRONIZED: os << "SYNCHRONIZED"; break;
		case BREAK_LOCK: os << "BREAK_LOCK"; break;
		case PLANET_DEATH: os << "PLANET_DEATH"; break;
		case WIND_SATURATION: os << "WIND_SATURATION"; break;
		case EXTERNAL: os << "EXTERNAL";
	}
	return os;
}

std::valarray<double> SynchronizedCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double a=(evol_mode==LOCKED_TO_DISK || evol_mode==NO_PLANET ?
            __initial_semimajor : std::pow(orbit[0], 1.0/6.5)*Rsun_AU),
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

std::valarray<double> BreakLockCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system, 
		std::valarray<double> &stop_deriv,
		EvolModeType) const
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
		//Accepts but never uses derivatives
		const std::valarray<double> &,
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
		//Accepts but never uses derivatives
		const std::valarray<double> &,
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

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const CombinedStoppingCondition &rhs)
{
	__sub_conditions.insert(__sub_conditions.end(),
			rhs.__sub_conditions.begin(), rhs.__sub_conditions.end());
	__num_subconditions+=rhs.__num_subconditions;
	if(__delete_subcond)
		const_cast<CombinedStoppingCondition &>(rhs).__delete_subcond=false;
	__types.insert(__types.end(), rhs.__types.begin(), rhs.__types.end());
	__interpolation_ranges.insert(__interpolation_ranges.end(),
			rhs.__interpolation_ranges.begin(),
			rhs.__interpolation_ranges.end());
	return *this;
}

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const StoppingCondition *rhs)
{
	size_t num_new_subcond=rhs->num_subconditions();
	__num_subconditions+=num_new_subcond;
	__sub_conditions.push_back(rhs);
	__types.reserve(__types.size()+num_new_subcond);
	__interpolation_ranges.reserve(
			__interpolation_ranges.size()+num_new_subcond);
	for(size_t i=0; i<num_new_subcond; i++) {
		__types.push_back(rhs->type(i));
		__interpolation_ranges.push_back(rhs->interpolation_range(i));
	}
	return *this;
}

CombinedStoppingCondition::~CombinedStoppingCondition()
{
	for(std::vector<const StoppingCondition *>::const_iterator
			i=__sub_conditions.begin(); i!=__sub_conditions.end(); i++)
		delete *i;
}

std::valarray<double> CombinedStoppingCondition::operator()(double age,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	std::valarray<double> result(__num_subconditions);
	stop_deriv.resize(__num_subconditions, NaN);
	size_t i=0;
	for(std::vector<const StoppingCondition *>::const_iterator
			cond=__sub_conditions.begin(); cond!=__sub_conditions.end();
			cond++) {
		std::valarray<double> sub_stop_deriv(1),
			temp_array=(**cond)(age, orbit, derivatives, system,
					sub_stop_deriv, evol_mode);
		for(size_t subcond_ind=0; subcond_ind<temp_array.size();
				subcond_ind++) {
			result[i]=temp_array[subcond_ind];
			stop_deriv[i]=sub_stop_deriv[subcond_ind];
			i++;
		}
	}
	return result;
}
