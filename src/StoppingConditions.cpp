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
		const std::valarray<double> &,
		const std::valarray<double> &derivatives,
		const TidalDissipation &tidal_dissipation,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType) const
{
	double worb=tidal_dissipation.orbital_frequency(),
		   Lspin=tidal_dissipation.spin_angular_momentum(__body_index),
		   wspin=(__body_index==0 ? 
				   system.get_star().spin_frequency(age, envelope, Lspin) :
				   system.get_planet().spin_frequency(age, Lspin));

	if(derivatives.size()>0) {
		double dwspin_dt=
			tidal_dissipation(__body_index, Dissipation::TORQUEZ)*
			(__body_index==0
			 ? 
			 system.get_star().spin_frequency_angmom_deriv(age, convective,
				 Lspin) 
			 :
			 system.get_planet().spin_frequency_angmom_deriv(age, Lspin)),
			   dworb_dt=tidal_dissipation(0, Dissipation::ORBIT_SPINUP)
				   +
				   tidal_dissipation(1, Dissipation::ORBIT_SPINUP);
		stop_deriv.resize(1, 
				((wspin*dworb_dt-dwspin_dt*worb)/std::pow(worb, 2)*
				 __spin_freq_mult)/__orbital_freq_mult);
	} else stop_deriv.resize(1, NaN);
	return std::valarray<double>(
			(worb - (wspin*__spin_freq_mult)/__orbital_freq_mult)/worb, 1);
}

std::valarray<double> BreakLockCondition::operator()(double,
		const std::valarray<double> &,
		const std::valarray<double> &derivatives,
		const TidalDissipation &tidal_dissipation,
		const StellarSystem &, 
		std::valarray<double> &stop_deriv,
		EvolModeType
#ifdef DEBUG
		evol_mode
#endif
		) const
{
#ifdef DEBUG
	assert(evol_mode==BINARY);
	assert(derivatives.size()==4);
#endif
	double below_da_dt=0, above_da_dt=0, da_dt=derivatives[0];
	for(short body_ind=0; body_ind<=1; ++body_ind) {
		double non_locked_da_dt=tidal_dissipation(body_ind,
				Dissipation::SEMIMAJOR_DECAY);
		below_da_dt+=tidal_dissipation(body_ind,
				Dissipation::SEMIMAJOR_DECAY, Dissipation::NO_DERIV, -1)
			+
			non_locked_da_dt;
		above_da_dt=tidal_dissipation(body_ind,
				Dissipation::SEMIMAJOR_DECAY, Dissipation::NO_DERIV,
				1)
			+
			non_locked_da_dt;
	}
	stop_deriv.resize(2, NaN);
	if(da_dt==0) return std::valarray<double>(Inf, 2);
	std::valarray<double> result(2);
	result[0]=(std::max(below_da_dt, above_da_dt)-da_dt)/da_dt;
	result[1]=(da_dt-std::min(below_da_dt, above_da_dt))/da_dt;
	return result;
}

std::valarray<double> PlanetDeathCondition::operator()(double age,
		const std::valarray<double> &,
		const std::valarray<double> &,
		const TidalDissipation &tidal_dissipation,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType
#ifdef DEBUG
		evol_mode
#endif
		) const
{
#ifdef DEBUG
	assert(evol_mode==BINARY);
#endif

	double min_semimajor=system.get_planet().minimum_semimajor(age);
	stop_deriv.resize(1, NaN);
	if(std::isinf(min_semimajor)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>(
			(tidal_dissipation.semimajor()-min_semimajor)/min_semimajor, 1);
}

double convective_frequency(double age, const StellarSystem &system, 
		const std::valarray<double> &orbit, EvolModeType evol_mode)
{
	if(evol_mode==BINARY) {
		if(orbit.size()==4)
			return orbital_angular_velocity(system.get_star().mass(),
					system.get_planet().mass(), orbit[0]*Rsun_AU);
		else {
#ifdef DEBUG
			assert(orbit.size()==5);
#endif
			return system.get_star().spin_frequency(age, convective,
				orbit[2]);
		}
	} else if(evol_mode==NO_PLANET)
		return system.get_star().spin_frequency(age, convective,
				orbit[0]);
	else if(evol_mode==LOCKED_TO_DISK)
		throw Error::BadFunctionArguments("The convective spin frequency"
				" is not defined for disk locked evolution mode.");
	else throw Error::BadFunctionArguments("Unrecognized wind saturation "
			"argument");
}

std::valarray<double> WindSaturationCondition::operator()(double age,
		const std::valarray<double> &orbit,
		//Accepts but never uses derivatives
		const std::valarray<double> &,
		//Accepts but never uses tidal dissipation
		const TidalDissipation &,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType evol_mode) const
{
	double wsat=system.get_star().wind_saturation_frequency(),
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
		const TidalDissipation &dissipation,
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
			temp_array=(**cond)(age, orbit, derivatives, dissipation, system,
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
