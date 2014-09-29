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

std::valarray<double> SynchronizedCondition::operator()(
		const std::valarray<double> &,
		const std::valarray<double> &derivatives,
		const TidalDissipation &tidal_dissipation,
		const StellarSystem &system,
		std::valarray<double> &stop_deriv,
		EvolModeType) const
{
	const Star &star=system.get_star();
	const Planet &planet=system.get_planet();
	double worb=tidal_dissipation.orbit_frequency(),
		   Lspin=tidal_dissipation.spin_angular_momentum(__body_index),
		   wspin=(__body_index==0 ? star.spin_frequency(envelope, Lspin) :
				   planet.spin_frequency(system.age(), Lspin));

	if(derivatives.size()>0) {
		double dwspin_dt=
			tidal_dissipation(__body_index, Dissipation::TORQUEZ)*
			(__body_index==0
			 ? star.spin_frequency_angmom_deriv(convective, Lspin) 
			 : planet.spin_frequency_angmom_deriv(system.age(), Lspin)),
			   dworb_dt=tidal_dissipation(0, Dissipation::ORBIT_SPINUP)
				   +
				   tidal_dissipation(1, Dissipation::ORBIT_SPINUP);
		stop_deriv.resize(1, 
				((wspin*dworb_dt-dwspin_dt*worb)/std::pow(worb, 2)*
				 __spin_freq_mult)/__orbital_freq_mult);
	} else stop_deriv.resize(1, NaN);
	return std::valarray<double>(
			(__orbital_freq_mult*worb - wspin*__spin_freq_mult)/worb, 1);
}

std::valarray<double> BreakLockCondition::operator()(
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
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
	assert(derivatives.size()==4);
#endif
	const Star &star=system.get_star();
	double mplanet=system.get_planet().mass(),
		   wconv=orbital_angular_velocity(star.mass(), mplanet,
				   orbit[0]),
		   Lconv=wconv*star.moment_of_inertia(convective),
		   wind_torque=star.wind_torque(wconv, 0),
		   coupling_z_torque=star.differential_rotation_torque_angmom(Lconv,
				   std::complex<double>(orbit[2], orbit[3])).real(),
		   above_fraction=system.above_lock_fraction(tidal_dissipation,
				   coupling_z_torque - wind_torque);
	stop_deriv.resize(2, NaN);
	std::valarray<double> result(2);
	result[0]=above_fraction;
	result[1]=above_fraction-1.0;
	return result;
}

std::valarray<double> PlanetDeathCondition::operator()(
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

	double min_semimajor=system.get_planet().minimum_semimajor(system.age());
	stop_deriv.resize(1, NaN);
	if(std::isinf(min_semimajor)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>(
			(tidal_dissipation.semimajor()-min_semimajor)/min_semimajor, 1);
}

double convective_frequency(const StellarSystem &system, 
		const std::valarray<double> &orbit, EvolModeType evol_mode)
{
	if(evol_mode==BINARY) {
		if(orbit.size()==4)
			return orbital_angular_velocity(system.get_star().mass(),
					system.get_planet().mass(), orbit[0]);
		else {
#ifdef DEBUG
			assert(orbit.size()==5);
#endif
			return system.get_star().spin_frequency(convective, orbit[2]);
		}
	} else if(evol_mode==NO_PLANET)
		return system.get_star().spin_frequency(convective, orbit[0]);
	else if(evol_mode==LOCKED_TO_DISK)
		throw Error::BadFunctionArguments("The convective spin frequency"
				" is not defined for disk locked evolution mode.");
	else throw Error::BadFunctionArguments("Unrecognized wind saturation "
			"argument");
}

std::valarray<double> WindSaturationCondition::operator()(
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
		   wconv=convective_frequency(system, orbit, evol_mode);

	stop_deriv.resize(1, NaN);
	if(std::isinf(wsat)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>((wconv-wsat)/wsat, 1);
}

void CombinedStoppingCondition::update_meta_information(
		const StoppingCondition *rhs)
{
	size_t num_new_subcond=rhs->num_subconditions();
	__num_subconditions+=num_new_subcond;
	__types.reserve(__types.size()+num_new_subcond);
	__interpolation_ranges.reserve(
			__interpolation_ranges.size()+num_new_subcond);
	for(size_t i=0; i<num_new_subcond; i++) {
		__types.push_back(rhs->type(i));
		__interpolation_ranges.push_back(rhs->interpolation_range(i));
	}
}

void CombinedStoppingCondition::add_subcondition_values(
		const StoppingCondition *cond, 
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		const TidalDissipation &dissipation,
		const StellarSystem &system,
		EvolModeType evol_mode,
		size_t &first_index,
		std::valarray<double> &values,
		std::valarray<double> &derivs) const
{
	std::valarray<double> sub_stop_deriv(cond->num_subconditions()),
		temp_array=(*cond)(orbit, derivatives, dissipation, system,
				sub_stop_deriv, evol_mode);
	for(size_t subcond_ind=0; subcond_ind<temp_array.size();
			subcond_ind++) {
		values[first_index]=temp_array[subcond_ind];
		derivs[first_index]=sub_stop_deriv[subcond_ind];
		++first_index;
	}
}

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const CombinedStoppingCondition &rhs)
{
	update_meta_information(&rhs);
	__generic_sub_conditions.insert(__generic_sub_conditions.end(),
			rhs.__generic_sub_conditions.begin(),
			rhs.__generic_sub_conditions.end());
	__sync_sub_conditions.insert(__sync_sub_conditions.end(),
			rhs.__sync_sub_conditions.begin(),
			rhs.__sync_sub_conditions.end());
	if(__delete_subcond) {
#ifdef DEBUG
		assert(rhs.__delete_subcond);
#endif
		const_cast<CombinedStoppingCondition &>(rhs).__delete_subcond=false;
	}
	return *this;
}

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const StoppingCondition *rhs)
{
	update_meta_information(rhs);
	__generic_sub_conditions.push_back(rhs);
	return *this;
}

CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
		const SynchronizedCondition *rhs)
{
	update_meta_information(rhs);
	__sync_sub_conditions.push_back(rhs);
	return *this;
}

CombinedStoppingCondition::~CombinedStoppingCondition()
{
	if(!__delete_subcond) return;
	for(std::vector<const StoppingCondition *>::const_iterator
			i=__generic_sub_conditions.begin();
			i!=__generic_sub_conditions.end(); ++i)
		delete *i;
	for(std::vector<const SynchronizedCondition *>::const_iterator
			i=__sync_sub_conditions.begin();
			i!=__sync_sub_conditions.end(); ++i)
		delete *i;
}

std::valarray<double> CombinedStoppingCondition::operator()(
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
			cond=__generic_sub_conditions.begin();
			cond!=__generic_sub_conditions.end(); ++cond)
		add_subcondition_values(*cond, orbit, derivatives, dissipation,
				system, evol_mode, i, result, stop_deriv);
	for(std::vector<const SynchronizedCondition *>::const_iterator
			cond=__sync_sub_conditions.begin();
			cond!=__sync_sub_conditions.end(); ++cond) 
		add_subcondition_values(*cond, orbit, derivatives, dissipation,
				system, evol_mode, i, result, stop_deriv);
	return result;
}
