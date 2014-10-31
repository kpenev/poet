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
		EvolModeType
#ifdef DEBUG
		evol_mode
#endif
		, const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode==BINARY);
	if(__system.number_locked_zones())
		assert(orbit[0]==__system.semimajor());
	else assert(std::pow(orbit[0], 1.0/6.5)==__system.semimajor());
	assert(orbit.size()==1 + 3*__system.number_zones() -
						 __system.number_locked_zones());
	assert(orbit.size()==derivatives.size());
#endif
	double m1=__system.primary().mass(), 
		   m2=__system.secondary().mass(),
		   semimajor=__system.semimajor(),
		   worb=orbital_angular_velocity(m1, m2, semimajor),
		   wspin=__zone.spin_frequency(),
		   dworb_dt=orbital_angular_velocity(m1, m2, semimajor, true)
			   *
			   derivatives[0];

	unsigned angmom_ind=1+2*__system.number_zones();
	if(!__primary) angmom_ind+=__system.primary().number_zones()
		-
			__system.primary.number_locked_zones();
	for(unsigned i=0; i<__zone_index; ++i)
		if(!__system.secondary().zone(i).locked()) ++angmom_ind;

	double dwspin_dt=(derivatives[angmom_ind]
					  -
					  __zone.moment_of_inertia(1)*wspin)
					/__zone.moment_of_inertia();
	if(__system.number_locked_zones()) 
		dworb_dt/=(6.5*orbit[0]/semimajor);
	stop_deriv.resize(1, 
			((wspin*dworb_dt-dwspin_dt*worb)/std::pow(worb, 2)*
			 __spin_freq_mult)/__orbital_freq_mult);
	return std::valarray<double>(
			(__orbital_freq_mult*worb - wspin*__spin_freq_mult)/worb, 1);
}

std::valarray<double> BreakLockCondition::operator()(
		EvolModeType
#ifdef DEBUG
		evol_mode
#endif
		,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode==BINARY);
	assert(orbit.size()==1 + 3*__system.number_zones() -
						 __system.number_locked_zones());
	assert(orbit.size()==derivatives.size());
#endif
	double frac=__system.above_lock_fraction(__locked_zone_index),
		   dfrac_dt=__system.above_lock_fraction(__locked_zone_index,
				   								 Dissipation::AGE)
			   		+
					__system.above_lock_fraction(__locked_zone_index,
												 Dissipation::RADIUS, 0,
												 false)
					*__system.primary().radius(1)
					+
					__system.above_lock_fraction(__locked_zone_index,
												 Dissipation::RADIUS, 0,
												 true)
					*__system.secondary().radius(1);
	unsigned deriv_zone_ind=0, unlocked_zone_ind=0;
	double *inclination_rate=derivatives+2,
		   *periapsis_rate=inclination_rate+__system.number_zones(),
		   *angmom_rate=periapsis_rate+__system.number_zones()-1;
	for(int i=0; i<2; ++i) {
		DissipatingBody &body=(i==0 ? __system.primary()
									: __system.secondary());
		for(unsigned zone_ind=0; zone_ind<body.number_zones(); ++i) {
			double dangmom_dt;
			DissipatingZone &zone=body.zone(zone_ind);
			if(!zone.locked()) dangmom_dt=angmom_rate[unlocked_zone_ind++];
			else {
				double dworb_dt=orbital_angular_velocity(m1, m2, semimajor,
														 true)
								*derivatives[0];
				dangmom_dt=zone.moment_of_inertia(1)*zone.spin_frequency()
						   +
						   zone.moment_of_inertia()
						   *zone.lock_held().spin(dworb_dt);
			}
			dfrac_dt+=__system.above_lock_fraction(
					__locked_zone_index, Dissipation::MOMENT_OF_INERTIA,
					deriv_zone_ind)*zone.moment_of_inertia(1)
				+
				__system.above_lock_fraction(
						__locked_zone_index, Dissipation::INCLINATION,
						deriv_zone_ind)*inclination_rate[deriv_zone_ind];
				+
				__system.above_lock_fraction(
						__locked_zone_index, Dissipation::SPIN_ANGMOM,
						deriv_zone_ind)*inclination_rate[deriv_zone_ind]
				*dangmom_dt;
			if(deriv_zone_ind) dfract_dt+=__system.above_lock_fraction(
					__locked_zone_index, Dissipation::PERIAPSIS,
					deriv_zone_ind)*periapsis_rate[deriv_zone_ind-1];
			++deriv_zone_ind;
		}

	}
	stop_deriv.resize(2, dfrac_dt);
	std::valarray<double> result(2);
	result[0]=frac;
	result[1]=frac-1.0;
	return result;
}

std::valarray<double> SecondaryDeathCondition::operator()(
		EvolModeType
#ifdef DEBUG
		evol_mode
#endif
		,
		const std::valarray<double> &orbit, 
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode==BINARY);
	assert(orbit.size()==1 + 3*__system.number_zones() -
						 __system.number_locked_zones());
	assert(orbit.size()==derivatives.size());
#endif
	double min_semimajor=__system.minimum_semimajor(),
		   semimajor=__system.semimajor(),
		   dsemimajor_dt=derivatives[0];
	if(__system.number_locked_zones()==0) 
		dsemimajor_dt*=semimajor/(6.5*orbit[0]);
	stop_deriv.resize(1,(dsemimajor_dt*min_semimajor
						  -semimajor*__system.minimum_semimajor(true))
						/std::pow(min_semimajor, 2));
	return std::valarray<double>(
			(semimajor-min_semimajor)/min_semimajor, 1);
}

std::valarray<double> WindSaturationCondition::operator()(
		EvolModeType evol_mode,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
#ifdef DEBUG
	assert(evol_mode!=LOCKED_SURFACE_SPIN);
	if(evole_mode!=BINARY) assert(__primary);
#endif 
	unsigned num_zones=__body.number_zones();
	if(evol_mode==BINARY) num_zones+=__other_body.number_zones();
	double wsurf=__body.spin_frequency(),
		   surf_angmom_deriv=
			   stop_deriv[1 + 2*num_zones
			   			  +
						  (__primary ? 0
						   			 : __other_body.number_zones()
									   -__other_body.number_locked_zones())
						 ];
	stop_deriv.resize(1,
		   	(surf_angmom_deriv - __body.zone(0).moment_of_inertia(1)*wsurf)
			/(__body.zone(0).moment_of_inertia()*wsat));
	if(std::isinf(wsat)) return std::valarray<double>(-1.0, 1);
	return std::valarray<double>((wsurf-wsat)/wsat, 1);
}

void CombinedStoppingCondition::update_meta_information(
		const StoppingCondition *rhs)
{
	size_t num_new_subcond=rhs->num_subconditions();
	__num_subconditions+=num_new_subcond;
	__types.reserve(__types.size()+num_new_subcond);
	for(size_t i=0; i<num_new_subcond; i++) __types.push_back(rhs->type(i));
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
	__sub_conditions.insert(__sub_conditions.end(),
			rhs.__sub_conditions.begin(), rhs.__sub_conditions.end());
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
	__sub_conditions.push_back(rhs);
	return *this;
}

std::valarray<double> CombinedStoppingCondition::operator()(
		EvolModeType evol_mode,
		const std::valarray<double> &orbit,
		const std::valarray<double> &derivatives,
		std::valarray<double> &stop_deriv) const
{
	std::valarray<double> result(__num_subconditions);
	stop_deriv.resize(__num_subconditions, NaN);
	size_t i=0;
	for(std::vector<const StoppingCondition *>::const_iterator
			cond=__sub_conditions.begin(); cond!=__sub_conditions.end();
			++cond)
		add_subcondition_values(*cond, orbit, derivatives, dissipation,
				system, evol_mode, i, result, stop_deriv);
	return result;
}

void CombinedStoppingCondition::reached(unsigned index, short deriv_sign)
{
	std::vector<const StoppingCondition *>::iterator sc_iter=
		__sub_conditions.begin();
	while(index>=sc_iter->num_subconditions()) {
#ifdef DEBUG
		assert(sc_iter!=__sub_conditions.end());
#endif
		index-=sc_iter->num_subconditions();
		++sc_iter;
	}
	sc_iter->reached(index, deriv_sign);
}

CombinedStoppingCondition::~CombinedStoppingCondition()
{
	if(!__delete_subcond) return;
	for(std::vector<const StoppingCondition *>::const_iterator
			i=__sub_conditions.begin(); i!=__sub_conditions.end(); ++i)
		delete *i;
}


