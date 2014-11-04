#include "SynchronizedCondition.h"
#include "BinarySystem.h"
#include "DissipatingZone.h"

SynchronizedCondition::SynchronizedCondition(int orbital_freq_mult,
		int spin_freq_mult, short deriv_sign, bool primary,
		unsigned zone_index, BinarySystem &system) :
		__orbital_freq_mult(orbital_freq_mult),
		__spin_freq_mult(spin_freq_mult), __primary(primary),
		__zone_index(zone_index), 
		__expected_crossing_deriv_sign(deriv_sign),
		__zone((primary 
				? system.primary() : system.secondary()).zone(zone_index)),
		__system(system)
{}

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
			__system.primary().number_locked_zones();
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

void SynchronizedCondition::reached(short deriv_sign, unsigned index)
{
#ifdef DEBUG
	assert(deriv_sign==__expected_crossing_deriv_sign);
#endif
	StoppingCondition::reached(deriv_sign, index);
	__system.check_for_lock(__orbital_freq_mult, __spin_freq_mult,
							(__primary ? 0 : 1), __zone_index);
}
