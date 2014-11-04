#include "YRECEnvelope.h"

void YRECEnvelope::reset_current_quantities()
{
	for(size_t i=0; i<__current_age_quantities.size(); ++i)
		if(__current_age_quantities[i]) {
			delete __current_age_quantities[i];
			__current_age_quantities[i]=NULL;
		}
}

void YRECEnvelope::reset()
{
	reset_current_quantities();
	if(__outer_radius) {
		delete __outer_radius;
		__outer_radius=NULL;
	}
	if(__moment_of_inertia) {
		delete __moment_of_inertia;
		__moment_of_inertia=NULL;
	}
}

double YRECEnvelope::current_age_quantity(CurrentAgeQuantities quantity,
		unsigned deriv_order) const
{
	if(__current_age_quantities[quantity]==NULL) {
		switch(quantity) {
			case RADIUS :
				__current_age_quantities[RADIUS]=
					__outer_radius->deriv(__current_age);
				break;
			case INERTIA :
				__current_age_quantities[INERTIA]=
					__moment_of_inertia->deriv(__current_age);
				break;
			default :
#ifdef DEBUG
				assert(false)
#endif
					;
		}
	}
	return __current_age_quantities[quantity]->order(deriv_order);
}

void YRECEnvelope::create(double stellar_mass, 
		const EvolvingStellarQuantity *outer_radius,
		const EvolvingStellarQuantity *moment_of_inertia)
{
	reset();
	__stellar_mass=stellar_mass;
	__outer_radius=outer_radius;
	__moment_of_inertia=moment_of_inertia;
}

void YRECEnvelope::configure(double age, double orbital_frequency,
		double eccentricity, double orbital_angmom, double spin_angmom,
		double inclination, double periapsis)
{
	DissipatingZone::configure(age, orbital_frequency, eccentricity,
							   orbital_angmom, spin_angmom, inclination,
							   periapsis);
	if(__current_age!=age) {
		__current_age=age;
		reset_current_quantities();
	}
}
