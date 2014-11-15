#include "YRECCore.h"

void YRECCore::reset_current_quantities()
{
	for(size_t i=0; i<__current_age_quantities.size(); ++i)
		if(__current_age_quantities[i]) {
			delete __current_age_quantities[i];
			__current_age_quantities[i]=NULL;
		}
}

void YRECCore::reset()
{
	reset_current_quantities();
	if(__mass) {
		delete __mass;
		__mass=NULL;
	}
	if(__radius) {
		delete __radius;
		__radius=NULL;
	}
	if(__moment_of_inertia) {
		delete __moment_of_inertia;
		__moment_of_inertia=NULL;
	}
}

double YRECCore::current_age_quantity(CurrentAgeQuantities quantity,
		unsigned deriv_order) const
{
	if(__current_age<__formation_age) return 0;
	if(__current_age_quantities[quantity]==NULL) {
		switch(quantity) {
			case MASS :
				__current_age_quantities[MASS]=
					__mass->deriv(__current_age);
				break;
			case RADIUS :
				__current_age_quantities[RADIUS]=
					__radius->deriv(__current_age);
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

void YRECCore::configure(double age, double orbital_frequency,
		double eccentricity, double orbital_angmom, double spin_angmom,
		double inclination, double periapsis)
{
	if(__current_age!=age) {
		__current_age=age;
		reset_current_quantities();
	}
	DissipatingZone::configure(age, orbital_frequency, eccentricity,
							   orbital_angmom, spin_angmom, inclination,
							   periapsis);
}
