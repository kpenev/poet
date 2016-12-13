#include "EvolvingStellarZone.h"

void EvolvingStellarZone::reset_current_quantities()
{
	for(size_t i=0; i<__current_age_quantities.size(); ++i)
		if(__current_age_quantities[i]) {
			delete __current_age_quantities[i];
			__current_age_quantities[i]=NULL;
		}
}

double EvolvingStellarZone::current_age_quantity(size_t quantity,
                                                 unsigned deriv_order) const
{
    assert(quantity < __evolving_quantities.size());

	if(__current_age_quantities[quantity]==NULL) {
        __current_age_quantities[quantity] = 
            __evolving_quantities[quantity]->deriv(__current_age);
	}
	return __current_age_quantities[quantity]->order(deriv_order);
}

double EvolvingStellarZone::any_age_quantity(size_t quantity,
                                             double age, 
                                             unsigned deriv_order=0) const
{
    if(deriv_order == 0)
        return (*(__evolving_quantities[quantity]))(age);
    else {
        FunctionDerivatives 
            *deriv = __evolving_quantities[quantity]->deriv(age);
        double result = deriv->order(deriv_order);
        delete deriv;
        return result;
    }
}

void EvolvingStellarZone::configure(double age,
                                    double orbital_frequency,
                                    double eccentricity,
                                    double orbital_angmom,
                                    double spin,
                                    double inclination,
                                    double periapsis,
                                    bool spin_is_frequency)
{
	if(__current_age != age) {
		__current_age = age;
		reset_current_quantities();
	}
	DissipatingZone::configure(age,
                               orbital_frequency,
                               eccentricity,
							   orbital_angmom,
                               spin,
                               inclination,
                               periapsis,
							   spin_is_frequency);
}

EvolvingStellarZone::~EvolvingStellarZone() 
{
    reset_current_quantities();

    for(size_t i = 0; i < __evolving_quantities.size(); ++i)
        if(__evolving_quantities[i]) delete __evolving_quantities[i];
}

void EvolvingStellarZone::reached_critical_age(double age)
{
    for(size_t i = 0; i < __evolving_quantities.size(); ++i)
        __evolving_quantities[i]->reached_critical_age(age);
}

double EvolvingStellarZone::next_stop_age() const
{
    double result = Inf;
    for(size_t i = 0; i < __evolving_quantities.size(); ++i)
        result = std::min(result,
                          __evolving_quantities[i]->next_stop_age(age));
    return result;
}