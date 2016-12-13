#ifndef __EVOLVING_STELLAR_ZONE_H
#define __EVOLVING_STELLAR_ZONE_H

#ifdef TWO_QS
#include "TwoPhaseLagZone.h"
#else
#include "PowerlawPhaseLagZone.h"
#endif

#include "Functions.h"
#include "EvolvingStellarQuantity.h"

#include <initializer_list>

///\brief Base class for zones of low mass evolving stars.
class EvolvingStellarZone :
#ifdef TWO_QS
    virtual public TwoPhaseLagZone
#else
    virtual public PowerlawPhaseLagZone
#endif
{
private:
	///The age for the last configure() call.
	double __current_age;

	///\brief Pre-computed values and derivatives for quantities which only
	///depend on age at the current age.
	mutable std::vector< const FunctionDerivatives* >
		__current_age_quantities;

    ///The quantities describing the stellar zone.
    const std::vector< EvolvingStellarQuantity* > __evolving_quantities;

	///Forgets any previously calculated quantities for the current age.
	void reset_current_quantities();

public:
    ///Create an evolving stellar zone described by the given quantities.
    ///
    ///WARNING: All quantities are destroyed by the desctructor.
    EvolvingStellarZone(
        ///The quantities describing the zone.
        std::initializer_list<EvolvingStellarQuantity*> evolving_quantities
    ) : 
        __current_age(NaN),
        __current_age_quantities(evolving_quantities.size(), NULL),
        __evolving_quantities(evolving_quantities)
    {}

	///\brief Defines the current orbit, triggering re-calculation of all
	///quantities.
	virtual void configure(
			///The age to set the zone to.
			double age,

			///The angular velocity of the orbit in rad/day.
			double orbital_frequency,

			///The eccentricity of the orbit
			double eccentricity,
			
			///The absolute value of the angular momentum of the orbit.
			double orbital_angmom,

			///The angular (momentum/velocity) of the spin of the zone if
			///the zone is not in a spin-orbit lock (ignored it if is).
			double spin,
			
			///The inclination of the zone relative to the orbit.
			double inclination,
			
			///The argument of periapsis of the orbit in the equatorial
			///planet of the zone.
			double periapsis,
			
			///Is spin an angular velocity instead of angular momentum?
			bool spin_is_frequency
    );

	///\brief The current age value of the given quantity (or its 
	///derivative).
	///
	///Computes the value if necessary or retrieves it from
	//__current_age_quantities if already present.
	double current_age_quantity(
        ///The index of the quantity to return from within the list provided
        ///at construction.
        size_t quantity,

        ///The order of the derivitive to return
        unsigned deriv_order=0
    ) const;

	///\brief The value of the given quantity (or its derivative) at an
    ///arbitrary age.
	///
	///Always computes the value/derivatives.
    double any_age_quantity(
        ///The index of the quantity to return from within the list provided
        ///at construction.
        size_t quantity,

        ///The age at which to evaluate the quantity/derivative.
        double age, 

        ///The order of the derivitive to return
        unsigned deriv_order=0
    ) const;

    ///Return the last age with which ::configure() was called.
    double current_age() {return __current_age;}

	///\brief Change the body as necessary at the given age.
	///
	///Handles things like interpolation discontinuities. 
	void reached_critical_age(double age);

	///\brief The next age when the evolution needs to be stopped for a
	///change in one of the bodies.
	double next_stop_age() const;
};

#endif