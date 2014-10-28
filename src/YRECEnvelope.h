#ifndef __YREC_ENVELOPE_H
#define __YREC_ENVELOPE_H

/**\file
 *
 * \brief Declares a class representing convective zones in low mass YREC
 * stars or the whole star for high mass stars.
 */

#include "DissipatingZone.h"

///\brief Surface convective zone for low mass YREC stars or the entire star
///for high mass stars.
class YRECEnvelope : public DissipatingZone {
private:
	///The age for the last configure() call.
	double __current_age, __stellar_mass;

	///\brief Identifiers for the various age dependent values which are only
	///computed once per fixed age.
	///
	///All quantities also refer to their derivatives.
	enum CurrentAgeQuantities {
		RADIUS, ///< The radius
		INERTIA,///< The moment of inertia
		NUM_CURRENT_AGE_QUANTITIES///< Number of quantities tracked.
	};

	///\brief Pre-computed values and derivatives for quantities which only
	///depend on age at the current age.
	mutable std::vector< const FunctionDerivatives* >
		__current_age_quantities;

	const EvolvingStellarQuantity
		///\brief The radius of the star in \f$R_\odot\f$ as a function of
		///age in Gyr.
		*__outer_radius,

		///\brief The moment of inertia of the stellar convective envelope in
		///units of \f$M_\odot \cdot R_\odot^2\f$ as a function of age in
		///Gyrs.
		*__moment_of_inertia;

	///Forgets any previously calculated quantities for the current age.
	void reset_current_quantities();

	///Resets the zone to what the default constructor creates.
	void reset();

	///\brief The current age value of the given quantity (or its 
	///derivative).
	///
	///Computes the value if necessary or retrieves it from
	//__current_age_quantities if already present.
	double current_age_quantity(CurrentAgeQuantities quantity,
			unsigned deriv_order=0) const;
public:
	YRECEnvelope(
			///The mass of the star this zone is part of.
			double stellar_mass,

			///The radius of the star.
			const EvolvingStellarQuantity *outer_radius=NULL,

			///The moment of inertia of the zone.
			const EvolvingStellarQuantity *moment_of_inertia=NULL) :
		__current_age(NaN), __stellar_mass(stellar_mass),
		__current_age_quantities(NUM_CURRENT_AGE_QUANTITIES, NULL),
		__outer_radius(outer_radius), __moment_of_inertia(moment_of_inertia)
	{}

	///Set the properties of the zone.
	void create(double stellar_mass,
			const EvolvingStellarQuantity *outer_radius,
			const EvolvingStellarQuantity *moment_of_inertia);

	///\brief Defines the current orbit, triggering re-calculation of all
	///quantities.
	void configure(
			///The age to set the zone to.
			double age,

			///The angular velocity of the orbit in rad/day.
			double orbital_frequency,

			///The eccentricity of the orbit
			double eccentricity,
			
			///The absolute value of the angular momentum of the orbit.
			double orbital_angmom,

			///The angular momentum of the spin of the zone if the zone is
			///not locked (ignored it if is).
			double spin_angmom,
			
			///The inclination of the zone relative to the orbit.
			double inclination,
			
			///The argument of periapsis of the orbit in the equatorial
			///planet of the zone.
			double periapsis);

	///See DissipatingZone::moment_of_inertia.
	double moment_of_inertia(int deriv_order=0) const
	{return current_age_quantity(INERTIA, deriv_order);}

	///See DissipatingZone::outer_radius.
	double outer_radius(int deriv_order=0) const
	{return current_age_quantity(RADIUS, deriv_order);}

	///See DissipatingZone::outer_mass.
	double outer_mass(int deriv_order=0) const
	{return (deriv_order ? 0 : stellar_mass);}

	~YRECEnvelope() {reset();}
};

#endif
