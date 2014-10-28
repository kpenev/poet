#ifndef __YREC_STAR_H
#define __YREC_STAR_H

/**\file
 *
 * \brief Declares the class for stars that user the YREC tracks.
 */

#include "DissipatingBody.h"

class YRECStar : public DissipatingBody {
private:
	///\brief Is this a low mass star according to the stellar evolution with
	///which it was constructed?
	bool __low_mass;

	const EvolvingStellarQuantity
		///\brief The luminosity of the star in \f$L_\odot\f$ as a function
		///of age in Gyr.
		///
		///Since the luminosity does not enter in the equations solved it
		///might not be defined. If that is the case, this member should be
		///NULL.
		*__luminosity;

	///The age at which the star leaves the main sequence in Gyrs.
	double __lifetime,

		   ///\brief The age at which the core first forms in Gyr.
		   __core_formation;

	///The surface zone of the star (the entire star if high mass).
	YRECEnvelope __envelope;

	///The core of the star (NULL if high mass).
	YRECCore __core;

public:
	YRECStar(
			///Mass of the star
			double mass,

			///A StellarEvolution interpolator.
			const StellarEvolution &evolution) :
		__low_mass(mass<evolution.get_mass_break()),
		__luminosity(evolution.interpolate_luminosity(mass)),
		__lifetime(9*std::pow(mass, -3)),
		__core_formation(evolution.core_formation_age())
		__envelope(mass, evolution.interpolate_radius(mass)
				(__low_mass
				 ? evolution.interpolate_moment_of_inertia(mass, convective)
				 : evolution.interpolate_moment_of_inertia(mass, total))),
		__core(__core_formation, 
				(__low_mass
				 ? evolution.interpolate_zone_mass(mass, radiative)
				 : NULL),
				(__low_mass
				 ? evolution.interpolate_core_boundary(mass)
				 : NULL),
				(__low_mass
				 ? evolution.interpolate_moment_of_inertia(mass, radiative)
				 : NULL))
		{}

	///The number of zones the body consists of.
	unsigned number_zones() const {return (__low_mass ? 2 : 1);}

	///See DissipatingBody::zone().
	DissipatingZone &zone(unsigned zone_index)
	{
#ifdef DEBUG
		assert(zone_index<=1);
#endif
		return (zone_index ? __core : __envelope);
	}

	///Cleanup after the star.
	~YRECStar() {delete __luminosity;}
};

#endif
