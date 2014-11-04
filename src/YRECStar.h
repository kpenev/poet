#ifndef __YREC_STAR_H
#define __YREC_STAR_H

/**\file
 *
 * \brief Declares the class for stars that user the YREC tracks.
 */

#include "SaturatingSkumanichWindBody.h"
#include "ExponentialDecayDiffRotBody.h"
#include "StellarEvolution.h"
#include "YRECCore.h"
#include "YRECEnvelope.h"

class YRECStar : public SaturatingSkumanichWindBody,
				 public ExponentialDecayDiffRotBody {
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

			///The strength of the wind.
			double wind_strength,

			///The frequency at which the wind loss saturates in rad/day.
			double wind_saturation_frequency,

			///The timescale for differential rotation coupling.
			double diff_rot_coupling_timescale,

			///A StellarEvolution interpolator.
			const StellarEvolution &evolution) :
		SaturatingSkumanichWindBody(wind_strength,
									wind_saturation_frequency),
		ExponentialDecayDiffRotBody(diff_rot_coupling_timescale),
		__low_mass(mass<evolution.get_mass_break()),
		__luminosity(evolution.interpolate_luminosity(mass)),
		__lifetime(9*std::pow(mass, -3)),
		__core_formation(evolution.core_formation_age()),
		__envelope(mass, evolution.interpolate_radius(mass),
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
		if(zone_index==0) return __envelope;
		else return __core;
	}

	///The envelope of the star.
	const YRECEnvelope &envelope() const {return __envelope;}

	///The core of the star.
	const YRECCore &core() const
	{
#ifdef DEBUG
		assert(__low_mass);
#endif
		return __core;
	}

	///See DissipatingBody::zone().
	const DissipatingZone &zone(unsigned zone_index) const
	{
#ifdef DEBUG
		assert(zone_index<=1);
#endif
		if(zone_index==0) return __envelope;
		else return __core;
	}

	///The lifetime of the star (where tracks end).
	double lifetime() {return __lifetime;}

	///The luminosity of the star at the given age.
	double luminosity(double age) const {return (*__luminosity)(age);}

	///Cleanup after the star.
	~YRECStar() {delete __luminosity;}

	///Is this a low mass star?
	bool is_low_mass() {return __low_mass;}

	///The age when the core forms.
	double core_formation_age() {return __core_formation;}
};

#endif
