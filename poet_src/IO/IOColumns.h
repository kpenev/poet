/**\file
 *
 * \brief Declares enumerations that define the possible input/outut
 * quantities.
 */

#ifndef __IO_COLUMNS_H
#define __IO_COLUMNS_H

#include <vector>
#include <string>

///Isolates the tags for the input columns.
namespace InCol {
	///\brief Tags for the quantities required to fully specify the system to
	///evolve.
	enum InputColumns {
		///\brief Wind strength in 
		/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
		/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
		WINDK=0, 

		LOW_MASS_WINDK=WINDK,///< The wind strength for low mass stars
		HIGH_MASS_WINDK,///< The wind strength for high mass stars
		WIND_SAT_W, ///< Wind saturation frequency in rad/day.

		///Wind saturation frequency in rad/day for low mass stars.
		LOW_MASS_WIND_SAT_W=WIND_SAT_W,

		///Wind saturation frequency in rad/day for high mass stars.
		HIGH_MASS_WIND_SAT_W,

		WIND_SAT_P, ///< Wind saturation period in days.

		///Wind saturation period in days for low mass stars.
		LOW_MASS_WIND_SAT_P=WIND_SAT_P, 

		///Wind saturation period in days for high mass stars.
		HIGH_MASS_WIND_SAT_P,

		///\brief Core-envelope coupling timescale in Myr, ignore for high
		///mass stars.
		CORE_ENV_COUPLING_TIMESCALE,

#ifdef TWO_QS
		LGQ,///< Lg(Q*) outside the inertial mode range.
		LGQ_INERTIAL,///< Lg(Q*) in the inertial mode range.
#else
		LGQ,///< Lg(Q*) at tidal frequency of 1 day.
		LGQ_POWERLAW,///< powerlaw index of Q* dependence.
		LGQ_MIN,///< Lg(Q*) outside the inertial mode range.
#endif
		MSTAR, ///< Mass of the star in \f$M_\odot\f$.
		MPLANET, ///< Mass of the planet in Jupiter masses.
		RPLANET, ///< Radius of the planet in Jupiter radii.
		PLANET_FORMATION_AGE,///< Age when planet appears in Gyr.
		WDISK,///< Stellar surface spin while disk is present in rad/day.

		///Stellar surface rotation period in days while disk is present.
		PDISK,
		
		TDISK, ///< Age in Myr when disk dissipates.
		A_FORMATION, ///< Semimajor axis at which the planet forms in AU.
		P_FORMATION, ///< Orbital period in days at which the planet forms.
		E_FORMATION, ///< Eccentricity at which the planet forms.
		INCLINATION_FORMATION, ///< Inclination with which the planet forms.
		TSTART, ///< The minimum age to start evolution at in Gyr.
		TEND, ///< The maximum age to stop the evolution at in Gyr.

		///\brief Initial rotation of the radiative core in rad/day if the
		///evolution starts after the core has formed.
		START_WRAD, 

		///\brief Surface rotation of the star in rad/day if the evolution
		///starts after the disk has dissipated.
		START_WSURF, 

		MAX_STEP,///< The maximum timestep to take.

		///The number of significant figures require of the evolution.
		PRECISION,

		OUT_FNAME,///< The name of the file to write the evolution to.

		///The number of real values quantities.
		NUM_REAL_INPUT_QUANTITIES=OUT_FNAME,

		///\brief Should the evolution start with the stellar surface
		///spinning synchronously with the orbit? 
		START_LOCKED,	

		///A list of ages guaranteed to be included in the tabulated orbit.
		REQUIRED_AGES,

		SKIP,///< A column which is not needed to calculate the evolution.

		///The number of different input quantities supported.
		NUM_INPUT_QUANTITIES=SKIP
	};
};

///Isolates the tags for the output columns.
namespace OutCol {
	///Tags for the possible columns to output.
	enum OutputColumns {
		AGE,///< Age of the system in Gyr.

		///\brief Moment of inertia of the convective zone of the star (low
		///mass stars only) in \f$M_\odot R_\odot^2\f$.
		ICONV,

		///\brief Moment of inertia of the radiative zone of the star (low
		///mass stars only) in \f$M_\odot R_\odot^2\f$.
		IRAD,

		///Total moment of inertia of the star in \f$M_\odot R_\odot^2\f$.
		ITOT,

		RSTAR,///< Radius of the star in \f$R_\odot\f$.
		LSTAR,///< Luminosity of the star in \f$L_\odot\f$.

		///Radius of the stellar core in \f$R_\odot\f$ (low mass stars only).
		RRAD,

		///Mass of the stellar core in \f$M_\odot\f$ (low mass stars only).
		MRAD,

		///\brief Age derivative of the convective moment of inertia in
		/// \f$M_\odot R_\odot^2/Gyr\f$.
		ICONV_DERIV,

		///\brief Age derivative of the radiative moment of inertia in
		/// \f$M_\odot R_\odot^2/Gyr\f$.
		IRAD_DERIV,

		///\brief Age derivative of the moment of inertia of the entire star
		///in \f$M_\odot R_\odot^2/Gyr\f$.
		ITOT_DERIV,

		///Age derivative of the stellar radius in \f$R_\odot/Gyr\f$.
		RSTAR_DERIV,

		///Age derivative of the radius of the radiative core in
		/// \f$R_\odot/Gyr\f$.
		RRAD_DERIV,

		///\brief Age derivative of the mass of the radiative core in
		/// \f$M_\odot/Gyr\f$.
		MRAD_DERIV,

		///Second derivative of the moment of ivertia of the convective zone.
		ICONV_SECOND_DERIV,

		///Second derivative of the moment of ivertia of the radiative zone.
		IRAD_SECOND_DERIV,

		///Second derivative of the moment of inertia of the entire star.
		ITOT_SECOND_DERIV,

		///Second derivative of the core-envelope boundary.
		RRAD_SECOND_DERIV,

		///The index of the last quantity requiring no orbital evolution.
		LAST_NO_ORBIT=RRAD_SECOND_DERIV,

		SEMIMAJOR,///< Semimajor axis of the orbit in AU.
		ECCENTRICITY,///< The eccentricity of the orbit.
		WORB,///< The orbital frequency in rad/day.
		PORB,///< The orbital period days.
		LORB,//< The orbital angular momentum

		///\brief The angle between the stellar surface spin and orbital
		///angular momentum in radians
		CONV_INCLINATION,

		///\brief The angle between the stellar core spin and orbital angular
		///momentum in radians
		RAD_INCLINATION,

		///\brief The orbital periapsis in the reference frame of the stellar
		///convective zone in radians
		CONV_PERIAPSIS,

		///\brief The orbital periapsis in the reference frame of the stellar
		///radiative zone in radians
		RAD_PERIAPSIS,

		///\brief Angular momentum of the convective zone of the star in
		/// \f$ M_\odot R_\odot^2 \mathrm{rad}/\mathrm{day}\f$ (low mass
		///stars only)
		LCONV,

		///\brief Angular momentum of the radiative zone of the star in
		/// \f$ M_\odot R_\odot^2 \mathrm{rad}/\mathrm{day}\f$ (low mass
		///stars only)
		LRAD,

		///\brief Angular momentum of the entire star in
		/// \f$ M_\odot R_\odot^2 \mathrm{rad}/\mathrm{day}\f$
		LTOT,

		WSURF, ///< Angular velocity of the stellar surface in rad/day.

		///\brief Angular velocity of the stellar core in rad/day (low mass
		///stars only).
		WRAD,

		PSURF, ///< Spin period of the stellar surface in days.

		///Spin period of the stellar core in days (low mass stars only).
		PRAD,

		///The number of real valued output quantities.
		NUM_REAL_OUTPUT_QUANTITIES,

		///The evolution mode for the step that starts at this age.
		EVOL_MODE=NUM_REAL_OUTPUT_QUANTITIES,

		WIND_STATE,///< The saturation state of the wind.

		///The number of different output quantities supported.
		NUM_OUTPUT_QUANTITIES
	};
};

///Define the names of the output columns.
const std::vector<std::string> &output_column_names();

#endif
