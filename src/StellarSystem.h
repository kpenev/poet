/**\file
 *
 * \brief Defines the StellarSystem class.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __STELLAR_SYSTEM_H
#define __STELLAR_SYSTEM_H

#include "StellarDissipation.h"
#include "Planet.h"
#include "AstronomicalConstants.h"
#include "TidalDissipation.h"
#include "Common.h"
#include "Error.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_siman.h>
#include <string>
#include <limits>
#include <iostream>


///\brief Describes a planet-star system.
///
///The following variable are referred to throughout:
/// - age: Gyr
/// - a: semimajor axis in \f$R_\odot\f$
/// - inclination: angle between surface spin and orbital angular momentum in
///   radians
/// - Lconv: Convective zone angular momentum in
///   \f$M_\odot\cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
/// - Lrad_parallel: Radiative zone angular momentum component along Lconv in
///   the same units as Lconv
/// - Lrad_perpendicular: Radiative zone angular momentum component along a
///   direction rotated from Lconv by 90 degrees counter-clockwise in the 
///   same units as Lconv.
///
///All rates of change are per Gyr
///
///\ingroup StellarSystem_group
class StellarSystem {
private:
	///The name of the stellar system (e.g. "HAT-P-20")
	std::string __name;
	
	///The present age of the stellar system in Gyrs.
	double __age;

	///The central star in the system.
	Star &__star;

	///The planet in the system.
	const Planet &__planet;

	///\brief The differential equation governing the rotation of the star's
	///radiative zone with the convective zone locked to a disk.
	int locked_conv_differential_equations(
			///A pointer to Lrad_parallel (Lrad_perpendicular assumed 0).
			const double *parameters,
			
			///On output is set to the rate of change of Lrad.
			double *differential_equations);

	///\brief The differential equation governing the evolution of the
	///stellar rotation if no planet or disk is present.
	int no_planet_differential_equations(
			///Should contain in the following order:
			/// -# Lconv
			/// -# Lrad_parallel
			/// -# Lrad_perpendicular
			const double *parameters,

			///The saturation state to assume for the stellar wind.
			WindSaturationState assume_wind_saturation,
			
			///On outputs is set to the rate of change of the orbital
			///parameters.
			double *evolution_rates)
		const;

	///\brief The differential equation for the evolution while the
	///planet is still present and locked to the stellar spin.
	int locked_to_orbit_differential_equations(
			///Should contain in the following order:
			/// -# a
			/// -# inclination
			/// -# Lrad_parallel
			/// -# Lrad_perpendicular
			const double *parameters,

			///The saturation state to assume for the stellar wind.
			WindSaturationState assume_wind_saturation,

			///See description in differential_equations method
			const SpinOrbitLockInfo &star_lock,

			///The evolution rates of various quantities due to tidal
			///dissipation.
			const TidalDissipation &dissipation,

			///On outputs is set to the rate of change of the orbital
			///parameters.
			double *evolution_rates) const;

	///\brief The differential equation for the evolution while the
	///planet is still present and not locked to the stellar spin.
	int not_locked_differential_equations(
			///The orbital parameters. Should be:
			/// -# a^6.5
			/// -# inclination
			/// -# Lconv
			/// -# Lrad_parallel
			/// -# Lrad_perpendicular
			const double *parameters,

			///The wind saturation state to assume. If it is UNKNOWN,
			///compares the stellar spin frequency with the wind saturation
			///frequency to determine the wind angular momentum loss rate,
			///otherwise, regardless of the saturation frequency, the
			///indicated assumption is used.
			WindSaturationState assume_wind_saturation,

			///The evolution rates of various quantities due to tidal
			///dissipation.
			const TidalDissipation &dissipation,

			///The rate of change of the orbital parameters per Gyr from the
			///differential equation.
			double *evolution_rates) const;

public:
	///Construct a stellar system.
	StellarSystem(
			Star &star,///< The star in the system.
			const Planet &planet,///< The planet in the system.
			double age=NaN,///< The present age of the system in Gyr.
			const std::string &system_name=""///< The name of the system.
			);

	///Returns the name of the system.
	const std::string get_name() const {return __name;}

	///Returns the present age of the system in Gyr.
	double current_age() const {return __age;}

	///Returns the star in the system (immutable).
	const Star &get_star() const {return __star;}

	///Returns the star in the system (immutable).
	Star &get_star() {return __star;}

	///Returns the planet in the system (immutable).
	const Planet &get_planet() const {return __planet;}

	///\brief The differential equation and jacobian for the evolution of the
	///system.
	int differential_equations(
			///The system age in Gyr.
			double age, 

			///Contains the variables being evolved. The expected content
			///depends on the values of evolution_mode and star_lock.
			///
			///If evolution_mode is LOCKED_TO_DISK (Lconv and Lrad are 
			///assumed aligned):
			/// -# Lrad_parallel
			///
			///If evolution mode is NO_PLANET:
			/// -# Lconv
			/// -# Lrad_parallel
			/// -# Lrad_perpendicular
			///
			///The TABULATION evolution mode is not allowed and for BINARY
			///the content, depending on star_lock is:
			///
			///If star_lock converts to true (i.e. some lock is held):
			/// -# a
			/// -# inclination
			/// -# Lrad_parallel
			/// -# Lrad_perpendicular
			///
			///If star is not locked
			/// -# a^6.5
			/// -# inclination
			/// -# Lconv
			/// -# Lrad_parallel
			/// -# Lrad_perpendicular
			const double *parameters,

			///The evolution mode to assume for the system.
			EvolModeType evolution_mode, 

			///The saturation state to assume for the stellar wind.
			WindSaturationState assume_wind_saturation,

			///Whether a harmonic of the stellar rotation frequency is locked
			///to a harmonic of the orbital period, and which harmonics.
			const SpinOrbitLockInfo &star_lock,

			///The pre-computed tidal dissipation evolution rates.
			const TidalDissipation &dissipation,

			///On outputs gets filled  with the rates at which the entries in
			///parameters evolve. It is assumed that sufficient space has
			///been allocated to hold the results.
			double *differential_equations);

/*
	///\brief The jacobian of the differential equation governing the orbital 
	///evolution while the planet is still present.
	int orbit_jacobian(
			///The system age in Gyr.
			double age,
			
			///See the description of this parameter for the
			///orbit_differential_equation method.
			const double* orbital_parameters,

			///The output jacobian. All parameters are in the same units as
			///in the orbital_parameters argument.
			double* param_derivs,
			
			///The output partial age derivatives of the orbital differential
			///equations. Orbital parameters are in the same units as in the
			///orbital_parameters method and age is in Gyr.
			double* age_derivs,
			
			///See the description of this parameter for the
			///orbit_differential_equation method.
			int assume_sign=0,

			///See the description of this parameter for the
			///orbit_differential_equation method.
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///\brief The jacobian of the differential equation governing the orbital
	///evolution in the case of the star's rotation locked to the orbit.
	int locked_orbit_jacobian(\
			///The system age in Gyr.
			double age,
			
			///See the description of this parameter for the
			///locked_orbit_differential_equation method.
			const double* orbital_parameters,

			///The output jacobian. All parameters are in the same units as
			///in the orbital_parameters argument.
			double* param_derivs,
			
			///The output partial age derivatives of the orbital differential
			///equations. Orbital parameters are in the same units as in the
			///orbital_parameters method and age is in Gyr.
			double* age_derivs,

			///See the description of this parameter for the
			///orbit_differential_equation method.
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///\brief The jacobian of the differential equation governing the
	///evolution of the stellar rotation if no planet or disk is present.
	int no_planet_jacobian(
			///The system age in Gyr.		
			double age,
			
			///See the description of this parameter for the
			///no_planet_differential_equation method.
			const double *orbital_parameters,

			///The output jacobian. All parameters are in the same units as
			///in the orbital_parameters argument.
			double *param_derivs,
			
			///The output partial age derivatives of the orbital differential
			///equations. Orbital parameters are in the same units as in the
			///orbital_parameters method and age is in Gyr.
			double *age_derivs,

			///See the description of this parameter for the
			///orbit_differential_equation method.
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///\brief The jacobian of the differential equation governing the
	///rotation of the star's radiative zone with the convective zone locked
	///to a disk.
	int locked_conv_jacobian(
			///The system age in Gyr.		
			double age,
			
			///See the description of this parameter for the
			///locked_conv_differential_equation method.
			const double* orbital_parameters,

			///The output jacobian. All parameters are in the same units as
			///in the orbital_parameters argument.
			double* param_derivs,
			
			///The output partial age derivatives of the orbital differential
			///equations. Orbital parameters are in the same units as in the
			///orbital_parameters method and age is in Gyr.
			double *age_derivs) const;*/

	///\brief Writes the evolution of the system (orbital and stellar) to a
	///file with the given name.
	void output_evolution(std::string &filename) const;
};

///Not sure what it does but is not used.
bool stop_evolution(double age, const double* y, void *stellar_system);

///\brief The rate of change of the angular velocity of the stellar 
///convective zone that would occur if no planet were present in rad/day/Gyr.
double no_planet_dwconv_dt(
		///The age of the system in Gyr.
		double age,
		
		///The present orbit for a spin-orbit locked system.
		const std::valarray<double> &orbit,

		///The stellar system
		const StellarSystem &system);

#endif
