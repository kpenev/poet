/**\file
 *
 * \brief Defines the StellarSystem class.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __STELLAR_SYSTEM_H
#define __STELLAR_SYSTEM_H

#include "StellarQ.h"
#include "Planet.h"
#include "AstronomicalConstants.h"
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
///\ingroup StellarSystem_group
class StellarSystem {
private:
	///\brief Age in Gyr at which evolution starts when searching for initial
	///conditions
//	const static double MIN_AGE = 0.005;

	///Precision to which orbit should be solved
//	const static double PRECISION=1e-3;

	///\brief What energy to assign if the planet does not survive to the
	///present age when searching for initial conditions.
//	const static double INVALID_ENERGY = 0;

	///The name of the stellar system (e.g. "HAT-P-20")
	std::string name;
	
	///The present age of the stellar system in Gyrs.
	double age;

	///The central star in the system.
	Star *star;

	///The planet in the system.
	Planet *planet;

	///Minimum semimajor axis when solving for initial conditions.
	double min_a0,

		   ///Maximum semimajor axis when solving for initial conditions.
		   max_a0,
		   
		   ///Minimum stellar spin when solving for initial conditions.
		   min_w0,
		   
		   ///Maximum stellar spin when solving for initial conditions.
		   max_w0;

	///Initial semimajor axis found by the initial condition solver.
	double a0,
		   
		   ///Initial stellar spin found by the initial condition solver.
		   w0;
public:
	///Construct a stellar system.
	StellarSystem(
			Star *system_star,///< The star in the system.
			Planet *system_planet,///< The planet in the system.
			double age=NaN,///< The present age of the system in Gyr.
			const std::string &system_name=""///< The name of the system.
			);

	///\brief Transforms planet into Earth.
	///
	///\deprecated It was used to calculate the evolution after the planet
	///died. Now the planet dying is handled by changing the equations.
	void transform_into_earth(); 

	//Make a stellar system from the given star and planet having the
	//given age and name.
	//unused for now
	//StellarSystem(Star* system_star, Planet* system_planet,
	//		OrbitSolver solve_ode);

	///Returns the name of the system.
	const std::string get_name() const {return name;}

	///Returns the present age of the system in Gyr.
	double current_age() const {return age;}

	///Returns the star in the system (immutable).
	const Star &get_star() const {return *star;}

	///Returns the planet in the system (immutable).
	const Planet &get_planet() const {return *planet;}

	///\brief The differential equation for the orbital evolution while the
	///planet is still present and not locked to the stellar spin.
	int orbit_differential_equation(
			///The system age in Gyr.
			double age, 

			///The orbital parameters. Should be:
			/// -# a^6.5 (in units of \f$R_\odot^{6.5}\f$
			/// -# Lconv (in units of \f$M_\odot\cdot R_\odot^2 \cdot
			///    \mathrm{rad}/\mathrm{day}\f$
			/// -# Lrad  (same units as Lconv)
			const double* orbital_parameters,


			///The rate of change of the orbital parameters per Gyr from the
			///differential equation.
			double* orbital_derivatives,
			
			///If assume_sign>0, assumes that the stellar spin period is
			///shorter than the orbital period, regardless of its actual
			///value, if it is <0 assumes it is longer, and if assume_sign is
			///zero actually compares the stellar spin and orbital frequency
			///to determine whether the orbit shrinks or exands.
			int assume_sign=0,

			///The wind saturation state to assume. If it is UNKNOWN,
			///compares the stellar spin frequency with the wind saturation
			///frequency to determine the wind angular momentum loss rate,
			///otherwise, regardless of the saturation frequency, the
			///indicated assumption is used.
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

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

/*	int orbit_diffeq_prescribed_rot(
			double age, const double* orbital_parameters,
			double* orbital_derivatives) const;

	int orbit_jacobian_prescribed_rot(
			double age, const double* orbital_parameters,
			double* param_derivs, double* age_derivs) const;*/

	///\brief The differential equation for the orbital evolution while the
	///planet is still present and locked to the stellar spin.
	///
	
	///The differential equation governing the orbital evolution in the case
	///of the star's rotation locked to the orbit: takes age and dependent
	///variables on input and returns the derivatives of the dependent
	///variables on output.
	///The orbital parameters should be: 
	/// 0. a^6.5 (in units of Rsun^6.5), 
	/// 1. Lconv (in units of Msun*Rsun^2/day, 
	/// 2. Lrad  (same units as Lconv) 
	int locked_orbit_differential_equation(
			///The system age in Gyr.		
			double age,

			///The orbital parameters. Should be:
			/// -# a (in units of \f$R_\odot\f$
			/// -# Lrad (in units of \f$M_\odot\cdot R_\odot^2 \cdot
			///    \mathrm{rad}/\mathrm{day}\f$
			const double* orbital_parameters,

			///The rate of change of the orbital parameters per Gyr from the
			///differential equation.
			double* orbital_derivatives,

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

	///\brief The differential equation governing the evolution of the
	///stellar rotation if no planet or disk is present.
	int no_planet_differential_equation(
			///The system age in Gyr.		
			double age,

			///The orbital parameters. Should be:
			/// -# Lconv (in units of \f$M_\odot\cdot R_\odot^2 \cdot
			///   \mathrm{rad}/\mathrm{day}\f$
			/// -# Lrad  (same units as Lconv)
			const double *orbital_parameters,
			
			///The rate of change of the orbital parameters per Gyr from the
			///differential equation.
			double *orbital_derivatives,

			///See the description of this parameter for the
			///orbit_differential_equation method.
			WindSaturationState assume_wind_saturation=UNKNOWN)
		const;

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

	///\brief The differential equation governing the rotation of the star's
	///radiative zone with the convective zone locked to a disk.
	int locked_conv_differential_equation(
			///The system age in Gyr.		
			double age,

			///A pointer to the radiative core angular momentum in units of
			/// \f$M_\odot\cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			const double *orbital_parameters,
			
			///The rate of change of the radiative core angular momentum per
			///Gyr from the differential equation.
			double *orbital_derivatives);

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
			double *age_derivs) const;

	///\brief The energy definition for the simulated annealing initial
	///condition solver.
//	double annealing_energy();

	///\brief Performs a single step of the simulated annealing initial
	///condition solver.
//	void annealing_step(const gsl_rng *r, double step_size);

	///\brief A definition of the distance between the initial conditions
	///between this system and that.
//	double annealing_metric(StellarSystem* that);

	///\brief Print the current guess for the initial semimajor axis and
	///rotation frequency.
//	void annealing_print();

	///\brief Solve for the initial conditions that reproduce the present
	///conditions using simulated annealing.
/*	bool anneal_solve_IC(
			///Whether to output the progress of solving to stdout.
			bool verbose=true,

			///The minimum initial semimajor axis to consider in AU.
			double min_a0=0.01,

			///The maximum initial semimajor axis to consider in AU.
			double max_a0=0.1,
			
			///The minimum initial stellar rotation to consider in rad/day.
			double min_w0=0.63,
			
			///The maximum initial stellar rotation to consider in rad/day.
			double max_w0=6.3,

			///If an error less than this cannot be achieved it is declared
			///that no initial conditions can reproduced the observed present
			///configuration.
			double max_err=0.1);*/

//	void solve_init(double start_age, double curr_age, double curr_spin);

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
