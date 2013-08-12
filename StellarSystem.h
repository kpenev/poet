#ifndef __STELLAR_SYSTEM_H
#define __STELLAR_SYSTEM_H

#include "Star.h"
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


class StellarSystem {
private:
	//age at which orbit solver starts
	const static double MIN_AGE = 0.005;
	//precision to which orbit should be solved
	const static double PRECISION=1e-3;

	const static double INVALID_ENERGY = 0;

	///The name of the stellar system (e.g. "HAT-P-20")
	std::string name;
	
	///The present age of the stellar system in Gyrs.
	double age;

	///The central star in the system (only one star systems supported at
	///this time).
	Star *star;

	///The planet in the system (only single planet systems supported at
	///this time).
	Planet *planet;

	//assumed bounds on parameters for solving initial conditions
	double min_a0, max_a0, min_w0, max_w0;

	double a0, w0; //initial conditions
public:
	StellarSystem(
			Star *system_star,///< The star in the system
			Planet *system_planet,///< The planet in the system
			double age=NaN,///< The present age of the system
			const std::string &system_name="");///< The name of the system

	void transform_into_earth(); //transforms planet into Earth

	///Make a stellar system from the given star and planet having the
	///given age and name.
	//unused for now
	//StellarSystem(Star* system_star, Planet* system_planet,
	//		OrbitSolver solve_ode);

	///Returns the name of the system.
	const std::string get_name() const {return name;}

	///Returns the present age of the system.
	double current_age() const {return age;}

	///Returns the star in the system (immutable).
	const Star &get_star() const {return *star;}

	///Returns the planet in the system (immutable).
	const Planet &get_planet() const {return *planet;}

	///The differential equation governing the orbital evolution: takes
	///age and dependent variables on input and returns the derivatives
	///of the dependent variables on output.
	///The orbital parameters should be: 
	/// 0. a^6.5 (in units of Rsun^6.5), 
	/// 1. Lconv (in units of Msun*Rsun^2/day, 
	/// 2. Lrad  (same units as Lconv) 
	int orbit_differential_equation(double age, 
			const double* orbital_parameters,
			double* orbital_derivatives, int assume_sign=0,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///The jacobian of the differential equation governing the orbital 
	///evolution: takes age and dependent variables on input and returns
	///the jacobian of the derivatives of the dependent variables on output.
	///The orbital parameters should be: 
	/// * a^6.5 (in units of Rsun^6.5), 
	/// * Lconv (in units of Msun*Rsun^2/day), 
	/// * Lrad  (same units as Lconv) 
	int orbit_jacobian(double age, const double* orbital_parameters,
			double* param_derivs, double* age_derivs, int assume_sign=0,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

/*	int orbit_diffeq_prescribed_rot(
			double age, const double* orbital_parameters,
			double* orbital_derivatives) const;

	int orbit_jacobian_prescribed_rot(
			double age, const double* orbital_parameters,
			double* param_derivs, double* age_derivs) const;*/

	///The differential equation governing the orbital evolution in the case
	///of the star's rotation locked to the orbit: takes age and dependent
	///variables on input and returns the derivatives of the dependent
	///variables on output.
	///The orbital parameters should be: 
	/// 0. a^6.5 (in units of Rsun^6.5), 
	/// 1. Lconv (in units of Msun*Rsun^2/day, 
	/// 2. Lrad  (same units as Lconv) 
	int locked_orbit_differential_equation(double age,
			const double* orbital_parameters,
			double* orbital_derivatives,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///The jacobian of the differential equation governing the orbital
	///evolution in the case of the star's rotation locked to the orbit:
	///takes age and dependent variables on input and returns the jacobian
	///and partial time derivatives of the evolution equations. 
	///The orbital parameters should be: 
	/// 0. a^6.5 (in units of Rsun^6.5), 
	/// 1. Lconv (in units of Msun*Rsun^2/day, 
	/// 2. Lrad  (same units as Lconv) 
	int locked_orbit_jacobian(double age, const double* orbital_parameters,
			double* param_derivs, double* age_derivs,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///The differential equation governing the evolution of the stellar
	///rotation if no planet or disk is present. The orbital parameters
	///should be:
	/// 0. Lconv (in units of Msun*Rsun^2/day)
	/// 1. Lrad (same units as Lconv)
	int no_planet_differential_equation(double age,
			const double *orbital_parameters, double *orbital_derivatives,
			WindSaturationState assume_wind_saturation=UNKNOWN)
		const;

	///The jacobian of the differential equation governing the evolution of
	///the stellar rotation if no planet or disk is present. The orbital
	///parameters should be:
	/// 0. Lconv (in units of Msun*Rsun^2/day)
	/// 1. Lrad (same units as Lconv)
	int no_planet_jacobian(double age, const double *orbital_parameters,
			double *param_derivs, double *age_derivs,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///The differential equation governing the rotation of the star's
	///radiative zone with the convective zone locked to a disk: takes age
	///and pointer to the radiative antgular momentum (in units of
	///Msun*Rsun^2/day) dependent variables on input and updates the
	///derivative of Lrad on output.
	int locked_conv_differential_equation(double age,
			const double *orbital_parameters, double *orbital_derivatives);

	///The jacobian of the differential equation governing the rotation of
	///the star's radiative zone with the convective zone locked to a disk:
	///takes age and pointer to the radiative antgular momentum (in units of
	///Msun*Rsun^2/day) dependent variables on input and updates the jacobian
	///and partial time derivative of the Lrad equation on output.
	int locked_conv_jacobian(double age, const double* orbital_parameters,
			double* param_derivs, double *age_derivs) const;


	///The energy definition for the simulated annealing initial condition
	///solver
	double annealing_energy();

	///Performs a single step of the simulated annealing initial condition
	///solver
	void annealing_step(const gsl_rng *r, double step_size);

	///A definition of the distance between the initial conditions between
	///this system and that.
	double annealing_metric(StellarSystem* that);

	///Print the current guess for the initial semimajor axis and rotation
	///frequency.
	void annealing_print();

	///Solve for the initial conditions that reproduce the present conditions
	///using simulated annealing.
	bool anneal_solve_IC(bool verbose=true, double min_a0=0.01,
			double max_a0=0.1, double min_w0=0.63, double max_w0=6.3,
			double max_err=0.1);

//	void solve_init(double start_age, double curr_age, double curr_spin);

	///Writes the evolution of the system (orbital and stellar) to a file
	///with the given name.
	void output_evolution(std::string &filename) const;
};

bool stop_evolution(double age, const double* y, void *stellar_system);

#endif
