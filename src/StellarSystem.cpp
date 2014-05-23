/**\file
 *
 * \brief The definition of some of the methods of the StellarSystem class.
 *
 * \ingroup StellarSystem_group
 */

#include "StellarSystem.h"
#include "Error.h"
#include "OrbitSolver.h"
#include <limits>
#include <iostream>

#include <gsl/gsl_matrix.h>

int StellarSystem::locked_conv_differential_equations(
		const double *parameters, double *evolution_rates)
{
	evolution_rates[0]=-__star.differential_rotation_torque_angmom(
			__star.disk_lock_frequency()*
			__star.moment_of_inertia(convective), parameters[0]).real();
	return GSL_SUCCESS;
}

int StellarSystem::no_planet_differential_equations(const double *parameters,
		WindSaturationState assume_wind_saturation, double *evolution_rates)
	const
{
	double Lconv=parameters[0], 
		   wconv=__star.spin_frequency(convective, Lconv);
	std::complex<double> Lrad(parameters[1], parameters[2]),
		coupling_torque=__star.differential_rotation_torque_angmom(Lconv,
				Lrad);
	double rotation=coupling_torque.imag()/Lconv;
	evolution_rates[0]=coupling_torque.real() -
		__star.wind_torque(wconv, assume_wind_saturation);
	evolution_rates[1]=-coupling_torque.real()+rotation*Lrad.imag();
	evolution_rates[2]=-coupling_torque.imag()-rotation*Lrad.real();
	return GSL_SUCCESS;
}

int StellarSystem::locked_to_orbit_differential_equations(
		const double *parameters, WindSaturationState assume_wind_saturation,
		const SpinOrbitLockInfo &star_lock,
		const TidalDissipation &dissipation, double *evolution_rates) const
{
#ifdef DEBUG
	assert(parameters[0]>0);
#endif

	std::complex<double> Lrad(parameters[2], parameters[3]);
	double wconv=dissipation.orbital_frequency(),
		   Iconv=__star.moment_of_inertia(convective),
		   Lconv=wconv*Iconv,
		   dIconv_dt=__star.moment_of_inertia_deriv(convective),
		   wind_torque=__star.wind_torque(wconv, assume_wind_saturation);
	std::complex<double>
		coupling_torque=__star.differential_rotation_torque_angmom(Lconv,
				Lrad);
	throw Error::NotImplemented("Locked orbital evolution");

/*
		   semimajor_rsun=orbital_parameters[0],
		   semimajor_meters=semimajor_rsun*AstroConst::solar_radius,
		   semimajor_au=semimajor_rsun*Rsun_AU,
		   mstar=__star.mass()*AstroConst::solar_mass,
		   mplanet=__planet.mass()*AstroConst::jupiter_mass,
		   mtotal=mstar+mplanet,
		   torque_coef=std::sqrt(semimajor_meters/mtotal/AstroConst::G)*
			   AstroConst::solar_radius/AstroConst::day;
	orbital_derivatives[0]=2.0*(torque_coef*(coupling_torque-wind_torque)
		- dIconv_dt/semimajor_rsun)/
		(mstar*mplanet/(mtotal*AstroConst::solar_mass) -
		 3.0*Iconv/std::pow(semimajor_rsun,2));
	orbital_derivatives[1]=-coupling_torque;*/
	return GSL_SUCCESS;
}

int StellarSystem::not_locked_differential_equations(
		const double *parameters,
		WindSaturationState assume_wind_saturation,
		const TidalDissipation &dissipation,
		double *evolution_rates) const
{
	double da_dt=-dissipation(0, Dissipation::SEMIMAJOR_DECAY)*Rsun_AU;
	evolution_rates[0]=6.5*std::pow(parameters[0], 11.0/13.0)*da_dt;
	evolution_rates[1]=-dissipation(0, Dissipation::INCLINATION_DECAY);
	no_planet_differential_equations(parameters, assume_wind_saturation,
		evolution_rates+2);
	evolution_rates[2]+=dissipation(0, Dissipation::TORQUEZ);
	return GSL_SUCCESS;
}

StellarSystem::StellarSystem(Star &star, const Planet &planet,
		double age, const std::string &system_name) :
	__name(system_name), __age(age), __star(star), __planet(planet)
{
}

int StellarSystem::differential_equations(
		double age, const double *parameters,
		EvolModeType evolution_mode, 
		WindSaturationState assume_wind_saturation,
		const SpinOrbitLockInfo &star_lock,
		const TidalDissipation &dissipation,
		double *evolution_rates)
{
	double Lconv, inclination;
	switch(evolution_mode) {
		case LOCKED_TO_DISK :
			Lconv=__star.disk_lock_frequency()*
				__star.moment_of_inertia(age, convective);
			inclination=0;
			break;
		case NO_PLANET : Lconv=parameters[0]; inclination=0; break;
		case BINARY :
			Lconv=(star_lock ? dissipation.orbital_frequency()*
					__star.moment_of_inertia(age, convective) :
					parameters[2]);
			inclination=parameters[1];
			break;
		default :
#ifdef DEBUG
			assert(false)
#endif
				;
	}
	if(__star.age()!=age) __star(age, Lconv, inclination);
	if(evolution_mode==LOCKED_TO_DISK)
		return locked_conv_differential_equations(parameters,
				evolution_rates);
	else if(evolution_mode==NO_PLANET)
		return no_planet_differential_equations(parameters, 
				assume_wind_saturation, evolution_rates);
	else if(evolution_mode==BINARY) {
		if(star_lock) 
			return locked_to_orbit_differential_equations(parameters,
					assume_wind_saturation, star_lock, dissipation,
					evolution_rates);
		else return not_locked_differential_equations(parameters,
				assume_wind_saturation, dissipation, evolution_rates);
	} else
		throw Error::BadFunctionArguments("Unrecognized evolution mode in "
				"StellarSystem::differential_equations");
}

/*
int StellarSystem::no_planet_jacobian(double age,
		const double *orbital_parameters, double *param_derivs,
		double *age_derivs,
		WindSaturationState assume_wind_saturation) const
{
	double Lconv=orbital_parameters[0],
		   Lrad=orbital_parameters[1],
		   wconv=__star.spin_frequency(age, convective, Lconv),
		   dwconv_dLconv=__star.spin_frequency_angmom_deriv(age, convective,
				Lconv),
		   dcoupling_torque_dLconv=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective),
		   dcoupling_torque_dLrad=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, radiative),
		   dcoupling_torque_dage=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad);
	param_derivs[2]=-dcoupling_torque_dLconv;
	param_derivs[3]=-dcoupling_torque_dLrad;
	param_derivs[0]=-param_derivs[2]
		-__star.wind_torque_freq_deriv(age, wconv, assume_wind_saturation)*
		dwconv_dLconv;
	param_derivs[1]=-param_derivs[3];
	age_derivs[0]=dcoupling_torque_dage -
		__star.wind_torque_age_deriv(age, Lconv, true,
				assume_wind_saturation);
	age_derivs[1]=-dcoupling_torque_dage;
	return GSL_SUCCESS;
}

///The jacobian of the differential equation governing the rotation of
///the star's radiative zone with the convective zone locked to a disk:
///takes age and pointer to the radiative antgular momentum (in units of
///Msun*Rsun^2/day) dependent variables on input and updates the jacobian
///and partial time derivative of the Lrad equation on output.
int StellarSystem::locked_conv_jacobian(double age, const double* orbital_parameters,
		double* param_derivs, double *age_derivs) const
{
	double wconv=__star.disk_lock_frequency(),
		   Lconv=wconv*__star.moment_of_inertia(age, convective),
		   Lrad=orbital_parameters[0];
	param_derivs[0]=-__star.differential_rotation_torque_deriv(
			age, Lconv, Lrad, radiative),
	age_derivs[0] = -__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective)*wconv*
		__star.moment_of_inertia_deriv(age, convective) -
		__star.differential_rotation_torque_deriv(age, Lconv, Lrad);
	return GSL_SUCCESS;
}

int StellarSystem::orbit_jacobian(double age,
		const double* orbital_parameters, double* param_derivs, 
		double* age_derivs, int assume_sign,
		WindSaturationState assume_wind_saturation) const
{
	double semimajor=std::pow(std::abs(orbital_parameters[0]), 1.0/6.5)*
		Rsun_AU,
		   Lconv=orbital_parameters[1],
		   Lrad=orbital_parameters[2],
		   wconv=__star.spin_frequency(age, convective, Lconv),
		   dwconv_dLconv=__star.spin_frequency_angmom_deriv(age, convective,
				Lconv),
		   dcoupling_torque_dLconv=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective),
		   dcoupling_torque_dLrad=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, radiative),
		   dcoupling_torque_dage=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad),
		   wconv_tidal=(assume_sign==0 ? wconv :
				   (assume_sign>0 ? Inf : -Inf)),
		   semi_deriv = __planet.tidal_decay(age, semimajor, wconv_tidal,
				   false);
	param_derivs[0]= __planet.tidal_decay_semimajor_deriv(age, semimajor,
			wconv_tidal, true);
	param_derivs[2]=0.0;
	param_derivs[1]=__planet.tidal_decay_star_spin_deriv(age, semimajor,
			wconv_tidal, true)*-dwconv_dLconv/std::pow(Rsun_AU, 6.5);

	param_derivs[6]=0;
	param_derivs[7]=-dcoupling_torque_dLconv;
	param_derivs[8]=-dcoupling_torque_dLrad;

	param_derivs[3]=-__planet.orbital_angmom_deriv_semimajor_deriv(semimajor,
			semi_deriv, wconv_tidal)/(6.5*std::pow(semimajor, 5.5))*
			std::pow(Rsun_AU, 6.5);
	param_derivs[4]=-param_derivs[7]
		-__star.wind_torque_freq_deriv(age, wconv, assume_wind_saturation)*
		dwconv_dLconv +
		__planet.tidal_torque_star_spin_deriv(age, semimajor, semi_deriv,
				wconv_tidal)*
		-dwconv_dLconv;
	param_derivs[5]=-param_derivs[8];

	age_derivs[0]=__planet.tidal_decay_age_deriv(age, semimajor, wconv_tidal,
			true)/pow(Rsun_AU, 6.5);
	age_derivs[2]= -dcoupling_torque_dage;
	age_derivs[1]=-age_derivs[2] -
		__star.wind_torque_age_deriv(age, Lconv, true, assume_wind_saturation)
			-__planet.orbit_angmom_deriv_age_deriv(age, semimajor, semi_deriv,
					wconv_tidal);
	return GSL_SUCCESS;
}

int StellarSystem::locked_orbit_jacobian(double age,
		const double* orbital_parameters, double* param_derivs,
		double* age_derivs,
		WindSaturationState assume_wind_saturation) const
{
	double semimajor_rsun=orbital_parameters[0],
		   semimajor_meters=semimajor_rsun*AstroConst::solar_radius,
		   semimajor_au=semimajor_meters/AstroConst::AU,
		   wconv=__planet.orbital_angular_velocity_semimajor(semimajor_au),
		   Iconv=__star.moment_of_inertia(age, convective),
		   dIconv_dt=__star.moment_of_inertia_deriv(age, convective),
		   d2Iconv_dt2=__star.moment_of_inertia_deriv(age, convective, 2),
		   Lconv=wconv*Iconv,
		   Lrad=orbital_parameters[1],
		   dwconv_da=__planet.orbital_angular_velocity_semimajor_deriv(
				   semimajor_au)*Rsun_AU,
		   dLconv_da=dwconv_da*Iconv,
		   dLconv_dt=wconv*dIconv_dt,
		   coupling_torque=__star.differential_rotation_torque_angmom(age,
				   Lconv, Lrad),
		   dcoupling_torque_dLconv=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective),
		   dcoupling_torque_dLrad=__star.differential_rotation_torque_deriv(
				   age, Lconv, Lrad, radiative),
		   dcoupling_torque_da=dcoupling_torque_dLconv*dLconv_da,
		   dcoupling_torque_dt=dcoupling_torque_dLconv*dLconv_dt +
			   __star.differential_rotation_torque_deriv(age, Lconv, Lrad),
		   wind_torque=__star.wind_torque(age, wconv, assume_wind_saturation),
		   dwind_torque_da=__star.wind_torque_freq_deriv(age, wconv,
				   assume_wind_saturation)*dwconv_da,
		   dwind_torque_dt=__star.wind_torque_age_deriv(age, Lconv, false,
				   assume_wind_saturation),
		   total_torque=(coupling_torque-wind_torque),
		   mstar=__star.mass()*AstroConst::solar_mass,
		   mplanet=__planet.mass()*AstroConst::jupiter_mass,
		   mtotal=mstar+mplanet,
		   torque_coef=std::sqrt(semimajor_meters/mtotal/AstroConst::G)/
			   AstroConst::day*AstroConst::solar_radius,
		   dtorque_coef_da=0.5*std::pow(AstroConst::solar_radius,2)/
			   AstroConst::day/
			   std::sqrt(semimajor_meters*mtotal*AstroConst::G),
		   denominator=(mstar*mplanet/(mtotal*AstroConst::solar_mass) -
				   3.0*Iconv/std::pow(semimajor_rsun,2));

	param_derivs[0]=2.0*(dtorque_coef_da*total_torque +
			torque_coef*(dcoupling_torque_da-dwind_torque_da) + 
			dIconv_dt/std::pow(semimajor_rsun, 2))/denominator - 
		12.0*(torque_coef*total_torque - dIconv_dt/semimajor_rsun)/
		std::pow(denominator, 2)*Iconv/std::pow(semimajor_rsun, 3);

	param_derivs[1]=2.0*(torque_coef*dcoupling_torque_dLrad)/denominator;

	param_derivs[2]=-dcoupling_torque_da;
	param_derivs[3]=-dcoupling_torque_dLrad;

	age_derivs[0]=2.0*(torque_coef*(dcoupling_torque_dt-dwind_torque_dt) - 
			d2Iconv_dt2/semimajor_rsun)/denominator + 
		6.0*(torque_coef*total_torque - dIconv_dt/semimajor_rsun)/
		std::pow(denominator, 2)*dIconv_dt/std::pow(semimajor_rsun, 2);

	age_derivs[1]=-dcoupling_torque_dt;

	return GSL_SUCCESS;
}*/

void StellarSystem::output_evolution(std::string &) const
{
	throw Error::Runtime("Method output_evolution of StellarSystem is not implemented yet.");
}

double no_planet_dwconv_dt(double age, const std::valarray<double> &orbit,
		const StellarSystem &system)
{
	throw Error::NotImplemented("StellarSystem::no_planet_dwconv_dt");
/*
	double no_planet_diff_eq[2], no_planet_orbit[2];
	double semimajor=orbit[0]*Rsun_AU,
		   worb=system.get_planet().orbital_angular_velocity_semimajor(
				   semimajor),
		   Iconv=system.get_star().moment_of_inertia(age, convective),
		   dIconv_dt=system.get_star().moment_of_inertia_deriv(age,
				   convective);
	no_planet_orbit[0]=worb*Iconv;
	no_planet_orbit[1]=orbit[1];
	system.no_planet_differential_equation(age, no_planet_orbit,
			no_planet_diff_eq);
	return (no_planet_diff_eq[0] - dIconv_dt*worb)/Iconv;*/
}
