/**\file
 *
 * \brief Functions that evaluate to the various expressions from Lai 2012
 * used for the tests.
 *
 * \ingroup UnitTests_group
 */

#ifndef __LAI_EXPRESSIONS_H
#define __LAI_EXPRESSIONS_H

#include "ConstPhaseLagDissipatingZone.h"
#include "../AstronomicalConstants.h"
#include "../OrbitalExpressions.h"

///Evaluates Lai 2012 equation (35) without the normalization constant.
double dimensionless_torque_x_Lai(double inclination, const Lags &lags);

///Evaluates Lai 2012 equation (27) without the normalization constant.
double dimensionless_torque_z_Lai(double inclination, const Lags &lags);

///Evaluates Lai 2012 equation (28) without the normalization constant.
double dimensionless_power_Lai(double inclination, const Lags &lags);

///Evalutaes Lai 2012 equation (23) in units of
/// \f$M_\odot R_\odot^2 \mathrm{day}^{-1} \mathrm{Gyr}^{-1}\f$.
double torque_norm_Lai(
		///The mass of the body perturbing the dissipating one in
		/// \f$M_\odot\f$
		double perturber_mass, 
		
		///The radius of the body doing the dissipating in \f$R_\odot\f$
		double dissipator_radius, 
		
		///The semimajor axis in \f$R_\odot\f$
		double semimajor);

///Retuns the normalization constant of Lai 2012 equation (28) in units of
/// \f$M_\odot R_\odot^2 \mathrm{day}^{-2} \mathrm{Gyr}^{-1}\f$.
double power_norm_Lai(
		///The mass of the body doing the dissipating in \f$M_\odot\f$
		double dissipator_mass, 

		///The mass of the body doing the perturbing in \f$M_\odot\f$
		double perturber_mass,
		
		///The radius of the body doing the dissipating in \f$R_\odot\f$
		double dissipator_radius, 
		
		///The semimajor axis in \f$R_\odot\f$
		double semimajor);

#endif
