/**\file
 *
 * \brief A collection of functions which calculate various quantities for
 * two body orbits.
 *
 */

#ifndef __ORBITAL_EXPRESSIONS_H
#define __ORBITAL_EXPRESSIONS_H

#include <cmath>
#include "AstronomicalConstants.h"

///Returns the orbital angular velocity of the given orbit in rad/day.
double orbital_angular_velocity(
		///The mass of the first body in \f$M_\odot\f$.
		double m1,

		///The mass of the second body in \f$M_\odot\f$.
		double m2,
		
		///The semimajor axis in AU
		double semimajor,
		
		///Whether to return the derivative with respect to the semimajor
		///axis instead of the value.
		bool deriv=false);

///\brief The energy of the orbit (assuming 0 gravitational potential at
///infinity) in \f$\frac{M_\odot R_\odot^2 rad^2}{day^2}\f$.
double orbital_energy(
		///The mass of the first body in \f$M_\odot\f$.
		double m1,

		///The mass of the second body in \f$M_\odot\f$.
		double m2,
		
		///The semimajor axis in AU
		double semimajor);

///\brief The angular momentum of the orbit in 
/// \f$\frac{M_\odot R_\odot^2 rad}{day}\f$.
double orbital_angular_momentum(
		///The mass of the first body in \f$M_\odot\f$.
		double m1,

		///The mass of the second body in \f$M_\odot\f$.
		double m2,
		
		///The semimajor axis in AU
		double semimajor,

		///The eccentricity.
		double eccentricity);

#endif
