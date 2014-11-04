/**\file
 *
 * \brief Functions and classes of general use for all unit tests.
 *
 * \ingroup UnitTests_group
 */

#ifndef __UNITTEST_COMMON_H
#define __UNITTEST_COMMON_H

#include "../AstronomicalConstants.h"
#include <cpptest.h>
#include <valarray>
#include <sstream>

///\brief Returns true iff \f$|x-y|\leq\mathrm{abs\_tolerance} +
/// \mathrm{frac\_tolerance}\cdot\max(|x|,|y|)\f$.
bool check_diff(double x, double y, double frac_tolerance, 
		double abs_tolerance);

///\brief Returns true iff \f$ \forall i\ |x_i-y_i|\leq
/// \mathrm{abs\_tolerance}_i +
/// \mathrm{frac\_tolerance}_i\cdot\max(|x_i|,|y_i|)\f$.
bool check_diff(std::valarray<double> x, std::valarray<double> y,
		std::valarray<double> frac_tolerance,
		std::valarray<double> abs_tolerance);

///\todo Get rid of this function and use check_diff instead.
bool isEqual(double a, double b);

///\todo Get rid of this function and use check_diff instead.
double getError(double predicted, double actual);

///\todo Get rid of this function and use check_diff instead.
bool approxEqual(double predicted, double actual, double thres=0.02);

///The orbital angular momentum corresponding to the given frequency.
double orbital_angmom_from_freq(double m1, double m2, double freq, double e);

///The lowest stellar mass to use in tests in \f$M_\odot\f$.
const double min_stellar_mass=0.4,
	  
	  ///The boundary between high and low mass stars in \f$M_\odot\f$.
	  max_low_mass=1.075,

	  ///The highest stellar mass to use in tests in \f$M_\odot\f$.  
	  max_stellar_mass=1.3,

	  ///Most tests start at this age in Gyr.
	  min_age=1e-7,

	  ///Most tests end at this age in Gyr.  
	  max_age=10.0;

///The lower limit of the mass of random planets.
const double min_planet_mass=10,

	  ///The upper limit of the mass of random planets.
	  max_planet_mass=80,

	  ///The lower limit of the radius of random planets.
	  min_planet_radius=5,

	  ///The upper limit of the radius of random planets.
	  max_planet_radius=15;

///\brief Outputs a comma separated list of the values in the array to the
///given stream.
std::ostream &operator<<(std::ostream &os,
		const std::valarray<double> &array);

#endif
