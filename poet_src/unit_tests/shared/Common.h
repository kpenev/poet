/**\file
 *
 * \brief Functions and classes of general use for all unit tests.
 *
 * \ingroup UnitTests_group
 */

#ifndef __UNITTEST_COMMON_H
#define __UNITTEST_COMMON_H

#include "../Evolve/StopInformation.h"
#include "../Core/AstronomicalConstants.h"

#include <cpptest.h>
#include <valarray>
#include <sstream>
#include <cstdlib>
#include <string>

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
const double MIN_STELLAR_MASS = 0.4,

      ///The boundary between high and low mass stars in \f$M_\odot\f$.
      MAX_LOW_MASS = 1.075,

      ///The highest stellar mass to use in tests in \f$M_\odot\f$.  
      MAX_STELLAR_MASS = 1.3,

      ///Most tests start at this age in Gyr.
      MIN_AGE = 1e-7,

      ///Most tests end at this age in Gyr.  
      MAX_AGE = 10.0,

      ///The lower limit of the mass of random planets.
      MIN_PLANET_MASS = 10,

      ///The upper limit of the mass of random planets.
      MAX_PLANET_MASS = 80,

      ///The lower limit of the radius of random planets.
      MIN_PLANET_RADIUS = 5,

      ///The upper limit of the radius of random planets.
      MAX_PLANET_RADIUS = 15;

///\brief Generates a uniformly distributed random number.
///
///Seeding the random number generator is the caller's responsibility.
double uniform_rand(double min, double max);

///Create a string with a description of the given stop info.
std::string stop_info_to_str(const Evolve::StopInformation &stop);

///\brief Outputs the mass and age polynomial defined by the given polynomial
///coefficients array
std::ostream &operator<<(
    std::ostream &os, 
    const std::valarray< std::valarray<double> > &poly_coef
);

///A uniform random real value in the given range.
double rand_value(double min, double max);

///A uniform integer value in the given range.
int rand_value(int min, int max);

///Fills the given valarray with a random set of polynomial coefficients.
void rand_poly_coef(std::valarray< std::valarray<double> > &poly_coef,
                    double max_mass=-1);

///Returns a random set of polynomial coefficients
std::valarray< std::valarray<double> > rand_poly_coef(double max_mass=-1);

///\brief Returns new polynomial coefficienst such that output
///polynomial(mass, age+age_offset)=input polynomial(mass, age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset);

///\brief Returns new polynomial coefficienst such that output
///polynomial(age+age_offset)=input polynomial(age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset);

///Given n, m and (n)C(m) returns (n)C(m+1)
unsigned next_binom_coef(unsigned n, unsigned m, unsigned nCm);


#endif
