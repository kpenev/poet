/**\file
 *
 * \brief Declaration of some general purpose utilities.
 * 
 * \ingroup Utilities_group
 */

#ifndef __COMMON_DEFINITIONS_H
#define __COMMON_DEFINITIONS_H

#include "Error.h"
#include "AstronomicalConstants.h"
#include <list>
#include <valarray>
#include <limits>
#include <sstream>
#include <iostream>
#include <iterator>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

///Has no business being here.
///TODO: move it to a more appropriate file.
struct RotationScenario {
	double initSpin;
	double K;
	double Tc;
	double Tdisk;
	RotationScenario(double initSpin, double K, double Tc, double Tdisk):
		initSpin(initSpin), K(K), Tc(Tc), Tdisk(Tdisk) {}
};

///Not a number.
const double NaN=std::numeric_limits<double>::quiet_NaN();

///Infinity
const double Inf=std::numeric_limits<double>::infinity();

///A shortcut for the solar radius in astronomical units.
const double Rsun_AU=AstroConst::solar_radius/AstroConst::AU;

///\brief The various evolution modes.
///
///Each corresponds to a different system of differential equations with
///different variables.
enum EvolModeType {
	///The orbital period is shorter than the stellar spin period.
	FAST_PLANET=-1, 

	///The orbital and the stellar spin periods are locked to be the same.
	LOCKED_TO_PLANET,
	
	///The orbital period is longer than the stellar spin period.
	SLOW_PLANET,
	
	///There is no planet in the system (only the stellar rotation evolves).
	NO_PLANET,

	///The stellar surface rotation is locked to a protoplanetary disk.
	LOCKED_TO_DISK,
	
	///Used as the mode to transform to from all other modes when storing the
	///computed evolution.
	TABULATION};

///More civilized output for EvolModeType variables.
std::ostream &operator<<(std::ostream &os, const EvolModeType &evol_mode);

///Outputs a valarray as a sequence of ', ' separated values.
std::ostream &operator<<(std::ostream &os, std::valarray<double> &arr);

///Creates a valarray containing the values in the given list.
std::valarray<double> list_to_valarray(const std::list<double> &inlist);

///Solves the cubic equation \f$ \sum_{i=0}^3 c_i x^i=0 \f$
std::valarray<double> solve_cubic(double c0, double c1,
		double c2, double c3);

///\brief Finds a zero of a smooth curve defined by the given points and
///derivatives.
///
///If no derivative information is provided, returns the abscissa at which
///the unique straight line passing through two points is zero.
///
///If derivative informaiton is provided, returns one of the zeroes in the
///interval (x0, x1) of a cubic function passing through two points and
///having specified values of its derivatives at those points.
///
///The two function values (y0 and y1) must have opposite sign.
double estimate_zerocrossing(
		///The abscissa of the first point though which the function
		///should pass.
		double x0,

		///The oordinate of the first point though which the function
		///should pass.
		double y0,

		///The abscissa of the second point though which the function
		///should pass.
		double x1,
		
		///The oordinate of the second point though which the
		///function should pass.
		double y1,

		///The derivative of the function at the first point
		///(optional).
		double dy0=NaN,
		
		///The derivative of the function at the second point
		///(optional).
		double dy1=NaN);

///\brief Finds the abscissa of a zero of a quadratic defined by three
///points.
///
///The three points are (x0, y0), (x1, y1), (x2, y2). 
///Two of the function values must have opposite sign.
double quadratic_zerocrossing(
		double x0, double y0, double x1, double y1,
		double x2, double y2, double require_range_low=NaN,
		double require_range_high=NaN);

///\brief Finds the abscissa of a zero of a cubic defined by four points.
///
///The four points are (x0, y0), (x1, y1), (x2, y2), (x3, y3). 
///Two of the function values must have opposite sign.
double cubic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3,
		double require_range_low=NaN, double require_range_high=NaN);


///\brief Finds an extremum of a cubic defined by two points and derivatives.
///
///If the optional last argument is not NULL, the location it points to gets
///overwritten with the value of the function at the extremum.
///
///The two derivatives must have an opposite sign.
double estimate_extremum(
		///The abscissa of the first point though which the function
		///should pass.	
		double x0,
		
		///The oordinate of the first point though which the function
		///should pass.
		double y0,
		
		///The abscissa of the second point though which the function
		///should pass.
		double x1,

		///The oordinate of the second point though which the
		///function should pass.
		double y1,

		///The derivative of the function at the first point
		double dy0,
		
		///The derivative of the function at the second point
		double dy1,
		
		///A location to store the value the function takes at the
		///extremum.
		double *extremum_y=NULL);


///\brief Finds the abscissa of the extremum of the quadratic function
///passing through three points.
///
///The three points are (x0, y0), (x1, y1), (x2, y2).
///
///If the optional last argument is not NULL, the location it points to gets
///overwritten with the value of the function at the extremum.
double quadratic_extremum(double x0, double y0, double x1,
		double y1, double x2, double y2, double *extremum_y=NULL);

///\brief Finds the abscissa of an extremum of the cubic function passing
///through four points.
///
///The four points are (x0, y0), (x1, y1), (x2, y2), (x3, y3).
///
///The input y values must not be monotonic.
///
///If the optional last argument is not NULL, the location it points to gets
///overwritten with the value of the function at the extremum.
double cubic_extremum(double x0, double y0, double x1,
		double y1, double x2,double y2, double x3,	double y3,
		double *extremum_y=NULL, double require_range_low=NaN,
		double require_range_high=NaN);

///\brief Returns the polynomial coefficients of the cubic going through four
///points.
///
///The four points are (x0, y0), (x1, y1), (x2, y2), (x3, y3).
///
///The return value is a [GSL](http://www.gnu.org/software/gsl/) vector with
///the \f$i^{th}\f$ component being the coefficient in front of \f$x^i\f$.
///
///The vector must be freed when no longer needed.
gsl_vector *cubic_coefficients(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3);

///\brief Finds the polynomial that goes through the points iterated over by
///x_i and y_i.
///
///The return value is a [GSL](http://www.gnu.org/software/gsl/) vector with
///the \f$i^{th}\f$ component being the coefficient in front of \f$x^i\f$.
///
///The vector must be freed when no longer needed.
template<class ITERATOR>
gsl_vector *polynomial_coefficients(ITERATOR x_i, ITERATOR y_i,
		size_t num_points)
{
	gsl_vector *y_vec=gsl_vector_alloc(num_points);
	gsl_matrix *xpowers=gsl_matrix_alloc(num_points, num_points);
	for(size_t i=0; i<num_points; i++) {
		double x=*x_i, xpow=1.0;
		gsl_vector_set(y_vec, i, *y_i);
		for(size_t pow=0; pow<num_points; pow++) {
			gsl_matrix_set(xpowers, i, pow, xpow);
			xpow*=x;
		}
		x_i++; y_i++;
	}
	gsl_matrix *xpowers_LU=gsl_matrix_alloc(num_points, num_points);
	gsl_matrix_memcpy(xpowers_LU, xpowers);

	gsl_permutation *permutation=gsl_permutation_alloc(num_points);
	gsl_vector *coefficients=gsl_vector_alloc(num_points);
	int permutation_sign;
	gsl_linalg_LU_decomp(xpowers_LU, permutation, &permutation_sign);
	gsl_linalg_LU_solve(xpowers_LU, permutation, y_vec, coefficients);
	gsl_vector *residuals=gsl_vector_alloc(num_points);
	gsl_linalg_LU_refine(xpowers, xpowers_LU, permutation, y_vec,
			coefficients, residuals);
	gsl_permutation_free(permutation);
	gsl_vector_free(y_vec);
	gsl_vector_free(residuals);
	gsl_matrix_free(xpowers);
	gsl_matrix_free(xpowers_LU);
	return coefficients;
}

#endif
