#ifndef __COMMON_DEFINITIONS_H
#define __COMMON_DEFINITIONS_H

#include "Error.h"
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

struct RotationScenario {
	double initSpin;
	double K;
	double Tc;
	double Tdisk;
	RotationScenario(double initSpin, double K, double Tc, double Tdisk):
		initSpin(initSpin), K(K), Tc(Tc), Tdisk(Tdisk) {}
};

const double NaN=std::numeric_limits<double>::quiet_NaN();
const double Inf=std::numeric_limits<double>::infinity();

///Creates a valarray containing the values in the given list.
std::valarray<double> list_to_valarray(const std::list<double> &inlist);

///solves a cubic polynomial where c0 is the coefficient of x^0, c1 the
///coefficient of x^1, etc.
std::valarray<double> solve_cubic(double c0, double c1,
		double c2, double c3);

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will cross zero given values (and optionally derivatives)
///of the quantity at two locations. The two function values (y0 and y1)
///must have opposite sign.
double estimate_zerocrossing(double x0, double y0, double x1, double y1,
		double dy0=NaN, double dy1=NaN);

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will cross zero given values of the quantity at three locations.
///Two of the function values must have opposite sign.
double quadratic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2, double require_range_low=NaN,
		double require_range_high=NaN);

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will cross zero given values of the quantity at four locations.
///Two of the function values must have opposite sign.
double cubic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3,
		double require_range_low=NaN, double require_range_high=NaN);

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will have an extremum given values and derivatives of the 
///quantity at two locations. The two derivative values (dy0 and dy1)
///must have opposite sign. Additionally, if extremum_y is not NULL it
///gets filled with the best estimate of the value of the function at
///the extremum
double estimate_extremum(double x0, double y0, double x1,
		double y1, double dy0, double dy1, double *extremum_y=NULL);

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will have an extremum given values of the quantity at three
///locations. The input y values must not be monotonic. Additionally, if
///extremum_y is not NULL it gets filled with the best estimate of the value
///of the function at the extremum
double quadratic_extremum(double x0, double y0, double x1,
		double y1, double x2, double y2, double *extremum_y=NULL);

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will have an extremum given values of the quantity at four
///locations. The input y values must not be monotonic. Additionally, if
///extremum_y is not NULL it gets filled with the best estimate of the value
///of the function at the extremum. If two extrema exist in the range
///(x0,x2), the one with the smaller absicissa is returned.
double cubic_extremum(double x0, double y0, double x1,
		double y1, double x2,double y2, double x3,	double y3,
		double *extremum_y=NULL, double require_range_low=NaN,
		double require_range_high=NaN);

///Returns the polynomial coefficients of the cubic going through the given
///points as a GSL vector where C_i is the coefficient in front of x^i. The
///vector must be freed when no longer needed.
gsl_vector *cubic_coefficients(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3);

///Returns the coefficients of the polynomial that goes through the points
///iterated over by x_i and y_i.
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
	return coefficients;
}

#endif
