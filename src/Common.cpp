/**\file
 *
 * \brief The implementation of some of the utility functions.
 * 
 * \ingroup Utilities_group
 */

#include "Common.h"
#include "gsl/gsl_poly.h"
#include <assert.h>

std::ostream &operator<<(std::ostream &os, const EvolModeType &evol_mode)
{
	switch(evol_mode) {
		case FAST_PLANET: os << "FAST_PLANET"; break;
		case LOCKED_TO_PLANET: os << "LOCKED_TO_PLANET"; break;
		case SLOW_PLANET: os << "SLOW_PLANET"; break;
		case NO_PLANET: os << "NO_PLANET"; break;
		case LOCKED_TO_DISK: os << "LOCKED_TO_DISK"; break;
		case TABULATION : os << "TABULATION";
	}
	return os;
}

std::ostream &operator<<(std::ostream &os, std::valarray<double> &arr)
{
	for(size_t i=0; i<arr.size(); i++) {
		if(i) os << ", ";
		os << arr[i];
	}
	return os;
}

std::valarray<double> list_to_valarray(const std::list<double> &inlist)
{
	std::valarray<double> result(inlist.size());
	double *dest=&result[0];
	for(std::list<double>::const_iterator li=inlist.begin(); 
			li!=inlist.end(); li++) {
		*dest = *li;
		dest++;
	}
	return result;
}

std::valarray<double> solve_cubic(double c0, double c1,
		double c2, double c3) {
	double x[3];
	int n = 0;
	if (c3 == 0)
		n = gsl_poly_solve_quadratic(c2, c1, c0, &x[0], &x[1]);
	else
		n = gsl_poly_solve_cubic(c2/c3, c1/c3, c0/c3, &x[0], &x[1], &x[2]);
	std::valarray<double> result(x, n);
	return result;
}

double estimate_zerocrossing(double x0, double y0, double x1,
		double y1, double dy0,	double dy1)
{
	double dx=x1-x0, dy=y1-y0, slope_inv=dx/dy;
	double linear_solution=x0 - y0*slope_inv;
	if(linear_solution<x0 || linear_solution>x1)
		linear_solution=x1 - y1*slope_inv;
	if(std::isnan(dy0) || std::isnan(dy1)) {
#ifdef DEBUG
		std::cerr << "Linear zerocrossing between (" << x0 << ", " << y0
			<< ") and (" << x1 << ", " << y1 << ")=" << linear_solution 
			<< std::endl;
#endif
		return linear_solution;
	} else {
		double x1_2=std::pow(x1, 2), x0_2=std::pow(x0, 2),
			   a=(dx*(dy1+dy0) - 2.0*dy)/std::pow(dx, 3),
			   b=(dy1-dy0 - 3.0*a*(x1_2-x0_2))/(2.0*dx),
			   c=dy0 - 3.0*a*x0_2 - 2.0*b*x0,
			   d=y0 - a*x0*x0_2 - b*x0_2 - c*x0;
		std::valarray<double> solutions=solve_cubic(d, c, b, a);
#ifdef DEBUG
		std::cerr << "Cubic zerocrossing between (" << x0 << ", " << y0
		   << ", " << dy0 << ") and (" << x1 << ", " << y1 << ", " << dy1
		   << "), coef=(" << a << ", " << b << ", " << c << ", " << d
		   << "), solutions=(" << solutions << ")";
#endif
		for(size_t i=0; i<solutions.size(); i++) 
			if(solutions[i]>=x0 && solutions[i]<=x1) {
#ifdef DEBUG
				std::cerr << ", selected: " << solutions[i] << std::endl;
#endif
				return solutions[i];
			}
		if(y0*y1<=0) {
#ifdef DEBUG
			std::cerr << ", fallback to linear: " << linear_solution
				<< std::endl;
#endif
			return linear_solution;
		}
		std::ostringstream msg;
		msg.precision(16);
		msg << "No solution to " << a << "*x**3 + " << b << "*x**2 + " << c
			<< "*x + " << d << " (x0=" << x0 << ", y0=" << y0 << ", x1=" << x1
			<< ", y1=" << y1 << ", dy0=" << dy0 << ", dy1=" << dy1
			<< " found in the given range in estimate zerocrossing. "
			"Solutions(" << solutions.size() << "): ";
		if(solutions.size()>0) msg << solutions[0];
		for(size_t i=1; i<solutions.size(); i++) 
			msg << ", " << solutions[i];
		throw Error::BadFunctionArguments(msg.str());
	}
}

double quadratic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2, double require_range_low,
		double require_range_high)
{
	double d12=(y1-y2)/(x1-x2), d02=(y0-y2)/(x0-x2), dx10=(x1-x0),
		   a=d12/(x1-x0) - d02/dx10,
		   b=d12 + d02*(x1+x2)/dx10 - d12*(x1+x2)/dx10,
		   c=y2 + d12*x1*x2/dx10 - d12*x2 - d02*x1*x2/dx10,
		   xpre, ypre, xpost, ypost;
	if(y0*y1<0) {xpre=x0; xpost=x1; ypre=y0; ypost=y1;}
	else if(y1*y2<0) {xpre=x1; xpost=x2; ypre=y1; ypost=y2;}
	else {
		std::ostringstream msg;
		msg.precision(16);
		msg << "No sign change found in quadratic_zerocrossing(x0=" 
			<< x0 << ", y0=" << y0 << ", x1=" << x1 << ", y1=" << y1 
			<< ", x2=" << x2 << ", y2=" << y2 << ")";
		throw Error::BadFunctionArguments(msg.str());
	}
	std::valarray<double> solutions=solve_cubic(c, b, a, 0);
	if(std::isnan(require_range_low)) {
		assert(std::isnan(require_range_high));
		require_range_low=x0;
		require_range_high=x2;
	}
#ifdef DEBUG
	std::cerr << "Quadratic zerocrossing between (" << x0 << ", " << y0 
		<< "), (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2
		<< ") in range (" << require_range_low << ", " << require_range_high
		<< "), coef=(" << a << ", " << b << ", " << c << "), solutions=("
		<< solutions << ")";
#endif
	for(size_t i=0; i<solutions.size(); i++) 
		if(solutions[i]>=require_range_low &&
				solutions[i]<=require_range_high) {
#ifdef DEBUG
			std::cerr << ", selected: " << solutions[i] << std::endl;
#endif
			return solutions[i];
		}
#ifdef DEBUG
	std::cerr << ", fallback to linear: ";
#endif
	return estimate_zerocrossing(xpre, ypre, xpost, ypost);
}

double cubic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3,
		double require_range_low, double require_range_high)
{
	double xpre, ypre, xpost, ypost;
	if(y0*y1<=0) {xpre=x0; xpost=x1; ypre=y0; ypost=y1;}
	else if(y1*y2<=0) {xpre=x1; xpost=x2; ypre=y1; ypost=y2;}
	else if(y2*y3<=0) {xpre=x2; xpost=x3; ypre=y2; ypost=y3;}
	else {
		std::ostringstream msg;
		msg.precision(16);
		msg << "No sign change found in cubic_zerocrossing(x0=" 
			<< x0 << ", y0=" << y0 << ", x1=" << x1 << ", y1=" << y1 
			<< ", x2=" << x2 << ", y2=" << y2 
			<< ", x3=" << x3 << ", y3=" << y3 << ")";
		throw Error::BadFunctionArguments(msg.str());
	}
	gsl_vector *cubic_coef=cubic_coefficients(x0,y0, x1,y1, x2,y2, x3,y3);
	std::valarray<double> solutions=solve_cubic(
			gsl_vector_get(cubic_coef, 0),
			gsl_vector_get(cubic_coef, 1),
			gsl_vector_get(cubic_coef, 2),
			gsl_vector_get(cubic_coef, 3));
	if(std::isnan(require_range_low)) {
		assert(std::isnan(require_range_high));
		require_range_low=x0;
		require_range_high=x3;
	}
#ifdef DEBUG
	std::cerr << "Cubic zerocrossing between (" << x0 << ", " << y0 
		<< "), (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2
		<< "), (" << x3 << ", " << y3 << ") in range ("
		<< require_range_low << ", " << require_range_high << "), coef=(";
	for(int i=3; i>=0; i--) {
		std::cerr << gsl_vector_get(cubic_coef, i);
		if(i>0) std::cerr << ", ";
	}
	std::cerr << "), solutions=(" << solutions << ")";
#endif
	gsl_vector_free(cubic_coef);
	for(size_t i=0; i<solutions.size(); i++) 
		if(solutions[i]>=require_range_low &&
				solutions[i]<=require_range_high) {
#ifdef DEBUG
			std::cerr << ", selected: " << solutions[i] << std::endl;
#endif
			return solutions[i];
		}
#ifdef DEBUG
	std::cerr << ", fallback to linear: ";
#endif
	return estimate_zerocrossing(xpre, ypre, xpost, ypost);
}

double estimate_extremum(double x0, double y0, double x1,
		double y1, double dy0,	double dy1, double *extremum_y)
{
	assert(dy0*dy1<0);
	double dx=x1-x0;
	double x1_2=std::pow(x1, 2), x0_2=std::pow(x0, 2),
			   a=(dx*(dy1+dy0) - 2.0*(y1-y0))/std::pow(dx, 3),
			   b=(dy1-dy0 - 3.0*a*(x1_2-x0_2))/(2.0*dx),
			   c=dy0 - 3.0*a*x0_2 - 2.0*b*x0,
			   d=y0 - a*x0*x0_2 - b*x0_2 - c*x0,
			   linear_extremum_x=x0 - dy0*dx/(dy1-dy0),
			   linear_extremum_y=(linear_extremum_x-x0<x1-linear_extremum_x ?
					   y0+dy0*(linear_extremum_x-x0) :
					   y1-dy1*(x1-linear_extremum_x));
	std::valarray<double> solutions=solve_cubic(c, 2.0*b, 3.0*a, 0);
#ifdef DEBUG
	std::cerr << "Cubic extrema between (" << x0 << ", " << y0 << ", " << dy0
		<< ") and (" << x1 << ", " << y1 << ", " << dy1 << "), coef=("
		<< a << ", " << b << ", " << c << ", " << d
		<< "), solutions: (" << solutions[0] << ", "
		<< a*std::pow(solutions[0], 3) + b*std::pow(solutions[0], 2) + 
		   c*solutions[0] + d << "), (" << solutions[1] << ", "
		<< a*std::pow(solutions[1], 3) + b*std::pow(solutions[1], 2) + 
		   c*solutions[1] + d << ")";
#endif
	for(size_t i=0; i<solutions.size(); i++) {
		if(solutions[i]>=x0 && solutions[i]<=x1) {
			double x=solutions[i];
			if(extremum_y!=NULL) {
				double ax3=a*x*x*x, bx2=b*x*x, cx=c*x;
				*extremum_y=ax3 + bx2 + cx + d;
				if(std::abs(*extremum_y)/std::max(
							std::max(std::abs(ax3), std::abs(bx2)),
							std::max(std::abs(cx), std::abs(d)))<1e-10)
					*extremum_y=linear_extremum_y;
			}
#ifdef DEBUG
			std::cerr << ", selected: (" << x;
			if(extremum_y) std::cerr << ", " << *extremum_y;
			std::cerr << ")" << std::endl;
#endif
			return x;
		}
	}
	if(extremum_y!=NULL) *extremum_y=linear_extremum_y;
#ifdef DEBUG
	std::cerr << ", falling back to linear: (" << linear_extremum_x << ", "
		<< linear_extremum_y << ")" << std::endl;
#endif
	return linear_extremum_x;
}

double quadratic_extremum(double x0, double y0, double x1,
		double y1, double x2, double y2, double *extremum_y)
{
	assert((y1-y0)*(y2-y1)<=0);
	double s02=(y0-y2)/(x0-x2), s12=(y1-y2)/(x1-x2), a=(s02 - s12)/(x0-x1),
		   extremum_x=0.5*(x0 + x2 - s02/a);
	if(extremum_y)
		*extremum_y=y0 - s02*s02/(2.0*a) + s02*x2 
			- a*(x0*x0 + x2*x2)/2.0;
#ifdef DEBUG
	double b=-2.0*a*extremum_x;
	std::cerr << "Quadratic extremum between (" << x0 << ", " << y0 << "), ("
		<< x1 << ", " << y1 << "), (" << x2 << ", " << y2 << "), coef=("
		<< a << ", " << b << ", "
		<< *extremum_y - a*extremum_x*extremum_x - b*extremum_x
		<< "): (" << extremum_x << ", " << *extremum_y << ")" << std::endl;
#endif
	return extremum_x;
}

double cubic_extremum(double x0, double y0, double x1,
		double y1, double x2,double y2, double x3,	double y3,
		double *extremum_y, double require_range_low,
		double require_range_high)
{
	assert((y1-y0)*(y2-y1)<=0 || (y2-y1)*(y3-y2)<=0);
	gsl_vector *cubic_coef=cubic_coefficients(x0,y0, x1,y1, x2,y2, x3,y3);
	double a=3.0*gsl_vector_get(cubic_coef, 3),
		   b=2.0*gsl_vector_get(cubic_coef, 2),
		   c=gsl_vector_get(cubic_coef, 1), sqrtD=std::sqrt(b*b-4.0*a*c);
#ifdef DEBUG
	double extremum_x1=(-b-sqrtD)/(2.0*a), extremum_x2=(-b+sqrtD)/(2.0*a),
		   d=y0 - a*x0*x0*x0/3.0 - b*x0*x0/2.0 - c*x0;
	if(std::isnan(require_range_low)) {
		assert(std::isnan(require_range_high));
		require_range_low=x0;
		require_range_high=x3;
	}
	std::cerr << "Cubic extrema between (" << x0 << ", " << y0 << "), ("
		<< x1 << ", " << y1 << "), (" << x2 << ", " << y2 << "), (" 
		<< x3 << ", " << y3 << ") in range (" << require_range_low << ", "
		<< require_range_high << "), coef=("
		<< a/3.0 << ", " << b/2.0 << ", " << c << ", " << d
		<< "), solutions: (" << extremum_x1 << ", " << 
		a*std::pow(extremum_x1,3)/3.0 + b*std::pow(extremum_x1,2)/2.0 +
		c*extremum_x1 + d << ") and (" << extremum_x2 << ", " 
		<< a*std::pow(extremum_x2,3)/3.0 + b*std::pow(extremum_x2,2)/2.0 +
		c*extremum_x2 + d << ")";
#endif
	if(a<0) {a*=-1; b*=-1; c*=-1;}
	double extremum_x=(-b-sqrtD)/(2.0*a);
	if(extremum_x<require_range_low || extremum_x>require_range_high)
		extremum_x=(-b+sqrtD)/(2.0*a);
	if(extremum_x<require_range_low || extremum_x>require_range_high) {
#ifdef DEBUG
		std::cerr << ", fallback to quadratic: ";
#endif
		if((y1-y0)*(y2-y1)<=0)
			extremum_x=quadratic_extremum(x0, y0, x1, y1, x2, y2,
					extremum_y);
		else extremum_x=quadratic_extremum(x1, y1, x2, y2, x3, y3,
				extremum_y);
		if(extremum_x<require_range_low || extremum_x>require_range_high) {
			std::ostringstream msg;
			msg.precision(16);
			msg.setf(std::ios_base::scientific);
			msg << "No extremum found in the range (" << require_range_low
				<< ", " << require_range_high << ") from points ("
				<< x0 << ", " << y0 << "), ("
				<< x1 << ", " << y1 << "), ("
				<< x2 << ", " << y2 << ") and ("
				<< x3 << ", " << y3 << ")!";
			throw Error::BadFunctionArguments(msg.str());
		}
	} else if(extremum_y) {
		*extremum_y=0;
		double xpow=1.0;
		for(size_t i=0; i<4; i++) {
			*extremum_y+=xpow*gsl_vector_get(cubic_coef, i);
			xpow*=extremum_x;
		}
	}
	gsl_vector_free(cubic_coef);
#ifdef DEBUG
	std::cerr << ", selected: (" << extremum_x;
	if(extremum_y) std::cerr << ", " << *extremum_y;
	std::cerr << ")" << std::endl;
#endif
	return extremum_x;
}

gsl_vector *cubic_coefficients(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3)
{
	std::list<double> xs, ys;
	xs.push_back(x0); ys.push_back(y0);
	xs.push_back(x1); ys.push_back(y1);
	xs.push_back(x2); ys.push_back(y2);
	xs.push_back(x3); ys.push_back(y3);
	return polynomial_coefficients(xs.begin(), ys.begin(),4);
}
