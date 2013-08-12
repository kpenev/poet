#include "Common.h"
#include "gsl/gsl_poly.h"
#include <assert.h>

///Creates a valarray containing the values in the given list.
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

///solves a cubic polynomial where c0 is the coefficient of x^0, c1 the
///coefficient of x^1, etc.
std::valarray<double> solve_cubic(double c0, double c1,
		double c2, double c3) {
	/*c0 is the coefficient of x^0, c1 the coefficient of x^1, etc.*/
	double x[3];
	int n = 0;
	if (c3 == 0)
		n = gsl_poly_solve_quadratic(c2, c1, c0, &x[0], &x[1]);
	else
		n = gsl_poly_solve_cubic(c2/c3, c1/c3, c0/c3, &x[0], &x[1], &x[2]);
	std::valarray<double> result(x, n);
	return result;
}

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will cross zero given values (and optionally derivatives)
///of the quantity at two locations. The two function values (y0 and y1)
///must have opposite sign.
double estimate_zerocrossing(double x0, double y0, double x1,
		double y1, double dy0,	double dy1)
{
	double dx=x1-x0;
	double linear_solution=x0 - y0*dx/(y1-y0);
	if(std::isnan(dy0) || std::isnan(dy1)) return linear_solution;
	else {
		double x1_2=std::pow(x1, 2), x0_2=std::pow(x0, 2),
			   a=(dx*(dy1+dy0) - 2.0*(y1-y0))/std::pow(dx, 3),
			   b=(dy1-dy0 - 3.0*a*(x1_2-x0_2))/(2.0*dx),
			   c=dy0 - 3.0*a*x0_2 - 2.0*b*x0,
			   d=y0 - a*x0*x0_2 - b*x0_2 - c*x0;
		std::valarray<double> solutions=solve_cubic(d, c, b, a);
		for(size_t i=0; i<solutions.size(); i++) 
			if(solutions[i]>=x0 && solutions[i]<=x1) return solutions[i];
		if(y0*y1<0) return linear_solution;
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

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will cross zero given values of the quantity at three locations.
///Two of the function values must have opposite sign.
double quadratic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2)
{
	double d12=(y1-y2)/(x1-x2), d02=(y0-y2)/(x0-x2), dx10=(x1-x0),
		   a=d12/(x1-x0) - d02/dx10,
		   b=d12 + d02*(x1+x2)/dx10 - d12*(x1+x2)/dx10,
		   c=y2 + d12*x1*x2/dx10 - d12*x2 - d02*x1*x2/dx10,
		   xpre, ypre, xpost, ypost;
	if(y0*y1<0) {xpre=x0; xpost=x1; ypre=y0; ypost=y1;}
	if(y1*y2<0) {xpre=x2; xpost=x2; ypre=y1; ypost=y2;}
	else {
		std::ostringstream msg;
		msg.precision(16);
		msg << "No sign change found in quadratic_zerocrossing(x0=" 
			<< x0 << ", y0=" << y0 << ", x1=" << x1 << ", y1=" << y1 
			<< ", x2=" << x2 << ", y2=" << y2 << ")";
		throw Error::BadFunctionArguments(msg.str());
	}
	std::valarray<double> solutions=solve_cubic(c, b, a, 0);
	for(size_t i=0; i<solutions.size(); i++) 
		if(solutions[i]>=x0 && solutions[i]<=x1) return solutions[i];
	return xpre - ypre*(xpost-xpre)/(ypost-ypre);
}

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will cross zero given values of the quantity at four locations.
///Two of the function values must have opposite sign.
double cubic_zerocrossing(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3)
{
	assert(y0*y1<0 || y1*y2<0 || y2*y3<0);
	gsl_vector *cubic_coef=cubic_coefficients(x0,y0, x1,y1, x2,y2, x3,y3);
	std::valarray<double> solutions=solve_cubic(
			gsl_vector_get(cubic_coef, 0),
			gsl_vector_get(cubic_coef, 1),
			gsl_vector_get(cubic_coef, 2),
			gsl_vector_get(cubic_coef, 3));
	gsl_vector_free(cubic_coef);
	for(size_t i=0; i<solutions.size(); i++) 
		if(solutions[i]>=x0 && solutions[i]<=x3) return solutions[i];
	assert(false);
}

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will have an extremum given values and derivatives of the 
///quantity at two locations. The two derivative values (dy0 and dy1)
///must have opposite sign. Additionally, if extremum_y is not NULL it
///gets filled with the best estimate of the value of the function at
///the extremum
double estimate_extremum(double x0, double y0, double x1,
		double y1, double dy0,	double dy1, double *extremum_y)
{
	assert(dy0*dy1<0);
	double dx=x1-x0;
	double x1_2=std::pow(x1, 2), x0_2=std::pow(x0, 2),
			   a=(dx*(dy1+dy0) - 2.0*(y1-y0))/std::pow(dx, 3),
			   b=(dy1-dy0 - 3.0*a*(x1_2-x0_2))/(2.0*dx),
			   c=dy0 - 3.0*a*x0_2 - 2.0*b*x0,
			   d=y0 - a*x0*x0_2 - b*x0_2 - c*x0;
	std::valarray<double> solutions=solve_cubic(c, 2.0*b, 3.0*a, 0);
	for(size_t i=0; i<solutions.size(); i++) {
		if(solutions[i]>=x0 && solutions[i]<=x1) {
			double x=solutions[i];
			if(extremum_y!=NULL) *extremum_y=a*x*x*x + b*x*x + c*x + d;
			return x;
		}
	}
	if(dy0*dy1<0) {
		double extremum=x0 - dy0*dx/(dy1-dy0);
		if(extremum-x0<x1-extremum) *extremum_y=y0+dy0*(extremum-x0);
		else *extremum_y=y1-dy1*(x1-extremum);
		return extremum;
	}
	std::ostringstream msg;
	msg << "No solution to " << a << "*x**2 + " << b << "*x + " << c
		<< "(x0=" << x0 << ", y0=" << y0 << ", x1=" << x1
		<< ", y1=" << y1 << ", dy0=" << dy0 << ", dy1=" << dy1
		<< " found in the given range in estimate zerocrossing.";
	throw Error::BadFunctionArguments(msg.str());
}

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will have an extremum given values of the quantity at three
///locations. The input y values must not be monotonic. Additionally, if
///extremum_y is not NULL it gets filled with the best estimate of the value
///of the function at the extremum
double quadratic_extremum(double x0, double y0, double x1,
		double y1, double x2, double y2, double *extremum_y)
{
	assert((y1-y0)*(y2-y1)<=0);
	double s02=(y0-y2)/(x0-x2), a=(x0-x1)*(s02 - (y1-y2)/(x1-x2)),
		   extremum_x=0.5*(x0 + x2 - s02/a);
	if(extremum_y)
		*extremum_y=y0 - s02*s02/(2.0*a) + s02*x2 
			- a*(x0*x0 + x2*x2)/2.0;
	return extremum_x;
}

///Finds the best estimate of the abscissa at which a smoothly varying
///quantity will have an extremum given values of the quantity at four
///locations. The input y values must not be monotonic. Additionally, if
///extremum_y is not NULL it gets filled with the best estimate of the value
///of the function at the extremum. If two extrema exist in the range
///(x0,x2), the one with the smaller absicissa is returned.
double cubic_extremum(double x0, double y0, double x1,
		double y1, double x2,double y2, double x3,	double y3,
		double *extremum_y)
{
	assert((y1-y0)*(y2-y1)<=0 || (y2-y1)*(y3-y2)<=0);
	gsl_vector *cubic_coef=cubic_coefficients(x0,y0, x1,y1, x2,y2, x3,y3);
	double a=gsl_vector_get(cubic_coef, 3), b=gsl_vector_get(cubic_coef, 2),
		   c=gsl_vector_get(cubic_coef, 1), D=b*b-4.0*a*c;
	if(a<0) {a*=-1; b*=-1; c*=-1;}
	double extremum_x=(-b-std::sqrt(D))/(2.0*a);
	if(extremum_x<x0 || extremum_x>x2)
		extremum_x=(-b+std::sqrt(D))/(2.0*a);
	if(extremum_y) {
		*extremum_y=0;
		double xpow=1.0;
		for(size_t i=0; i<4; i++) {
			*extremum_y+=xpow*gsl_vector_get(cubic_coef, i);
			xpow*=extremum_x;
		}
	}
	gsl_vector_free(cubic_coef);
	return extremum_x;
}

///Returns the polynomial coefficients of the cubic going through the given
///points as a GSL vector where C_i is the coefficient in front of x^i. The
///vector must be freed when no longer needed.
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
