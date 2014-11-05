#include "Common.h"

bool check_diff(double x, double y, double frac_tolerance, 
		double abs_tolerance)
{
	if(std::isnan(x)) return std::isnan(y);
	return std::abs(x-y)<=(abs_tolerance+
		frac_tolerance*std::max(std::abs(x), std::abs(y)));
}

bool check_diff(std::valarray<double> x, std::valarray<double> y,
		std::valarray<double> frac_tolerance,
		std::valarray<double> abs_tolerance)
{
	bool result=true;
	for(size_t i=0; i<x.size(); i++)
		result=result && (std::abs(x[i]-y[i])<=
				(abs_tolerance[i]+
				frac_tolerance[i]*std::max(std::abs(x[i]), std::abs(y[i]))));
	return result;
}

std::ostream &operator<<(std::ostream &os,
		const std::valarray<double> &array)
{
	os << array[0];
	for(size_t i=1; i<array.size(); i++)
		os << ", " << array[i];
	return os;
}

bool isEqual(double a, double b) {
	return std::abs(a-b) < 1e-20;
}

double getError(double predicted, double actual) {
	return std::abs((actual - predicted)/actual);
}

bool approxEqual(double predicted, double actual, double thres) {
	if (isEqual(predicted, actual)) return true;
	return getError(predicted, actual) < 0.02;
}

double orbital_angmom_from_freq(double m1, double m2, double freq, double e)
{
	return m1*m2*AstroConst::day
		   /std::pow(AstroConst::solar_radius, 2)
		   *std::pow(std::pow(AstroConst::G*AstroConst::solar_mass, 2)
				   	 *AstroConst::day/((m1+m2)*freq), 1.0/3.0)
		   *std::sqrt(1.0-e*e);
}

double uniform_rand(double min, double max)
{
	return (max-min)*static_cast<double>(std::rand())/RAND_MAX+min;
}
