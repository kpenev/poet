#include "Common.h"

bool check_diff(double x, double y, double frac_tolerance, 
		double abs_tolerance)
{
	if(std::isnan(x)) return std::isnan(y);
	return std::abs(x-y)<=(abs_tolerance+
		frac_tolerance*std::min(std::abs(x), std::abs(y)));
}

bool check_diff(std::valarray<double> x, std::valarray<double> y,
		std::valarray<double> frac_tolerance,
		std::valarray<double> abs_tolerance)
{
	bool result=true;
	for(size_t i=0; i<x.size(); i++)
		result=result && (std::abs(x[i]-y[i])<=
				(abs_tolerance[i]+
				frac_tolerance[i]*std::min(std::abs(x[i]), std::abs(y[i]))));
	return result;
}

bool isEqual(double a, double b) {
	return std::abs(a-b) < 1e-20;
}

double getError(double predicted, double actual) {
	return (std::abs(actual - predicted)
            /
            std::max(std::abs(actual), std::abs(predicted)));
}

bool approxEqual(double predicted, double actual, double thres)
{
    if (isEqual(predicted, actual)) return true;
    return getError(predicted, actual) < thres;
}

double orbital_angmom_from_freq(double m1, double m2, double freq, double e)
{
    using namespace Core;
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

/*std::string stop_info_to_str(const Evolve::StopInformation &stop)
{
	std::ostringstream result;
	result << stop << std::endl;
	return result.str();
}*/

double rand_value(double min, double max)
{
	return (max-min)*rand()/RAND_MAX+min;
}

int rand_value(int min, int max)
{
	return rand()%(max-min+1)+min;
}

std::ostream &operator<<(
    std::ostream &os, 
    const std::valarray< std::valarray<double> > &poly_coef
)
{
	for(size_t age_i=0; age_i<poly_coef.size(); age_i++) {
		for(size_t mass_i=0; mass_i<poly_coef[age_i].size(); mass_i++) {
			os << poly_coef[age_i][mass_i];
			if(age_i) os << "*age";
			if(age_i>1) os << "**" << age_i;
			if(mass_i) os << "*m";
			if(mass_i>1) os << "**" << mass_i;
			if(age_i<poly_coef.size()-1 || mass_i<poly_coef.size()-1)
				os << " + ";
		}
	}
	return os;
}

void rand_poly_coef(std::valarray< std::valarray<double> > &poly_coef,
                    double max_mass)
{
	if (max_mass < 0) max_mass = MAX_STELLAR_MASS;
	poly_coef.resize(3, std::valarray<double>(3));
	double age_fac=1;
	for (unsigned age_i=0; age_i < poly_coef.size(); age_i++) {
		if (age_i > 1) age_fac *= age_i * MAX_AGE;
		double mass_fac = 1;
		for(unsigned m_i = 0; m_i < poly_coef[age_i].size(); ++m_i) {
			if(m_i>1) mass_fac *= m_i * max_mass;
			poly_coef[age_i][m_i] = (static_cast<double>(rand()) / RAND_MAX
                                     /
                                     mass_fac / age_fac);
		}
	}
}

std::valarray< std::valarray<double> > rand_poly_coef(double max_mass)
{
	std::valarray< std::valarray<double> > poly_coef(std::valarray<double>(3), 3);
	rand_poly_coef(poly_coef, max_mass);
	return poly_coef;
}

std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset)
{
	std::valarray< std::valarray<double> > result(
			std::valarray<double>(poly_coef[0].size()), poly_coef.size());
	for(size_t i=0; i<poly_coef.size(); i++) 
		for(size_t j=0; j<poly_coef[0].size(); j++) {
			double c_ij=poly_coef[i][j], offset_k=1.0;
			size_t binom_coef=1;
			for(size_t k=0; k<=i; k++) {
				result[i-k][j]+=c_ij*offset_k*binom_coef;
				offset_k*=age_offset;
				if(k<i) binom_coef = next_binom_coef(i, k, binom_coef);
			}
		}
	return result;
}

unsigned next_binom_coef(unsigned n, unsigned m, unsigned nCm)
{
	assert(n>=m+1);
	return (nCm*(n-m))/(m+1);
}

double lag_from_lgQ(double lgQ)
{
    return 15.0 / (16.0 * M_PI * std::pow(10.0, lgQ));
}

double lag_from_lgQ(double lgQ, double mass_ratio)
{
    return (
        15.0
        /
        (8.0 * M_PI * std::sqrt(1.0 + mass_ratio) * std::pow(10.0, lgQ))
    );
}
