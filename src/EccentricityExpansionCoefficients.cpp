#include "EccentricityExpansionCoefficients.h"

int EccentricityExpansionCoefficients::inner_index(int msign, int s,
		int epower)
{
#ifdef DEBUG
	assert(std::abs(msign)<2);
#endif
	return (epower-s+2*std::min(msign, s-msign))/2;
}

double EccentricityExpansionCoefficients::p_m2s(double e, int s, 
		unsigned max_e_power)
{
	double result=0, e2=std::pow(e, 2);
	int min_n=std::max(1, -s-1), gamma_ind1=s+__max_e_power+2;
	double e_pow=std::pow(e, s+2*min_n);
	for(int gamma_ind2=0; gamma_ind2<=inner_index(-1, s, max_e_power);
			++gamma_ind2) {
		result+=__gamma_minus.at(gamma_ind1).at(gamma_ind2)*e_pow;
		e_pow*=e2;
	}
	return result;
}

double EccentricityExpansionCoefficients::p_0s(double e, int s, 
		unsigned max_e_power)
{
	double result=0, e2=std::pow(e, 2);
	int min_n=std::max(0, -s), alpha_ind1=s+__max_e_power;
	double e_pow=std::pow(e, s+2*min_n);
	for(int alpha_ind2=0; alpha_ind2<=inner_index(0, s, max_e_power);
			++alpha_ind2) {
		result+=__alpha.at(alpha_ind1).at(alpha_ind2)*e_pow;
		e_pow*=e2;
	}
	return result;
}

double EccentricityExpansionCoefficients::p_p2s(double e, int s, 
		unsigned max_e_power)
{
	double result=0, e2=std::pow(e, 2);
	int min_n=std::max(-1, -s+1), gamma_ind1=s+__max_e_power-2;
	double e_pow=std::pow(e, s+2*min_n);
	for(int gamma_ind2=0; gamma_ind2<=inner_index(1, s, max_e_power);
			++gamma_ind2) {
		result+=__gamma_plus.at(gamma_ind1).at(gamma_ind2)*e_pow;
		e_pow*=e2;
	}
	return result;
}

EccentricityExpansionCoefficients::EccentricityExpansionCoefficients(
			const std::string &tabulated_pms_fname, int max_e_power) :
	__max_e_power(max_e_power), __alpha(2*max_e_power+1),
	__gamma_plus(2*max_e_power+1), __gamma_minus(2*max_e_power+1)
{
	std::ifstream tabulated_coef(tabulated_pms_fname.c_str());
	if(!tabulated_coef) throw Error::IO("Unable to open eccentricity "
			"expansion file: "+tabulated_pms_fname+"!");
	for(int epower=0; epower<=max_e_power; ++epower) {
		for(int s=-epower-2; s<=epower-2; s+=2) {
			if(s) {
				std::vector<double> 
					&destination=__gamma_minus.at(s+__max_e_power+2);
				if(destination.size()==0)
					destination.resize(inner_index(-1, s, max_e_power)+1);
				tabulated_coef >> destination.at(inner_index(-1, s, epower));
			}
		}
		for(int s=-epower; s<=epower; s+=2) {
			std::vector<double> &destination=__alpha.at(s+__max_e_power);
			if(destination.size()==0) 
				destination.resize(inner_index(0, s, max_e_power)+1);
			tabulated_coef >> destination.at(inner_index(0, s, epower));
		}
		for(int s=-epower+2; s<=epower+2; s+=2) {
			if(s) {
				std::vector<double> 
					&destination=__gamma_plus.at(s+__max_e_power-2);
				if(destination.size()==0)
					destination.resize(inner_index(1, s, max_e_power)+1);
				tabulated_coef >> destination.at(inner_index(1, s, epower));
			}
		}
	}
}

double EccentricityExpansionCoefficients::operator()(int m, int s, 
		double e, unsigned max_e_power)
{
	if(s<-static_cast<int>(max_e_power)+m || 
			s>static_cast<int>(max_e_power)+m) return 0;
	switch(m) {
		case -2 : return (s==0 ? 0 : p_m2s(e, s, max_e_power));
		case 0  : return p_0s(e, s, max_e_power);
		case 2  : return (s==0 ? 0 : p_p2s(e, s, max_e_power));
		default : throw Error::BadFunctionArguments(
						  "Asking for p_{m,s} with m other than +-2 and 0");
	};
}

int main(int argc, char **argv)
{
	std::cerr.precision(12);
	std::cerr << std::setw(20) << "e"
		<< std::setw(20) << "O(e)"
		<< std::setw(20) << "s"
		<< std::setw(20) << "p-2"
		<< std::setw(20) << "p0"
		<< std::setw(20) << "p+2"
		<< std::endl;
	EccentricityExpansionCoefficients coef(argv[1], 10);
	for(double e=0.3; e<1; e+=0.3) 
		for(int epower=0; epower<10; ++epower)
			for(int s=-epower-2; s<=epower+2; ++s) {
				std::cerr << std::setw(20) << e
					<< std::setw(20) << epower
					<< std::setw(20) << s;
				for(int m=-2; m<=2; m+=2) 
					std::cerr << std::setw(20) << coef(m, s, e, epower);
				std::cerr << std::endl;
			}
}
