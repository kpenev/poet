#include "SinCosIntegrals.h"

SinCosIntegrals::SinCosIntegrals(unsigned pre_allocate)
	: __pre_allocate(pre_allocate)
{
	__values.reserve(pre_allocate);
	__values.push_back(std::vector<double>(1, 2.0*M_PI));
	__values[0].reserve(pre_allocate);
}

double SinCosIntegrals::operator()(unsigned m, unsigned n)
{
	if(n>m) {
		unsigned temp=n;
		n=m;
		m=temp;
	}
	double twice_m;
	for(unsigned val_m_ind=__values.size(); val_m_ind<=m; ++val_m_ind) {
		twice_m=static_cast<double>(2*val_m_ind);
		__values.push_back((twice_m-1)/twice_m*__values.back());
	}
	std::vector<double> &val_m=__values[m];
	for(unsigned val_n_ind=val_m.size(); val_n_ind<=n; ++val_n_ind) {
		double twice_n=static_cast<double>(2*val_n_ind);
		val_m.push_back((twice_n-1)/(twice_n+twice_m)*val_m.back());
	}
	return __values[m][n];
}
