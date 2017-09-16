#include "ConstPhaseLagDissipatingZone.h"

double &Lags::operator()(int m, int mp)
{
	return (*this)[std::pair<int, int>(m, mp)];
}

double Lags::operator()(int m, int mp) const
{
	Lags::const_iterator it=this->find(std::pair<int, int>(m, mp));
	assert(it!=this->end());
	return it->second;
}

void ConstPhaseLagDissipatingZone::describe(std::ostream &os) const
{
	os << "I=" << inclination() << ", O(e)=" << eccentricity_order();
	for(Lags::const_iterator i=__lags.begin(); i!=__lags.end(); ++i) 
		os << ", D'_{" << i->first.first << "," << i->first.second << "}="
			<< i->second;
}


