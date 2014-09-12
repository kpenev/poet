#include "DissipationQuantities.h"

std::ostream &operator<<(std::ostream &os, 
		const Dissipation::Quantity &quantity)
{
	switch(quantity) {
		case Dissipation::POWER : os << "POWER"; break;
		case Dissipation::TORQUEX : os << "TORQUEX"; break;
		case Dissipation::TORQUEZ : os << "TORQUEZ"; break;
		case Dissipation::SEMIMAJOR_DECAY : os << "SEMIMAJOR_DECAY"; break;
		case Dissipation::ORBIT_SPINUP : os << "ORBIT_SPINUP"; break;
		case Dissipation::INCLINATION_DECAY : os << "INCLINATION_DECAY";
											  break;
		default :
#ifdef DEBUG
				  assert(false)
#endif
					  ;
	};
	return os;
}

///More civilized output for Dissipation::Derivative variables.
std::ostream &operator<<(std::ostream &os,
		const Dissipation::Derivative &deriv)
{
	switch(deriv) {
		case Dissipation::NO_DERIV : os << "NO_DERIV"; break;
		case Dissipation::AGE : os << "AGE"; break;
		case Dissipation::RADIUS : os << "RADIUS"; break;
		case Dissipation::MOMENT_OF_INERTIA :
								   os << "MOMENT_OF_INERTIA";
								   break;
		case Dissipation::SPIN_ANGMOM : os << "SPIN_ANGMOM"; break;
		case Dissipation::SEMIMAJOR : os << "SEMIMAJOR"; break;
		case Dissipation::INCLINATION : os << "INCLINATION"; break;
		default :
#ifdef DEBUG
				  assert(false)
#endif
					  ;
	};
	return os;
}
