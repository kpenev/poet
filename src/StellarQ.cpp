#include "StellarQ.h"

double Star::get_tidal_Q(double tidal_frequency) const
{
	if(tidal_frequency==Inf) 
		return tidal_Q;
	else if(tidal_frequency==-Inf)
		return -tidal_Q;
	return (std::abs(tidal_frequency)<Q_transition_width ?
		tidal_Q/std::sin(tidal_frequency/Q_transition_width*M_PI/2.0) :
		(tidal_frequency>0 ? tidal_Q : -tidal_Q));
}

double Star::get_tidal_Q_deriv(double tidal_frequency) const
{
	if(std::abs(tidal_frequency)<Q_transition_width) {
		double sin_arg=tidal_frequency/Q_transition_width*M_PI/2.0;
		return -tidal_Q/std::pow(std::sin(sin_arg),2)*
				std::cos(sin_arg)*M_PI/(2.0*Q_transition_width);
	} else return 0.0;
}
