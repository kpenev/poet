"Star.h"

///Converts lg(Q) to a tidal phase lag.
double lag_from_lgQ(double lgQ)
{
	return 15.0 / (16.0 * M_PI * std::pow(10.0, lgQ));
}
