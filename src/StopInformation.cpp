#include "StopInformation.h"

std::ostream &operator<<(std::ostream &os, const StopInformation &stop)
{
	std::streamsize orig_precision=os.precision();
	os.precision(16);
	std::ios_base::fmtflags orig_flags=os.flags();
	os.setf(std::ios_base::scientific);
	os << "Stop at t=" << stop.stop_age()
		<< (stop.is_crossing() ? ", crossing" : ", extremum")
		<< " of " << stop.stop_reason() << ", precision="
		<< stop.stop_condition_precision(); 
	os.precision(orig_precision);
	os.flags(orig_flags);
	return os;
}

