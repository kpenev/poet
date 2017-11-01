#define BUILDING_LIBRARY
#include "StopInformation.h"

namespace Evolve {

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
        if(stop.is_crossing())
            os << ", " << (stop.crossed_zero() ? "after" : "before")
                << " crossing with deriv sign=" << stop.deriv_sign_at_crossing();
        os.precision(orig_precision);
        os.flags(orig_flags);
        return os;
    }

}//End Evolve namespace.
