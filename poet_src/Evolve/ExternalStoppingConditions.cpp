#include "ExternalStoppingConditions.h"

namespace Evolve {

    std::valarray<double> RotFastCondition::operator()(
        Core::EvolModeType,
        const std::valarray<double> &,
        const std::valarray<double> &,
        std::valarray<double> &stop_deriv
    ) const
    {
        if(!std::isfinite(__spin_thres)) return std::valarray<double>(-1, 1);
        double spin_freq = __zone.spin_frequency();
        stop_deriv.resize(1, Core::NaN);
        return std::valarray<double>(
            (spin_freq - __spin_thres) / __spin_thres,
            1
        );
    }

} //End Evolve namespace.
