/**\file
 *
 * \brief Definitions of some of the methods of WindSaturationCondition.
 *
 * \ingroup Star_group
 */

#include "WindSaturationCondition.h"

namespace Star {
    std::valarray<double> WindSaturationCondition::operator()(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv) const
    {
#ifdef DEBUG
        assert(evol_mode != Core::LOCKED_SURFACE_SPIN);
        if(evol_mode != Core::BINARY) assert(__primary);
#endif 
        unsigned angmom_index = 1 + 2 * __body.number_zones();
        if(evol_mode == Core::BINARY) 
            angmom_index += 2 * __other_body.number_zones()
                            + 
                            (__primary ? 0
                                       : __other_body.number_zones()
                                         -__other_body.number_locked_zones());
        else angmom_index -= 3;
        assert(angmom_index <= derivatives.size());
        double wsurf = __body.spin_frequency(),
               surf_angmom_deriv,
               result = (
                   (std::abs(wsurf) - __saturation_freq)
                   / __saturation_freq
               ),
               wsurf_sign = (wsurf < 0 ? -1.0 : 1.0);
        if(std::isinf(__saturation_freq)) result=-1;
        if(__body.zone(0).locked()) return std::valarray<double>(result, 1);
        surf_angmom_deriv = derivatives[angmom_index];
        stop_deriv.resize(
                1
                ,
                (
                    wsurf_sign * surf_angmom_deriv 
                    - 
                    __body.zone(0).moment_of_inertia(1) * std::abs(wsurf)
                )
                /
                (__body.zone(0).moment_of_inertia() * __saturation_freq)
        );
        return std::valarray<double>(result, 1);
    }

}//End Star namespace.
