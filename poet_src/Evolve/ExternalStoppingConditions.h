/**\file
 *
 * \brief Users can define any stopping condition they wish the evolution to
 * search for in this file.
 *
 * \ingroup OrbitSolver_group
 */


#ifndef __EXTERNAL_STOPPING_CONDITIONS_H
#define __EXTERNAL_STOPPING_CONDITIONS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "StoppingCondition.h"
#include "DissipatingZone.h"

namespace Evolve {

    ///\brief Satisfied when a zone is rotating faster than a threshold.
    ///
    ///\ingroup OrbitSolver_group
    class LIB_LOCAL RotFastCondition : public ExternalStoppingCondition {
    private:
        ///The spin threshold in rad/day.
        double __spin_thres;

        ///Which zone's rotation to monitor.
        DissipatingZone &__zone;
    public:
        ///Create a condition tied to the given threshold in rad/day.
        RotFastCondition(double spin_thres, DissipatingZone &zone) :
            __spin_thres(spin_thres), __zone(zone) {}

        ///\brief Returns the difference between the convective zone spin and
        ///the threshold divided by the latter.
        std::valarray<double> operator()(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv
        ) const;

        ///See StoppingCondition::describe().
        virtual std::string describe(int index = 0) const;

    }; //End RotFastCondition class.

    ///\brief How to construct the external stopping condition(s).
    ///
    ///Leave undefined if no external conditions need to be tracked.
    //#define EXTERNAL_CONDITION RotFastCondition(4.0*M_PI)
    
} //End Evolve namespace.

#endif
