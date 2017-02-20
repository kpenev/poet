/**\file
 *
 * \brief Declares a stopping condition monitoring for critical spin
 * frequencies.
 *
 * \ingroup Evolve_group
 */

#ifndef __CRITICAL_SPIN_CONDITION_H
#define __CRITICAL_SPIN_CONDITION_H

#include "StoppingCondition.h"
#include "../Core/Error.h"

namespace Evolve {

    class DissipatingZone;

    ///\brief Satisfied when some zone reaches a critical spin.
    class CriticalSpinCondition : public StoppingCondition {
    private:
        ///The critical spin frequencies to watch for in rad/day.
        std::vector<double> __critical_spins;

        ///The zone being monitored (for more convenient access).
        const DissipatingZone &__zone;

        std::vector<double>::const_iterator
            ///\brief The __critical_spins entry immediately below the
            ///current spin of the monitored zone. If all entries are above,
            ///the value is __critical_spins.end()
            __critical_below_iter,

            ///\brief The __critical_spins entry immediately above the
            ///current spin of the monitored zone. If all entries are below,
            ///the value is __critical_spins.end()
            __critical_above_iter;

        const DissipatingBody 
            ///The body this condition is monitoring.
            &__body, 

            ///The primary body in the system if not __body.
            &__other_body;

        ///Is __body we are the primary in the system?
        bool __primary;

        ///\brief The index (within __body) of the zone whose spin is being
        ///monitored.
        unsigned __zone_index;

        ///Calculate the derivative of the stopping condition.
        double derivative(
            ///The derivative of the anguler momentum of the monitored zone
            ///in \f$M_\odot\,R_\odot^2\,rad\,day^{-1}\,Gyr^{-1}.
            double surf_angmom_deriv,

            ///The critical spin frequency in rad/day.
            double wcritical
        ) const;

        ///\brief Fill the derivatives argument of operator() for a locked
        ///zone.
        ///
        ///See StoppingCondition::operator()() for a description of the
        ///arguments.
        void fill_locked_derivs(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv
        ) const;

        ///\brief Fill the derivatives argument of operator() for an unlocked
        ///zone.
        ///
        ///See StoppingCondition::operator()() for a description of the
        ///arguments.
        void fill_unlocked_derivs(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv
        ) const;

    public:
        ///Create a critical spin condition for the given zone.
        CriticalSpinCondition(
			///The body whose spin to monitor.
			const DissipatingBody &body,

			///The other body in the system.
			const DissipatingBody &other_body,
			
			///Is the body we are monitoring the primary?
			bool primary

            ///The index (within body) of the zone for which to monitor the
            ///spin.
            unsigned zone_index,

            ///The critical spin frequency to watch for.
            std::vector<double> critical_spins
        );

        ///\brief Return the differences between the current spin of the zone
        ///monitored and the two closest critical frequencies divided by the
        ///corresponding critical frequency.
        ///
        ///See StoppingCondition::operator()() for a description of the
        ///arguments.
        ///
        ///All evolution modes allowed.
        std::valarray<double> operator()(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv
        ) const;

    };//End CriticalSpinCondition class.

}//End Evolve namespace.

#endif
