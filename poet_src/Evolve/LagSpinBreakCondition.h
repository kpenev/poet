/**\file
 *
 * \brief Declares a stopping condition monitoring for critical spin
 * frequencies.
 *
 * \ingroup Evolve_group
 */

#ifndef __CRITICAL_SPIN_CONDITION_H
#define __CRITICAL_SPIN_CONDITION_H

#include "../Core/SharedLibraryExportMacros.h"
#include "StoppingCondition.h"
#include "DissipatingBody.h"
#include "../Core/Error.h"
#include <vector>

namespace Evolve {

    class BrokenPowerlawPhaseLagZone;

    ///\brief Satisfied when some zone reaches a critical spin.
    class LIB_LOCAL LagSpinBreakCondition : public StoppingCondition {
    private:
        ///The zone being monitored (for more convenient access).
        BrokenPowerlawPhaseLagZone &__zone;

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

        ///\brief The index of the currently active powerlaw within
        /// __zone.__spin_frequency_powers.
        std::vector<double>::size_type __powerlaw_index;

        ///\brief See num_subcondition().
        unsigned __num_subconditions;

        ///Set the appropriate value of __num_subcondition.
        void set_num_subconditions();

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
        LagSpinBreakCondition(
            ///The zone being monitored.
            BrokenPowerlawPhaseLagZone &zone,

			///The body whose spin to monitor.
			const DissipatingBody &body,

			///The other body in the system.
			const DissipatingBody &other_body,
			
			///Is the body we are monitoring the primary?
			bool primary,

            ///The index (within body) of the zone for which to monitor the
            ///spin.
            unsigned zone_index
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

        ///\brief Adjust the above and below critical frequencies being
        ///monitored.
        ///
        ///See StoppingCondition::reached() for a description of the
        ///arguments.
        virtual void reached(short deriv_sign, unsigned index = 0);

        ///\brief The number of subconditions (1 if all critical spins are
        ///below or above the current spin, 2 otherwise).
        virtual size_t num_subconditions() const
        {return __num_subconditions;}

        ///Define stopping condition type as EXTERNAL.
        virtual StoppingConditionType type(unsigned =0) const
        {return Evolve::EXTERNAL;}

        ///See StoppingCondition::expected_crossing_deriv_sign().
        virtual short expected_crossing_deriv_sign(
            ///Which sub-condition.
            unsigned index = 0
        ) const;

        ///See StoppingCondition::describe().
        virtual std::string describe(int index = -1) const;

    };//End LagSpinBreakCondition class.

}//End Evolve namespace.

#endif
