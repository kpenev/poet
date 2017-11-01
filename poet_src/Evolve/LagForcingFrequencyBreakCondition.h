/**\file
 *
 * \brief Declares a stopping condition satisfied when a forcing frequency
 * reaches a critical value.
 *
 * \ingroup Evolve_group
 */

#ifndef __CRITICAL_FORCING_FREQUENCY_CONDITION_H
#define __CRITICAL_FORCING_FREQUENCY_CONDITION_H

#include "../Core/SharedLibraryExportMacros.h"
#include "StoppingCondition.h"
#include "DissipatingBody.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/Error.h"
#include <vector>

namespace Evolve {
    class BrokenPowerlawPhaseLagZone;

    ///\brief satisfied when a forcing frequency reaches a critical value.
    class LIB_LOCAL LagForcingFrequencyBreakCondition : 
        public StoppingCondition {
        private:
            int 
                ///The orbital frequency multiplier of the forcing term being
                ///monitored.
                __orbital_frequency_multiplier,

                ///The spin frequency multiplier of the forcing term being
                ///monitored.
                __spin_frequency_multiplier;

            ///The zone being monitored (for more convenient access).
            BrokenPowerlawPhaseLagZone &__zone;

            const DissipatingBody 
                ///The body this condition is monitoring.
                &__body, 

                ///The other body in the system.
                &__other_body;

            ///\brief The index of the monitored tidal term within
            /// __zone.__tidal_indices.
            std::vector< std::vector<double>::size_type >::size_type
                __term_index;

            ///\brief The index of the currently active powerlaw within
            /// __zone.__tidal_frequency_powers.
            std::vector<double>::size_type __powerlaw_index;

            ///\brief See num_subcondition().
            unsigned __num_subconditions;

            ///Set the appropriate value of __num_subcondition.
            void set_num_subconditions();
        public:
            ///\brief Monitor a single forcing term of a single zone for a
            ///number of critical forcing frequencies.
            LagForcingFrequencyBreakCondition(
                ///The zone being monitored.
                BrokenPowerlawPhaseLagZone &zone,

                ///The body whose spin to monitor.
                const DissipatingBody &body,

                ///The other body in the system.
                const DissipatingBody &other_body,

                ///The orbital frequency multiplier of the forcing term being
                ///monitored.
                int orbital_frequency_multiplier,

                ///The spin frequency multiplier of the forcing term being
                ///monitored.
                int spin_frequency_multiplier
            );

            ///\brief Return the differences between the current forcing
            ///frequency of the zone monitored and the two closest critical
            ///frequencies divided by the corresponding critical frequency.
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

            ~LagForcingFrequencyBreakCondition()
            {std::cerr << "Destroying: " << describe() << std::endl;}
        };//End LagForcingBreakFrequencyCondition class.

} //End Evolve namespace.

#endif
