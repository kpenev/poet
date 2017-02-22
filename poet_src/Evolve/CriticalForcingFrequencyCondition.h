/**\file
 *
 * \brief Declares a stopping condition satisfied when a forcing frequency
 * reaches a critical value.
 *
 * \ingroup Evolve_group
 */

#ifndef __CRITICAL_FORCING_FREQUENCY_CONDITION_H
#define __CRITICAL_FORCING_FREQUENCY_CONDITION_H

#include "StoppingCondition.h"
#include "DissipatingBody.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/Error.h"
#include <vector>

namespace Evolve {
    class DissipatingZone;

    ///\brief satisfied when a forcing frequency reaches a critical value.
    class CriticalForcingFrequencyCondition : public StoppingCondition {
    private:
        ///The critical forcing frequencies to watch for in rad/day.
        std::vector<double> __critical_frequencies;

        int 
            ///The orbital frequency multiplier of the forcing term being
            ///monitored.
            __orbital_frequency_multiplier,

            ///The spin frequency multiplier of the forcing term being
            ///monitored.
            __spin_frequency_multiplier;

        ///The zone being monitored (for more convenient access).
        const DissipatingZone &__zone;

        std::vector<double>::const_iterator
            ///\brief The __critical_frequencies entry immediately above the
            ///current spin of the monitored zone. If all entries are below,
            ///the value is __critical_spins.end()
            __critical_above_iter,

            ///\brief The __critical_frequencies entry immediately below the
            ///current spin of the monitored zone. If all entries are above,
            ///the value is __critical_spins.end()
            __critical_below_iter;

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

        ///\brief See num_subcondition().
        unsigned __num_subconditions;

        ///Set the appropriate value of __num_subcondition.
        void set_num_subconditions();
    public:
        ///\brief Monitor a single forcing term of a single zone for a number
        ///of critical forcing frequencies.
        CriticalForcingFrequencyCondition(
            ///The body whose spin to monitor.
			const DissipatingBody &body,

			///The other body in the system.
			const DissipatingBody &other_body,
			
			///Is the body we are monitoring the primary?
			bool primary,

            ///The index (within body) of the zone for which to monitor the
            ///spin.
            unsigned zone_index,

            ///The orbital frequency multiplier of the forcing term being
            ///monitored.
            int orbital_frequency_multiplier,

            ///The spin frequency multiplier of the forcing term being
            ///monitored.
            int spin_frequency_multiplier,

            ///The critical spin frequency to watch for.
            std::vector<double> critical_frequencies,

            ///The current orbital frequency for which to initialize this
            ///condition.
            double orbital_frequency
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
    };//End CriticalForcingFrequencyCondition class.
} //End Evolve namespace.

#endif
