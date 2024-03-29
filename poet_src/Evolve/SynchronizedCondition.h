/**\file
 *
 * \brief Declares a stopping condition monitoring spin-orbit
 * synchronization.
 *
 * \ingroup Evolve_group
 */

#ifndef __SYNCHRONIZED_CONDITION_H
#define __SYNCHRONIZED_CONDITION_H

#include "../Core/SharedLibraryExportMacros.h"
#include "StoppingCondition.h"
#include "SpinOrbitLockInfo.h"

namespace Evolve {

    class DissipatingZone;
    class BinarySystem;

    ///\brief Satisfied when some multiples of the orbit and stellar rotation are
    ///synchronized.
    class LIB_LOCAL SynchronizedCondition : public StoppingCondition{
    private:
        ///The mutiplier in front of the orbital frequency in the lock.
        int __orbital_freq_mult,

            ///The multiplier in front of the spin frequency in the lock.
            __spin_freq_mult;

        ///Which body's spin is checked for locking.
        bool __primary;

        ///Which zone is checked for locking.
        unsigned __zone_index;

        ///The zone whose spin is monitored.
        const DissipatingZone &__zone;

        ///The binary system this locking condition is attached to.
        BinarySystem &__system;

    public:
        ///Create the synchronization condition for the given object.
        SynchronizedCondition(
            ///The mutiplier in front of the orbital frequency in the lock.
            int orbital_freq_mult,

            ///The multiplier in front of the spin frequency in the lock.
            int spin_freq_mult,

            ///The sign the first derivative should have if a crossing
            ///occurs.
            short deriv_sign,

            ///Which body's spin is checked for locking.
            bool primary,

            ///Which zone's spin is checked for locking.
            unsigned zone_index,

            ///The binary system this locking condition is attached to
            BinarySystem &system
        );

        ///Create the synchronization condition for the given object from a lock
        SynchronizedCondition(
            ///The lock to base the condition on.
            const SpinOrbitLockInfo &monitored_lock,

            ///Which body's spin is checked for locking.
            bool primary,

            ///Which zone's spin is checked for locking.
            unsigned zone_index,

            ///The binary system this locking condition is attached to
            BinarySystem &system
        ) :
            SynchronizedCondition(
                monitored_lock.orbital_frequency_multiplier(),
                monitored_lock.spin_frequency_multiplier(),
                monitored_lock.lock_direction(),
                primary,
                zone_index,
                system
            )
        {}


        ///\brief Return the difference between the orbital and multiplier
        ///scaled stellar spin angular velocities divided by the orbital angular
        ///velocity.
        ///
        ///See StoppingCondition::operator()() for a description of the
        ///arguments.
        std::valarray<double> operator()(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv
        ) const;

        ///Identify this as a SYNCHRONIZED condition.
        StoppingConditionType type(unsigned =0) const {return SYNCHRONIZED;}

        ///The multiplier in front of the orbital frequency in the lock.
        int orbital_frequency_multiplier() const {return __orbital_freq_mult;}

        ///The multiplier in front of the spin frequency in the lock.
        int spin_frequency_multiplier() const {return __spin_freq_mult;}

        ///Which body's spin is checked for locking.
    //	short body_index() const {return __body_index;}

        ///See StoppingCondition::reached().
        void reached(short deriv_sign, unsigned index=0);

        ///See StoppingCondition::describe().
        virtual std::string describe(int index = -1) const;

    };//End SynchronizedCondition class.

}//End Evolve namespace.

#endif
