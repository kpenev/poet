#ifndef __SECONDARY_DEATH_CONDITION_H
#define __SECONDARY_DEATH_CONDITION_H

#include "StoppingCondition.h"

namespace Evolve {

    class BinarySystem;

    ///\brief Satisfied when the planet enters below either the roche sphere or
    ///the stellar photosphere.
    ///
    ///\ingroup OrbitSolver_group
    class SecondaryDeathCondition : public StoppingCondition {
    private:
        ///The system this condition is attached to.
        BinarySystem &__system;
    public:
        ///\brief Create a condition watching for the death of the secondary body
        ///in a system due to tidal disruption of engulfment.
        SecondaryDeathCondition(BinarySystem &system) :
            StoppingCondition(-1),
            __system(system) 
        {}

        ///\brief The difference between the semimajor axis and the larger of
        ///the roche radius and the stellar radius divided by the latter.
        ///
        ///See StoppingCondition::operator()() for a description of the
        ///arguments.
        ///
        ///The evolution mode must be FAST_PLANET, LOCKED_TO_PLANET or
        ///SLOW_PLANET.
        std::valarray<double> operator()(
                Core::EvolModeType evol_mode,
                const std::valarray<double> &orbit,
                const std::valarray<double> &derivatives,
                std::valarray<double> &stop_deriv) const;

        ///Identify this as a PLANET_DEATH condition.
        StoppingConditionType type(unsigned =0) const {return PLANET_DEATH;}

        ///See StoppingCondition::reached().
        void reached(short deriv_sign, unsigned index=0);

        ///See StoppingCondition::describe().
        virtual std::string describe(int index = -1) const;

    }; //End SecondaryDeathCondition class.

} //End Evolve namespace.

#endif
