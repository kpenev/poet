/**\file
 *
 * \brief Definitions of some of the methods of StaturatingSkumanichWindBody.
 *
 * \ingroup Star_group
 */

#include "SaturatingSkumanichWindBody.h"
#include "WindSaturationCondition.h"

namespace Star {

    double SaturatingSkumanichWindBody::angular_momentum_loss(
            Evolve::Dissipation::Derivative deriv
    ) const
    {
        double result = (__wind_strength*std::sqrt(radius() / mass())
                         *
                         (__saturated
                          ? spin_frequency()*std::pow(__saturation_freq, 2)
                          : std::pow(spin_frequency(), 3)));
        double freq_power = (__saturated ? 1.0 : 3.0);
        switch(deriv) {
            case Evolve::Dissipation::NO_DERIV :
                return result;
            case Evolve::Dissipation::SPIN_FREQUENCY :
                return freq_power * result / spin_frequency();
            case Evolve::Dissipation::RADIUS :
                return result / (2.0 * radius());
            case Evolve::Dissipation::MOMENT_OF_INERTIA :
                return -freq_power * result / zone(0).moment_of_inertia();
            case Evolve::Dissipation::SPIN_ANGMOM : 
                return freq_power * result / zone(0).angular_momentum();
            default :
                return 0;
        }
    }

    Evolve::CombinedStoppingCondition *
        SaturatingSkumanichWindBody::stopping_conditions(
            Evolve::BinarySystem &system,
            bool primary
        )
    {
#ifdef DEBUG
        if(primary) assert(this == &(system.primary()));
        else assert(this == &(system.secondary()));
#endif
        Evolve::CombinedStoppingCondition *result = 
            new Evolve::CombinedStoppingCondition();
        if(system.evolution_mode() != Core::LOCKED_SURFACE_SPIN)
            (*result) |= new WindSaturationCondition(
                *this,
                (primary ? system.secondary() : system.primary()),
                primary
            );
        (*result) |= Evolve::DissipatingBody::stopping_conditions(system,
                                                                  primary);
        return result;
    }

}//End Star namespace.
