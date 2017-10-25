/**\file
 *
 * \brief Declaration of a stopping condition monitoring for the surface spin
 * of a star crossing the wind saturation threshold.
 *
 * \ingroup Star_group
 */

#ifndef __WIND_SATURATION_CONDITION_H
#define __WIND_SATURATION_CONDITION_H

#include "../Evolve/StoppingCondition.h"
#include "SaturatingSkumanichWindBody.h"

namespace Star {
    ///\brief Satisfied when the surface zone of a body is spinning at exactly
    ///the wind saturation frequency.
    ///
    ///\ingroup OrbitSolver_group
    class WindSaturationCondition : public Evolve::StoppingCondition {
    private:
        ///The frequency at which the wind saturates.
        double __saturation_freq;

        ///The body this condition is monitoring.
        SaturatingSkumanichWindBody &__body;

        ///The primary body in the system if not __body.
        const Evolve::DissipatingBody &__other_body;

        ///Is __body we are the primary in the system?
        bool __primary;
    public:
        WindSaturationCondition(
            ///The body whose wind saturation we are monitoring.
            SaturatingSkumanichWindBody &body,

            ///The other body in the system.
            const Evolve::DissipatingBody &other_body,

            ///Is the body we are monitoring the primary?
            bool primary,

            ///Is the wind currently saturated or not?
            bool saturated
        ) :
            StoppingCondition((saturated ? -1 : 1)),
            __saturation_freq(body.saturation_frequency()),
            __body(body),
            __other_body(other_body),
            __primary(primary)
        {}

        ///\brief The difference between the convective and wind saturation
        ///angular velocities divided by the latter.
        std::valarray<double> operator()(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv
        ) const;

        ///Identify this as a WIND_SATURATION condition.
        Evolve::StoppingConditionType type(unsigned =0) const
        {return Evolve::WIND_SATURATION;}

        ///See StoppingCondition::reached().
        void reached(short deriv_sign, unsigned index=0)
        {
            Evolve::StoppingCondition::reached(deriv_sign, index);
            __body.saturation_freq_crossed(deriv_sign);
        }

        ///See StoppingCondition::describe().
        virtual std::string describe(int index = -1) const;

    };//End WindSaturationCondition class.

}//End Star namespace.
#endif
