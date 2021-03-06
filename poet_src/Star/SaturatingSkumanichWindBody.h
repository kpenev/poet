/**\file
 *
 * \brief Decrales a body subject to angular momentum loss 
 * \f$\propto\omega\min(\omega, \omega_{sat})^2\f$.
 *
 * \ingroup Star_group
 */


#ifndef __SATURATING_SKUMANICH_WIND_BODY_H
#define __SATURATING_SKUMANICH_WIND_BODY_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Evolve/DissipatingBody.h"
#include "../Evolve/BinarySystem.h"

namespace Star {

    ///\brief A DissipatingBody which loses angular momentum at a rate
    /// \f$\propto\omega\min(\omega, \omega_{sat})^2\f$
    ///
    ///\ingroup StellarSystem_group
    class LIB_PUBLIC SaturatingSkumanichWindBody : 
        virtual public Evolve::DissipatingBody {
    private:
        ///The strength of the magnetic wind
        double __wind_strength,

            ///The frequency at which the wind loss saturates in rad/day.
            __saturation_freq;

        ///Is the wind currently saturated?
        bool __saturated;

        ///The saturation states recorded by add_to_evolution() so far.
        std::list<bool> __saturation_evolution;

#ifndef NDEBUG
        bool __detected_saturation;
#endif
    public:
        SaturatingSkumanichWindBody(
            ///The strength of the wind.
            double wind_strength,

            ///The frequency at which the wind loss saturates in rad/day.
            double saturation_frequency
        ) :
            __wind_strength(wind_strength), 
            __saturation_freq(saturation_frequency)
#ifndef NDEBUG
            , __detected_saturation(false)
#endif
            {}

        ///See DissipatingBody::angular_momentum_loss().
        double angular_momentum_loss(
            Evolve::Dissipation::QuantityEntry entry = Evolve::Dissipation::NO_DERIV
        ) const;

        ///The frequency at which the wind loss saturates in rad/day.
        double saturation_frequency() {return __saturation_freq;}

        ///Sets the saturation based on the currently configured spin frequency.
        void detect_saturation() 
        {
            __saturated = (std::abs(spin_frequency()) > __saturation_freq);
#ifndef NDEBUG
            __detected_saturation = true;
#endif
        }

        ///Is the wind loss currently saturated?
        bool saturated() {return __saturated;}

        ///Called by the stopping condition monitoring wind saturation.
        void saturation_freq_crossed(
                ///The sign of the rate of change of the spin frequency when it
                ///was equal to the saturation frequency.
                short
#ifndef NDEBUG
                deriv_sign
#endif
                )
        {
            assert(deriv_sign == (__saturated ? -1 : 1));

            __saturated = !__saturated;
        }

        ///Appends the state defined by last configure(), to the evolution.
        virtual void add_to_evolution()
        {
            __saturation_evolution.push_back(__saturated); 
            Evolve::DissipatingBody::add_to_evolution();
        }

        ///Discards the last steps from the evolution.
        virtual void rewind_evolution(
                ///How many steps of evolution to discard.
                unsigned nsteps
        )
        {
            for(unsigned i = 0; i < nsteps; ++i)
                __saturation_evolution.pop_back();
            Evolve::DissipatingBody::rewind_evolution(nsteps);
        }

        ///Discards all evolution
        virtual void reset_evolution()
        {
            __saturation_evolution.clear();
            Evolve::DissipatingBody::reset_evolution();
        }

        ///\brief Conditions detecting the next possible discontinuities in the
        ///evolution due to this body.
        ///
        ///Must be deleted when no longer necessary.
        virtual Evolve::CombinedStoppingCondition *stopping_conditions(
                ///The system being evolved.
                Evolve::BinarySystem &system, 

                ///Is the body the primary in the system.
                bool primary
        );

        ///The tabulated wind saturation states so far.
        const std::list<bool> &wind_saturation_evolution() const
        {return __saturation_evolution;}

        ///Resets its saturation state after a discontinous spin jump.
        void spin_jumped()
        {
            detect_saturation();
            Evolve::DissipatingBody::spin_jumped();
        }

        virtual ~SaturatingSkumanichWindBody() {}
    };//End SaturatingSkumanichWindBody class.

}//End Star namespace.

#endif
