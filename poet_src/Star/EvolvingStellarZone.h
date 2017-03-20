/**\file
 *
 * \brief Declares a base class for all stellar zones.
 *
 * \ingroup Star_group
 */

#ifndef __EVOLVING_STELLAR_ZONE_H
#define __EVOLVING_STELLAR_ZONE_H

#include "../Evolve/BrokenPowerlawPhaseLagZone.h"
#include "../Core/Functions.h"
#include "../StellarEvolution/EvolvingStellarQuantity.h"

#include <initializer_list>

namespace Star {

    ///\brief Base class for zones of low mass evolving stars.
    class EvolvingStellarZone : virtual public Evolve::BrokenPowerlawPhaseLagZone
    {
    private:
        ///The age for the last configure() call.
        double __current_age;

        ///\brief Pre-computed values and derivatives for quantities which 
        ///only depend on age at the current age.
        mutable std::vector< const Core::FunctionDerivatives* >
            __current_age_quantities;

        ///The quantities describing the stellar zone.
        const std::vector< const StellarEvolution::EvolvingStellarQuantity* >
            __evolving_quantities;

        ///Forgets any previously calculated quantities for the current age.
        void reset_current_quantities();

    public:
        ///Create an evolving stellar zone described by the given quantities.
        ///
        ///WARNING: All quantities are destroyed by the desctructor.
        EvolvingStellarZone(
            ///The quantities describing the zone.
            std::initializer_list<
                const StellarEvolution::EvolvingStellarQuantity*
            > evolving_quantities
        ) : 
            __current_age(Core::NaN),
            __current_age_quantities(evolving_quantities.size(), NULL),
            __evolving_quantities(evolving_quantities)
        {}

        ///\brief Defines the current orbit, triggering re-calculation of all
        ///quantities.
        virtual void configure(
            ///Is this the first time configure() is invoked?
            bool initialize,

            ///The age to set the zone to.
            double age,

            ///The angular velocity of the orbit in rad/day.
            double orbital_frequency,

            ///The eccentricity of the orbit
            double eccentricity,

            ///The absolute value of the angular momentum of the orbit.
            double orbital_angmom,

            ///The angular (momentum/velocity) of the spin of the zone if
            ///the zone is not in a spin-orbit lock (ignored it if is).
            double spin,

            ///The inclination of the zone relative to the orbit.
            double inclination,

            ///The argument of periapsis of the orbit in the equatorial
            ///planet of the zone.
            double periapsis,

            ///Is spin an angular velocity instead of angular momentum?
            bool spin_is_frequency
        );

        ///\brief The current age value of the given quantity (or its 
        ///derivative).
        ///
        ///Computes the value if necessary or retrieves it from
        //__current_age_quantities if already present.
        double current_age_quantity(
            ///The index of the quantity to return from within the list provided
            ///at construction.
            size_t quantity,

            ///The order of the derivitive to return
            unsigned deriv_order=0
        ) const;

        ///\brief The value of the given quantity (or its derivative) at an
        ///arbitrary age.
        ///
        ///Always computes the value/derivatives.
        double any_age_quantity(
            ///The index of the quantity to return from within the list provided
            ///at construction.
            size_t quantity,

            ///The age at which to evaluate the quantity/derivative.
            double age, 

            ///The order of the derivitive to return
            unsigned deriv_order=0
        ) const;

        ///Return the last age with which ::configure() was called.
        double current_age() {return __current_age;}

        ///\brief Change the body as necessary at the given age.
        ///
        ///Handles things like interpolation discontinuities. 
        void reached_critical_age(double age);

        ///\brief The next age when the evolution needs to be stopped for a
        ///change in one of the bodies.
        double next_stop_age() const;

        ///The minimum age at wich zone quantities can be querried.
        double min_interp_age() const;

        ///\brief Prepare the zone quantities for interpolation around the
        ///given age.
        ///
        ///After calling this method, requesting values or derivatives
        ///outside the range of the continuous region containing this age
        ///(see ::discontinuities) fails an assert.
        virtual void select_interpolation_region(double age) const;

        ///Delete any dynamically allocated memory.
        ~EvolvingStellarZone();
    };//End EvolvingStellarZone class.

}//End Star namespace.

#endif
