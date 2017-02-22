/**\file
 *
 * \brief Declaration of the class representing stellar cores.
 *
 * \ingroup Star_group
 */

#ifndef __EVOLVING_STELLAR_CORE_H
#define __EVOLVING_STELLAR_CORE_H

#include "EvolvingStellarZone.h"

#include "../Core/Functions.h"
#include "../StellarEvolution/EvolvingStellarQuantity.h"

namespace Star {

    ///\brief Radiative core for low mass evolving stars.
    class EvolvingStellarCore : public EvolvingStellarZone {
    private:
        ///The age at which the core first forms.
        double __formation_age;

        ///\brief Identifiers for the various age dependent values which are only
        ///computed once per fixed age.
        ///
        ///All quantities also refer to their derivatives.
        enum QuantityIndex {
            MASS,   ///< The outer mass boundary of the zone.
            RADIUS, ///< The outer radius boundary of the zone.
            INERTIA,///< The moment of inertia
        };

    public:
        ///Create a stellar core with the specified properties.
        EvolvingStellarCore(
                ///The age at which the core forms.
                double formation_age=Core::Inf, 

                ///The mass of the core.
                const StellarEvolution::EvolvingStellarQuantity *mass=NULL,

                ///The radius of the core.
                const StellarEvolution::EvolvingStellarQuantity *radius=NULL,

                ///The moment of inertia of the zone.
                const StellarEvolution::EvolvingStellarQuantity 
                *moment_of_inertia=NULL
        ) :
            EvolvingStellarZone({mass, radius, moment_of_inertia}),
            __formation_age(formation_age)
        {}

        ///See DissipatingZone::moment_of_inertia(int).
        double moment_of_inertia(int deriv_order=0) const
        {return current_age_quantity(INERTIA, deriv_order);}

        ///See DissipatingZone::moment_of_inertia(double, int).
        double moment_of_inertia(double age, int deriv_order=0) const
        {return any_age_quantity(INERTIA, age, deriv_order);}

        ///See DissipatingZone::outer_radius(int).
        double outer_radius(int deriv_order=0) const
        {return current_age_quantity(RADIUS, deriv_order);}

        ///See DissipatingZone::outer_radius(double, int).
        double outer_radius(double age, int deriv_order=0) const
        {return any_age_quantity(RADIUS, age, deriv_order);}

        ///See DissipatingZone::outer_mass(int).
        double outer_mass(int deriv_order=0) const
        {return current_age_quantity(MASS, deriv_order);}

        ///See DissipatingZone::outer_mass(double, int).
        double outer_mass(double age, int deriv_order=0) const
        {return any_age_quantity(MASS, age, deriv_order);}

        ///No dissipation in the core.
        bool dissipative() const {return false;} 

        ///The age at which the core forms in Gyr.
        double formation_age() const {return __formation_age;}
    }; //End EvolvingStellarCore class.

}//End Star namespace.

#endif
