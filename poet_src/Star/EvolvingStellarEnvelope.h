#ifndef __EVOLVING_STELLAR_ENVELOPE_H
#define __EVOLVING_STELLAR_ENVELOPE_H

/**\file
 *
 * \brief Declares a class representing convective zones in low mass evolving
 * stars or the whole star for high mass stars.
 *
 * \ingroup Star_group
 */

#include "EvolvingStellarZone.h"

namespace Star {

    ///\brief Surface convective zone for low mass evolving stars or the entire
    ///star for high mass stars.
    class EvolvingStellarEnvelope : public EvolvingStellarZone {
    private:
        ///The age for the last configure() call.
        double __stellar_mass;

        ///\brief Identifiers for the various age dependent values which are only
        ///computed once per fixed age.
        ///
        ///All quantities also refer to their derivatives.
        enum CurrentAgeQuantities {
            RADIUS, ///< The radius
            INERTIA,///< The moment of inertia
        };

    public:
        EvolvingStellarEnvelope(
                ///The mass of the star this zone is part of.
                double stellar_mass,

                ///The radius of the star.
                const StellarEvolution::EvolvingStellarQuantity
                *outer_radius=NULL,

                ///The moment of inertia of the zone.
                const StellarEvolution::EvolvingStellarQuantity
                *moment_of_inertia=NULL
        ) :
            EvolvingStellarZone({outer_radius, moment_of_inertia}),
            __stellar_mass(stellar_mass)
        {
#ifndef NDEBUG
            std::cerr << "Created stellar envelope." << std::endl;
#endif
        }

        ///See DissipatingZone::moment_of_inertia(double, int)
        double moment_of_inertia(double age, int deriv_order=0) const
        {return any_age_quantity(INERTIA, age, deriv_order);}

        ///See DissipatingZone::moment_of_inertia(int).
        double moment_of_inertia(int deriv_order=0) const
        {return current_age_quantity(INERTIA, deriv_order);}

        ///See DissipatingZone::outer_radius(double, int).
        double outer_radius(double age, int deriv_order=0) const
        {return any_age_quantity(RADIUS, age, deriv_order);}

        ///See DissipatingZone::outer_radius(int).
        double outer_radius(int deriv_order=0) const
        {return current_age_quantity(RADIUS, deriv_order);}

        ///See DissipatingZone::outer_mass(double, int).
        double outer_mass(double, int deriv_order=0) const
        {return (deriv_order ? 0 : __stellar_mass);}

        ///See DissipatingZone::outer_mass.
        double outer_mass(int deriv_order=0) const
        {return (deriv_order ? 0 : __stellar_mass);}

    };//End EvolvingStellarEnvelope class.

}//End Star namespace.

#endif
