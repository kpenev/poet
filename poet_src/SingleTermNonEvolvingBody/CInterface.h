/**\file
 *
 * \brief Declare C-style functions for working with SingleTermNonEvolvingBody
 * instances.
 *
 * \ingroup Planet_group
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "SingleTermNonEvolvingBody.h"

extern "C" {
    ///\brief Opaque struct to cast
    //SingleTermNonEvolvingBody::SingleTermNonEvolvingBody instances to/from.
    struct CSingleTermNonEvolvingBody;

    ///\brief Create a an object with only one tidal term contributing to the
    ///dissipation and with non-elolving structure to use in a evolution calculation.
    ///
    ///The result must be de-allocated by the caller.
    LIB_PUBLIC CSingleTermNonEvolvingBody *create_single_term_non_evolving_body(
        ///The mass of the object to create.
        double mass,

        ///The radius of the object to create.
        double radius
    );

    ///\brief Destroy an object previously allocated using
    ///create_single_term_non_evolving_body
    LIB_PUBLIC void destroy_single_term_non_evolving_body(
        CSingleTermNonEvolvingBody *body
    );

    ///\brief Set the dissipative properties of the single tidal term non evolving
    ///object
    LIB_PUBLIC void set_single_term_non_evolving_body_dissipation(
        ///The body to set the dissipation for.
        CSingleTermNonEvolvingBody *body,

        ///The multiplier of the orbital frequency in the expression for
        ///the forcing frequency for the only dissipative term.
        int orbital_frequency_multiplier,

        ///The multiplier of the spin frequency in the expression for
        ///the forcing frequency for the only dissipative term.
        int spin_frequency_multiplier,

        ///The phase lag to assume for the only dissipative term.
        double phase_lag
    );

    ///Return the moment of inertia of the given non evolving object.
    LIB_PUBLIC double get_single_term_non_evolving_body_inertia(
        ///The object to get the moment of inertia for
        CSingleTermNonEvolvingBody *body
    );

} //End extern "C"
