/**\file
 *
 * \brief Declare C-style functions for working with LockedPlanet instances.
 *
 * \ingroup Planet_group
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "Planet.h"

extern "C" {
    ///Opaque struct to cast Planet::Planet instances to/from.
    struct CPlanet;

    ///\brief Create a planet to use in a evolution calculation.
    ///
    ///The result must be de-allocated by the caller.
    LIB_PUBLIC CPlanet *create_planet(
        ///The mass of the planet to create.
        double mass,

        ///The radius of the planet to create.
        double radius
    );

    ///Destroy a planet previously allocated using create_planet.
    LIB_PUBLIC void destroy_planet(CPlanet *planet);

    ///Set the dissipative properties of the planet
    LIB_PUBLIC void set_planet_dissipation(
        ///The planet to set the dissipation for.
        CPlanet *planet,

        ///See same name argument to set_zone_dissipation()
        unsigned num_tidal_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        unsigned num_spin_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        unsigned num_age_breaks,

        ///See same name argument to set_zone_dissipation()
        double *tidal_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        double *spin_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        double *tidal_frequency_powers,

        ///See same name argument to set_zone_dissipation()
        double *spin_frequency_powers,

        ///See same name argument to set_zone_dissipation()
        double *age_breaks,

        ///See same name argument to set_zone_dissipation()
        double *reference_phase_lags,

        ///See same name argument to set_zone_dissipation()
        double inertial_mode_enhancement,

        ///See same name argument to set_zone_dissipation()
        double inertial_mode_sharpness
    );

    ///Return the moment of inertia of the given planet.
    LIB_PUBLIC double get_planet_inertia(
        ///The planet to get the moment of inertia for
        CPlanet *planet
    );

} //End extern "C"
