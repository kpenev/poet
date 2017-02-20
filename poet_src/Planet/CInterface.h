/**\file
 *
 * \brief Declare C-style functions for working with LockedPlanet instances.
 *
 * \ingroup Planet_group
 */

#include "LockedPlanet.h"

extern "C" {
    ///Opaque struct to cast LockedPlanet instances to/from.
    struct LockedPlanet;

    ///\brief Create a planet to use in a evolution calculation.
    ///
    ///The result must be de-allocated by the caller.
    LockedPlanet *create_planet(
        ///The mass of the planet to create.
        double mass,

        ///The radius of the planet to create.
        double radius
    );

    ///Destroy a planet previously allocated using create_planet.
    void destroy_planet(LockedPlanet *planet);

}; //End extern "C"
