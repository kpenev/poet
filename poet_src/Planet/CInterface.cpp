/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup Planet_group
 */

#define BUILDING_LIBRARY
#include "CInterface.h"

LockedPlanet *create_planet(double mass, double radius)
{
    return reinterpret_cast<LockedPlanet*>(
        new Planet::LockedPlanet(mass, radius)
    );
}

void destroy_planet(LockedPlanet *planet)
{
    delete reinterpret_cast<Planet::LockedPlanet*>(planet);
}
