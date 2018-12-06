/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup Planet_group
 */

#define BUILDING_LIBRARY
#include "CInterface.h"

CPlanet *create_planet(double mass, double radius)
{
    return reinterpret_cast<CPlanet*>(
        new Planet::Planet(mass, radius)
    );
}

void destroy_planet(CPlanet *planet)
{
    delete reinterpret_cast<Planet::Planet*>(planet);
}

LIB_PUBLIC void set_planet_dissipation(CPlanet *planet,
                                       unsigned num_tidal_frequency_breaks,
                                       unsigned num_spin_frequency_breaks,
                                       double *tidal_frequency_breaks,
                                       double *spin_frequency_breaks,
                                       double *tidal_frequency_powers,
                                       double *spin_frequency_powers,
                                       double reference_phase_lag)
{
    Evolve::BrokenPowerlawPhaseLagZone *zone = &(
        reinterpret_cast<Planet::Planet*>(planet).zone()
    );

    set_zone_dissipation(
        reinterpret_cast<BrokenPowerlawPhaseLagZone*>(zone),
        num_tidal_frequency_breaks,
        num_spin_frequency_breaks,
        tidal_frequency_breaks,
        spin_frequency_breaks,
        tidal_frequency_powers,
        spin_frequency_powers,
        reference_phase_lag
    );
}
