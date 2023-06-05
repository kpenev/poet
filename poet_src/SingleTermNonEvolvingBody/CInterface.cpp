/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup Planet_group
 */

#define BUILDING_LIBRARY
#include "CInterface.h"
#include "../Evolve/CInterface.h"

CSingleTermNonEvolvingBody *create_single_term_non_evolving_body(double mass,
                                                                 double radius)
{
    return reinterpret_cast<CSingleTermNonEvolvingBody*>(
        new SingleTermNonEvolvingBody::SingleTermNonEvolvingBody(mass, radius)
    );
}

void destroy_single_term_non_evolving_body(CSingleTermNonEvolvingBody *body)
{
    delete reinterpret_cast<SingleTermNonEvolvingBody::SingleTermNonEvolvingBody*>(body);
}

void set_single_term_non_evolving_body_dissipation(
    CSingleTermNonEvolvingBody *body,
    int orbital_frequency_multiplier,
    int spin_frequency_multiplier,
    double phase_lag
)
{
    Evolve::SingleTermZone *zone = &(
        reinterpret_cast<SingleTermNonEvolvingBody::SingleTermNonEvolvingBody*>(
            body
        )->zone()
    );

    set_single_term_zone_dissipation(
        reinterpret_cast<SingleTermZone*>(zone),
        orbital_frequency_multiplier,
        spin_frequency_multiplier,
        phase_lag
    );
}

double get_single_term_non_evolving_body_inertia(
    CSingleTermNonEvolvingBody *body
)
{
    return reinterpret_cast<SingleTermNonEvolvingBody::SingleTermNonEvolvingBody*>(
        body
    )->zone().moment_of_inertia();
}
