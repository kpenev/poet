#include "../Planet/CInterface.h"
#include "../StellarEvolution/CInterface.h"
#include "../Star/CInterface.h"
#include "../Evolve/CInterface.h"

int main(int, char **)
{
    LockedPlanet *planet = create_planet(1.0, 1.0);
    double zero = 0.0;
    configure_body(reinterpret_cast<DissipatingBody*>(planet),
                   5e-3,
                   1.0,
                   4.0,
                   0.0,
                   &zero,
                   NULL,
                   NULL,
                   false,
                   true,
                   true);

    MESAInterpolator *interpolator = load_interpolator(
        "../stellar_evolution_interpolators/"
        "93878c27-baa3-4cbe-8f78-0ad2a21e121f"
    );
    EvolvingStar *star = create_star(1.0, 0.0, 0.15, 2.5, 5.0, interpolator);
    set_dissipation(star, 0, 0, 0, NULL, NULL, &zero, &zero, 1e-5);
    set_dissipation(star, 0, 0, 0, NULL, NULL, &zero, &zero, 0.0);

    DiskPlanetSystem *system = create_disk_planet_system(
        reinterpret_cast<DissipatingBody*>(star),
        reinterpret_cast<DissipatingBody*>(planet),
        4.0,
        0.0,
        0.0,
        1.0,
        5e-3,
        0.0
    );
    configure_system(system,
                     core_formation_age(star),
                     NaN,
                     NaN,
                     &zero,
                     NULL,
                     NULL,
                     LOCKED_SURFACE_SPIN_EVOL_MODE);

    OrbitSolver *solver = evolve_system(
        system,
        10.0,
        1e-3,
        1e-6,
        NULL,
        0
    );

    destroy_binary(system);
    destroy_planet(planet);
    destroy_interpolator(interpolator);
}
