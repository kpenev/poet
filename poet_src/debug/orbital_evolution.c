#include "../Planet/CInterface.h"
#include "../StellarEvolution/CInterface.h"
#include "../Star/CInterface.h"
#include "../Evolve/CInterface.h"
#include <iomanip>
#include <iostream>

int main(int, char **)
{
    read_eccentricity_expansion_coefficients(
        "eccentricity_expansion_coef.txt"
    );

    LockedPlanet *planet = create_planet(1.0, 1.0);
    double zero = 0.0;
    configure_planet(planet,
                     5e-3,
                     1.0,
                     5.0,
                     0.0,
                     &zero,
                     NULL,
                     NULL,
                     false,
                     true,
                     true);

    MESAInterpolator *interpolator = load_interpolator(
        "../stellar_evolution_interpolators/"
        "30b6d5cc-4ae8-43da-899d-740eefd18638"
    );
    EvolvingStar *star = create_star(1.0, 0.0, 0.15, 2.5, 5.0, interpolator);
    select_interpolation_region(star, core_formation_age(star));
    set_dissipation(star, 0, 0, 0, NULL, NULL, &zero, &zero, 3e-7);
    set_dissipation(star, 1, 0, 0, NULL, NULL, &zero, &zero, 0.0);

    StarPlanetSystem *system = create_star_planet_system(
        star,
        planet,
        5.0,
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
    detect_stellar_wind_saturation(star);

    OrbitSolver *solver = evolve_system(
        system,
        10.0,
        1e-3,
        1e-6,
        NULL,
        0
    );
    int num_steps = num_evolution_steps(solver);
    double *age = new double[num_steps],
           *semimajor = new double[num_steps];
    get_evolution(
        solver, system, star,
        age, semimajor, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
    );
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setw(25) << "Age[Gyr]"
              << std::setw(25) << "semimajor[Rsun]"
              << std::endl;
    for(int i = 0; i < num_steps; ++i)
        std::cout << std::setw(25) << age[i]
                  << std::setw(25) << semimajor[i]
                  << std::endl;

    destroy_binary(system);
    destroy_planet(planet);
    destroy_star(star);
    destroy_interpolator(interpolator);
    destroy_solver(solver);
    delete[] age;
    delete[] semimajor;
}
