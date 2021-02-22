#include "../Planet/CInterface.h"
#include "../StellarEvolution/CInterface.h"
#include "../Star/CInterface.h"
#include "../Evolve/CInterface.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/AstronomicalConstants.h"
#include <iomanip>
#include <iostream>
#include <valarray>

int main(int, char **)
{
    read_eccentricity_expansion_coefficients(
        "eccentricity_expansion_coef.txt"
    );

    double mstar = 1.0, mplanet = 1.0, a0 = 10.0, zero = 0.0;

    CPlanet *planet = create_planet(1.0, 1.0);
    configure_planet(planet,
                     5e-3,
                     mstar,
                     a0,
                     0.0,
                     &zero,
                     NULL,
                     NULL,
                     false,
                     true,
                     true);

    MESAInterpolator *interpolator = load_interpolator(
        "stellar_evolution_interpolators/"
        "90af7144-f918-4a1c-95a2-0b086a80d0a2"
    );
    EvolvingStar *star = create_star(mstar,
                                     0.0,
                                     0.13,
                                     2.78,
                                     1.0e-3,
                                     interpolator);
    select_interpolation_region(star, core_formation_age(star));

    double break_frequency = 2.0 * M_PI / 4.33;
    double powerlaws[] = {0.0, -3.1};
    set_star_dissipation(star, 0, 1, 0, &break_frequency, NULL, powerlaws, &zero, 3e-7, 1.0, 0.0);
    set_star_dissipation(star, 1, 0, 0, NULL, NULL, &zero, &zero, 0.0, 1.0, 0.0);

    DiskBinarySystem *system = create_star_planet_system(
        star,
        planet,
        a0,
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
        0,
        true,
        0
    );
    int num_steps = num_evolution_steps(solver);
    double *age = new double[num_steps],
           *semimajor = new double[num_steps],
           *lconv = new double[num_steps],
           *lrad = new double[num_steps];
    get_star_planet_evolution(
        solver,
        system,
        star,
        planet,
        age,
        semimajor,
        NULL, //eeccentricity
        NULL, //envelope inclination
        NULL, //core inclination
        NULL, //envelope periapsis
        NULL, //core periapsis
        lconv,//envelope angmom
        lrad, //core angmom
        NULL, //planet inclination
        NULL, //planet periapsis
        NULL, //planet angular momentum
        NULL, //evolution mode
        NULL, //wind saturation
        NULL, //semimajor rate
        NULL, //eeccentricity rate
        NULL, //envelope inclination rate
        NULL, //core inclination rate
        NULL, //envelope periapsis rate
        NULL, //core periapsis rate
        NULL, //envelope angmom rate
        NULL, //core angmom rate
        NULL, //planet inclination rate
        NULL, //planet periapsis rate
        NULL  //planet angular momentum rate
    );
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setw(25) << "Age[Gyr]"
              << std::setw(25) << "worb[rad/day]"
              << std::setw(25) << "wconv[rad/day]"
              << std::setw(25) << "wrad[rad/day]"
              << std::endl;
    for(int i = 0; i < num_steps; ++i)
        std::cout << std::setw(25) << age[i]
                  << std::setw(25) << Core::orbital_angular_velocity(
                      mstar,
                      (
                          mplanet
                          *
                          Core::AstroConst::jupiter_mass
                          /
                          Core::AstroConst::solar_mass
                      ),
                      semimajor[i]
                  )
                  << std::setw(25) << (lconv[i]
                                       /
                                       envelope_inertia(star, age[i]))
                  << std::setw(25) << lrad[i] / core_inertia(star, age[i])
                  << std::endl;

    destroy_binary(system);
    destroy_planet(planet);
    destroy_star(star);
    destroy_interpolator(interpolator);
    destroy_solver(solver);
    delete[] age;
    delete[] semimajor;
}
