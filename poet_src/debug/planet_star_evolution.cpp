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
    MESAInterpolator *interpolator = load_interpolator(
        "stellar_evolution_interpolators/"
        "660499e5-8692-4fd9-a5bd-9397bbe481e3"
    );

    prepare_eccentricity_expansion(
        "eccentricity_expansion_coef_O400.sqlite",
        1e-4,
        true,
        false
    );

    const double Mjup_to_Msun = (Core::AstroConst::jupiter_mass
                                 /
                                 Core::AstroConst::solar_mass),
                 Rjup_to_Msun = (Core::AstroConst::jupiter_radius
                                 /
                                 Core::AstroConst::solar_radius);

    double mstar = 1.0,
           mplanet = 1.0,//100.0 * Mjup_to_Msun,
           rplanet = 1.0 * Rjup_to_Msun,
           a0 = 25.0,
           e0 = 0.8,
           zero = 0.0,
           phase_lag = 3e-8;

    CPlanet *planet = create_planet(mplanet, rplanet);
    configure_planet(planet,
                     5e-3,
                     mstar,
                     a0,
                     e0,
                     &zero,
                     NULL,
                     NULL,
                     false,
                     true,
                     true);

    EvolvingStar *star = create_star(mstar,
                                     0.0,
                                     0.13,
                                     2.78,
                                     1.0e-3,
                                     interpolator);
    select_interpolation_region(star, core_formation_age(star));

    double break_frequency = 2.0 * M_PI / 20;
    double powerlaws[] = {1.0, 0.0};//-3.1};
    set_star_dissipation(star,
                         0,
                         1,
                         0,
                         &break_frequency,
                         NULL,
                         powerlaws,
                         &zero,
                         phase_lag,
                         1.0,
                         10.0);

    DiskBinarySystem *system = create_star_planet_system(
        star,
        planet,
        a0,
        e0,
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
        2.0,
        1e-3,
        1e-4,
        NULL,
        0,
        true,
        0
    );
    int num_steps = num_evolution_steps(solver);
    double *age = new double[num_steps],
           *semimajor = new double[num_steps],
           *eccentricity = new double[num_steps],
           *lconv = new double[num_steps],
           *lrad = new double[num_steps];
    get_star_planet_evolution(
        solver,
        system,
        star,
        planet,
        age,
        semimajor,
        eccentricity,
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
              << std::setw(25) << "e"
              << std::endl;
    for(int i = 0; i < num_steps; ++i)
        std::cout << std::setw(25) << age[i]
                  << std::setw(25) << Core::orbital_angular_velocity(
                      mstar,
                      mplanet,
                      semimajor[i]
                  )
                  << std::setw(25) << (lconv[i]
                                       /
                                       envelope_inertia(star, age[i]))
                  << std::setw(25) << lrad[i] / core_inertia(star, age[i])
                  << std::setw(25) << eccentricity[i]
                  << std::endl;

    destroy_binary(system);
    destroy_planet(planet);
    destroy_star(star);
    destroy_interpolator(interpolator);
    destroy_solver(solver);
    delete[] age;
    delete[] semimajor;
    delete[] eccentricity;
}
