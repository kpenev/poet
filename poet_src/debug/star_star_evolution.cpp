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

    MESAInterpolator *interpolator = load_interpolator(
        "../stellar_evolution_interpolators/"
//        "93878c27-baa3-4cbe-8f78-0ad2a21e121f"
        "30b6d5cc-4ae8-43da-899d-740eefd18638"
    );


    const double mprimary = 1.0,
                 msecondary = 1.0,
                 a0 = 10.0,
                 tdisk = 5e-3;
    double zero = 0.0;
    double phase_lag = lag_from_lgQ(6.0);
    double initial_secondary_angmom[] = {0.79097531984229608,
                                         7.6031508107941889e-02};

    EvolvingStar *primary = create_star(mprimary,   //mass
                                        0.0,        //[Fe/H]
                                        0.17,       //Kw
                                        2.78,       //wsat
                                        5.0e-3,     //tcoup
                                        interpolator);
    select_interpolation_region(primary, core_formation_age(primary));
    set_dissipation(primary, 0, 0, 0, NULL, NULL, &zero, &zero, phase_lag);
    set_dissipation(primary, 1, 0, 0, NULL, NULL, &zero, &zero, 0.0);

    EvolvingStar *secondary = create_star(msecondary,   //mass
                                          0.0,          //[Fe/H]
                                          0.17,   //Kw
                                          2.78,         //wsat
                                          5.0e-3,       //tcoup
                                          interpolator);
    select_interpolation_region(secondary, tdisk);
    set_dissipation(secondary, 0, 0, 0, NULL, NULL, &zero, &zero, 0.0);
    set_dissipation(secondary, 1, 0, 0, NULL, NULL, &zero, &zero, 0.0);

    configure_star(secondary,
                   tdisk,                       //formation age
                   mprimary,                    //companion mass
                   a0,                          //formation semimajor
                   0.0,                         //formation eccentricity
                   initial_secondary_angmom,    //spin angular momentum
                   &zero,                       //inclination
                   &zero,                       //periapsis
                   false,                       //locked surface?
                   true,                        //zero outer inclination?
                   true);                       //zero outer periapsis?
    detect_stellar_wind_saturation(secondary);

    DiskBinarySystem *system = create_star_star_system(
        primary,            //primary
        secondary,          //secondary
        a0,                 //initial semimajor
        0.0,                //initial eccentricity
        0.0,                //initial inclination
        2.0 * M_PI / 3.0,   //disk lock frequency
        tdisk,              //disk dissipation age
        tdisk               //secondary formation age
    );
    configure_system(system,
                     core_formation_age(primary),
                     NaN,                       //semimajor
                     NaN,                       //eccentricity
                     &zero,                     //spin angmom
                     NULL,                      //inclination
                     NULL,                      //periapsis
                     LOCKED_SURFACE_SPIN_EVOL_MODE);
    detect_stellar_wind_saturation(primary);

    OrbitSolver *solver = evolve_system(
        system,
        10.0,   //final age
        1e-3,   //max timestep
        1e-6,   //precision
        NULL,   //required ages
        0       //num required ages
    );
    int num_steps = num_evolution_steps(solver);
    double *age = new double[num_steps],
           *semimajor = new double[num_steps],
           *primary_lconv = new double[num_steps],
           *primary_lrad = new double[num_steps],
           *secondary_lconv = new double[num_steps],
           *secondary_lrad = new double[num_steps];
    get_star_star_evolution(
        solver,
        system,
        primary,
        secondary,
        age,
        semimajor,
        NULL,           //eeccentricity
        NULL,           //primary envelope inclination
        NULL,           //primary core inclination
        NULL,           //primary envelope periapsis
        NULL,           //primary core periapsis
        primary_lconv,  //primary envelope angmom
        primary_lrad,   //primary core angmom
        NULL,           //secondary envelope inclination
        NULL,           //secondary core inclination
        NULL,           //secondary envelope periapsis
        NULL,           //secondary core periapsis
        secondary_lconv,//secondary envelope angmom
        secondary_lrad, //secondary core angmom
        NULL,           //evolution mode
        NULL,           //primary wind saturation
        NULL            //secondary wind saturation
    );
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setw(25) << "Age[Gyr]"
              << std::setw(25) << "worb[rad/day]"
              << std::setw(25) << "prim_wconv[rad/day]"
              << std::setw(25) << "prim_wrad[rad/day]"
              << std::setw(25) << "sec_wconv[rad/day]"
              << std::setw(25) << "sec_wrad[rad/day]"
              << std::setw(25) << "prim_Lconv"
              << std::setw(25) << "prim_Lrad"
              << std::setw(25) << "sec_Lconv"
              << std::setw(25) << "sec_Lrad"
              << std::endl;
    for(int i = 0; i < num_steps; ++i)
        std::cout << std::setw(25) << age[i]
                  << std::setw(25) << Core::orbital_angular_velocity(
                      mprimary,
                      msecondary,
                      semimajor[i]
                  )
                  << std::setw(25) << (primary_lconv[i]
                                       /
                                       envelope_inertia(primary, age[i]))
                  << std::setw(25) << (primary_lrad[i]
                                       /
                                       core_inertia(primary, age[i]))
                  << std::setw(25) << (secondary_lconv[i]
                                       /
                                       envelope_inertia(secondary, age[i]))
                  << std::setw(25) << (secondary_lrad[i]
                                       /
                                       core_inertia(secondary, age[i]))
                  << std::setw(25) << primary_lconv[i]
                  << std::setw(25) << primary_lrad[i]
                  << std::setw(25) << secondary_lconv[i]
                  << std::setw(25) << secondary_lrad[i]
                  << std::endl;

    destroy_binary(system);
    destroy_star(primary);
    destroy_star(secondary);
    destroy_interpolator(interpolator);
    destroy_solver(solver);
    delete[] age;
    delete[] semimajor;
    delete[] primary_lconv;
    delete[] primary_lrad;
    delete[] secondary_lconv;
    delete[] secondary_lrad;
}
