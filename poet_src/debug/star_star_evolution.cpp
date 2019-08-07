#include "../Planet/CInterface.h"
#include "../StellarEvolution/CInterface.h"
#include "../Star/CInterface.h"
#include "../Evolve/CInterface.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/AstronomicalConstants.h"
#include "dirent_hacked.h"
#include <iomanip>
#include <iostream>
#include <valarray>

MESAInterpolator *get_interpolator(const std::string &interpolator_dir)
{
    MESAInterpolator *interpolator;
    DIR *dirstream = opendir(interpolator_dir.c_str());
    bool found = false;
    for(struct dirent *entry; (entry = readdir(dirstream));) {
        std::string fname(entry->d_name);
        std::cout << "Fname: " << fname << std::endl;
        if(
            fname[0] != '.'
            &&
            fname.substr(fname.size() - 7) != ".sqlite"
        ) {
            std::cout << "Fname tail: "
                      << std::string(fname.substr(fname.size() - 7))
                      << std::endl;

            assert(!found);
            found=true;
            interpolator = load_interpolator(
                (interpolator_dir + fname).c_str()
            );

        }
    }
    return interpolator;
}

int main(int, char **)
{

    const double PRIMARY_MASS = 1.02505108176533;
    const double SECONDARY_MASS = 0.9994248047211968;
    const double FEH = -0.06;
    const double INITIAL_PERIOD = 14.2;
    const double INITIAL_SEMIMAJOR = Core::semimajor_from_period(
        PRIMARY_MASS,
        SECONDARY_MASS,
        INITIAL_PERIOD
    );

    std::cerr << "Starting evolution with a0 = " << INITIAL_SEMIMAJOR << std::endl;

    const double DISK_PERIOD = 4.687960373696023;
    const double PRIMARY_PHASE_LAG =2.984155182973038e-09;
    const double SECONDARY_PHASE_LAG =2.984155182973038e-09;
    const double DISK_DISSIPATION_AGE = 5e-3;
    const double WIND_SATURATION_FREQUENCY = 2.54;
    const double DIFF_ROT_COUPLING_TIMESCALE = 5e-3;
    const double WIND_STRENGTH = 0.17;
    const double INCLINATION = 0.0;
    const double INITIAL_ECCENTRICITY = 0.0241;

    read_eccentricity_expansion_coefficients(
        "eccentricity_expansion_coef_O200.txt"
    );
    MESAInterpolator *primary_interpolator = get_interpolator(
        "stellar_evolution_interpolators/"
    );


    double zero = 0.0;
    double initial_secondary_angmom[] = {0.0, 0.0};

    EvolvingStar *primary = create_star(PRIMARY_MASS,
                                        FEH,
                                        WIND_STRENGTH,
                                        WIND_SATURATION_FREQUENCY,
                                        DIFF_ROT_COUPLING_TIMESCALE,
                                        primary_interpolator);
    select_interpolation_region(primary, core_formation_age(primary));
    set_star_dissipation(primary,
                         0,          //zone index
                         0,          //# tidal frequency breaks
                         0,          //# spin frequency breaks
                         NULL,       //tidal frequency breaks
                         NULL,       //spin frequency breaks
                         &zero,      //tidal frequency powers
                         &zero,      //spin frequency powers
                         PRIMARY_PHASE_LAG);

    EvolvingStar *secondary = create_star(SECONDARY_MASS,
                                          FEH,
                                          0.0,
                                          1e10,
                                          DIFF_ROT_COUPLING_TIMESCALE,
                                          primary_interpolator);
    select_interpolation_region(secondary, DISK_DISSIPATION_AGE);
    set_star_dissipation(secondary,
                         0,          //zone index
                         0,          //# tidal frequency breaks
                         0,          //# spin frequency breaks
                         NULL,       //tidal frequency breaks
                         NULL,       //spin frequency breaks
                         &zero,      //tidal frequency powers
                         &zero,      //spin frequency powers
                         SECONDARY_PHASE_LAG);

    configure_star(secondary,
                   DISK_DISSIPATION_AGE,        //formation age
                   PRIMARY_MASS,                //companion mass
                   INITIAL_SEMIMAJOR,           //formation semimajor
                   INITIAL_ECCENTRICITY,        //formation eccentricity
                   initial_secondary_angmom,    //spin angular momentum
                   &zero,                       //inclination
                   &zero,                       //periapsis
                   false,                       //locked surface?
                   true,                        //zero outer inclination?
                   true);                       //zero outer periapsis?
    detect_stellar_wind_saturation(secondary);

    DiskBinarySystem *system = create_star_star_system(
        primary,                    //primary
        secondary,                  //secondary
        INITIAL_SEMIMAJOR,          //initial semimajor
        INITIAL_ECCENTRICITY,       //initial eccentricity
        INCLINATION,                //initial inclination
        2.0 * M_PI / DISK_PERIOD,   //disk lock frequency
        DISK_DISSIPATION_AGE,       //disk dissipation age
        DISK_DISSIPATION_AGE        //secondary formation age
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
        8.705870938160976,    //final age
        1e-3,   //max timestep
        1e-6,   //precision
        NULL,   //required ages
        0,      //num required ages
        true    //Print stepping progress?
    );
    int num_steps = num_evolution_steps(solver);
    double *age = new double[num_steps],
           *semimajor = new double[num_steps],
           *eccentricity = new double[num_steps],
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
        eccentricity,   //eeccentricity
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
              << std::setw(25) << "eccentricity"
              << std::setw(25) << "prim_wconv[rad/day]"
              << std::setw(25) << "prim_wrad[rad/day]"
              << std::setw(25) << "sec_wconv[rad/day]"
              << std::setw(25) << "sec_wrad[rad/day]"
              << std::setw(25) << "prim_Lconv"
              << std::setw(25) << "prim_Lrad"
              << std::setw(25) << "sec_Lconv"
              << std::setw(25) << "sec_Lrad"
              << std::endl;
    double primary_Iconv, primary_Irad, secondary_Iconv, secondary_Irad;
    for(int i = 0; i < num_steps; ++i) {
        if(
            age[i] < core_formation_age(primary)
            ||
            age[i] > lifetime(primary)
        )
            primary_Iconv = primary_Irad = Core::NaN;
        else {
            primary_Iconv = envelope_inertia(primary, age[i]);
            primary_Irad = core_inertia(primary, age[i]);
        }


        secondary_Iconv = (age[i] <= 2e-3
                           ? Core::NaN
                           : envelope_inertia(secondary, age[i]));
        if(
            age[i] < core_formation_age(secondary)
            ||
            age[i] > lifetime(secondary)
        )
            secondary_Irad = Core::NaN;
        else {
            secondary_Irad = core_inertia(secondary, age[i]);
        }

        std::cout << std::setw(25) << age[i]
                  << std::setw(25) << Core::orbital_angular_velocity(
                      PRIMARY_MASS,
                      SECONDARY_MASS,
                      semimajor[i]
                  )
                  << std::setw(25) << eccentricity[i]
                  << std::setw(25) << primary_lconv[i] / primary_Iconv
                  << std::setw(25) << primary_lrad[i] / primary_Irad
                  << std::setw(25) << secondary_lconv[i] / secondary_Iconv
                  << std::setw(25) << secondary_lrad[i] / secondary_Irad
                  << std::setw(25) << primary_lconv[i]
                  << std::setw(25) << primary_lrad[i]
                  << std::setw(25) << secondary_lconv[i]
                  << std::setw(25) << secondary_lrad[i]
                  << std::endl;
    }

    destroy_binary(system);
    destroy_star(primary);
    destroy_star(secondary);
    destroy_interpolator(primary_interpolator);
    destroy_solver(solver);
    delete[] age;
    delete[] semimajor;
    delete[] primary_lconv;
    delete[] primary_lrad;
    delete[] secondary_lconv;
    delete[] secondary_lrad;
}
