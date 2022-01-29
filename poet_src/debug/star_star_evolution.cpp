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
    MESAInterpolator *interpolator=NULL;
    DIR *dirstream = opendir(interpolator_dir.c_str());
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

            if(interpolator)
                throw Core::Error::IO(
                    "Multiple candidate interpolators fund in "
                    +
                    interpolator_dir
                );
            interpolator = load_interpolator(
                (interpolator_dir + fname).c_str()
            );

        }
    }
    if(!interpolator)
        throw Core::Error::IO(
            "No interpolators found in "
            +
            interpolator_dir
        );
    return interpolator;
}

int main(int, char **)
{

    const double PRIMARY_MASS = 0.9920654905378268;
    const double SECONDARY_MASS = 0.7230138299345078;
    const double FEH = 0.23115300565169195;
    const double INITIAL_PERIOD = 16.145221014708316;
    const double LGQ_MIN = 5.186227042681173;
    const double LGQ_BREAK_PERIOD = 8.978669081328178;
    const double LGQ_POWERLAW = -3.1395911900486437;

    const double DISK_PERIOD = 10.0;
    const double DISK_DISSIPATION_AGE = 0.02;
    const double WIND_SATURATION_FREQUENCY = 2.54;
    const double WIND_STRENGTH = 0.17;
    const double DIFF_ROT_COUPLING_TIMESCALE = 5e-3;
    const double INITIAL_ECCENTRICITY = 0.8;
    const double INCLINATION = 0.0;

    const double INITIAL_SEMIMAJOR = Core::semimajor_from_period(
        PRIMARY_MASS,
        SECONDARY_MASS,
        INITIAL_PERIOD
    );
    const double DISK_FREQUENCY = 2.0 * M_PI / DISK_PERIOD;
    const double REF_PHASE_LAG = 15.0 / (16.0 * M_PI * std::pow(10.0, LGQ_MIN));

    std::cerr << "Starting evolution with a0 = " << INITIAL_SEMIMAJOR << std::endl;

    prepare_eccentricity_expansion(
        "eccentricity_expansion_coef_O400.sqlite",
        1e-4,
        true,
        false
    );
    MESAInterpolator *primary_interpolator = get_interpolator(
        "stellar_evolution_interpolators_bkp/"
    );


    double zero = 0.0;
    double initial_secondary_angmom[] = {0.2, 0.001};

    double tidal_frequency_breaks[] = {2.0 * M_PI / LGQ_BREAK_PERIOD};
    double tidal_frequency_powers[] = {
        std::max(0.0, LGQ_POWERLAW),
        std::min(0.0, LGQ_POWERLAW)
    };


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
                         tidal_frequency_breaks,       //tidal frequency breaks
                         NULL,       //spin frequency breaks
                         tidal_frequency_powers,      //tidal frequency powers
                         &zero,      //spin frequency powers
                         REF_PHASE_LAG,
                         1.0,
                         0.0);

    EvolvingStar *secondary = create_star(SECONDARY_MASS,
                                          FEH,
                                          WIND_STRENGTH,
                                          WIND_SATURATION_FREQUENCY,
                                          DIFF_ROT_COUPLING_TIMESCALE,
                                          primary_interpolator);
    select_interpolation_region(secondary, DISK_DISSIPATION_AGE);
    set_star_dissipation(secondary,
                         0,          //zone index
                         0,          //# tidal frequency breaks
                         0,          //# spin frequency breaks
                         tidal_frequency_breaks,       //tidal frequency breaks
                         NULL,       //spin frequency breaks
                         tidal_frequency_powers,      //tidal frequency powers
                         &zero,      //spin frequency powers
                         REF_PHASE_LAG,
                         1.0,
                         0.0);

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
        DISK_FREQUENCY,             //disk lock frequency
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
        7.0,    //final age
        1e-3,   //max timestep
        1e-6,   //precision
        NULL,   //required ages
        0,      //num required ages
        true,   //Print stepping progress?
        0
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
        NULL,           //secondary wind saturation
        NULL,           //semimajor rate
        NULL,           //eeccentricity rate
        NULL,           //primary envelope inclination rate
        NULL,           //primary core inclination rate
        NULL,           //primary envelope periapsis rate
        NULL,           //primary core periapsis rate
        NULL,           //primary envelope angmom rate
        NULL,           //primary core angmom rate
        NULL,           //secondary envelope inclination rate
        NULL,           //secondary core inclination rate
        NULL,           //secondary envelope periapsis rate
        NULL,           //secondary core periapsis rate
        NULL,           //secondary envelope angmom rate
        NULL            //secondary core angmom rate
    );
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setw(25) << "Age[Gyr]"            //1
              << std::setw(25) << "worb[rad/day]"       //2
              << std::setw(25) << "eccentricity"        //3
              << std::setw(25) << "prim_wconv[rad/day]" //4
              << std::setw(25) << "prim_wrad[rad/day]"  //5
              << std::setw(25) << "sec_wconv[rad/day]"  //6
              << std::setw(25) << "sec_wrad[rad/day]"   //7
              << std::setw(25) << "prim_Lconv"          //8
              << std::setw(25) << "prim_Lrad"           //9
              << std::setw(25) << "sec_Lconv"           //10
              << std::setw(25) << "sec_Lrad"            //11
              << std::setw(25) << "prim_Rconv"          //12
              << std::setw(25) << "prim_Rrad"           //13
              << std::setw(25) << "sec_Rconv"           //14
              << std::setw(25) << "sec_Rrad"            //15
              << std::endl;
    double primary_Iconv,
           primary_Irad,
           secondary_Iconv,
           secondary_Irad,
           primary_rconv,
           primary_rrad,
           secondary_rconv,
           secondary_rrad;
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
            primary_rconv = star_radius(primary, age[i]);
            secondary_rconv = star_radius(secondary, age[i]);
            primary_rrad = core_radius(primary, age[i]);
            secondary_rrad = core_radius(secondary, age[i]);
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
                  << std::setw(25) << primary_rconv
                  << std::setw(25) << primary_rrad
                  << std::setw(25) << secondary_rconv
                  << std::setw(25) << secondary_rrad
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
