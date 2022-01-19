#include "/home/kpenev/projects/git/poet/poet_src/Planet/CInterface.h"
#include "/home/kpenev/projects/git/poet/poet_src/StellarEvolution/CInterface.h"
#include "/home/kpenev/projects/git/poet/poet_src/Star/CInterface.h"
#include "/home/kpenev/projects/git/poet/poet_src/Evolve/CInterface.h"
#include "/home/kpenev/projects/git/poet/poet_src/Core/OrbitalExpressions.h"
#include "/home/kpenev/projects/git/poet/poet_src/Core/AstronomicalConstants.h"
#include "dirent_hacked.h"
#include <iomanip>
#include <iostream>
#include <valarray>

/******************************************************************************
 * PRIMARY PARAMETERS                                                         *
 ******************************************************************************/
/* Physical properties */
const double PRIMARY_MASS = 1.1391152261289241e+00;
const double FEH = -3.5825605105145464e-01;

/* Wind */
const double PRIMARY_WIND_STRENGTH = 1.7000000000000001e-01;
const double PRIMARY_WIND_SATURATION_FREQUENCY =
    2.5400000000000000e+00;
const double PRIMARY_DIFF_ROT_COUPLING_TIMESCALE =
    5.0000000000000001e-03;

/* Dissipation */
#if 1
    const double PRIMARY_PHASE_LAG = 2.9841551829730380e-06;
    double PRIMARY_TIDAL_FREQUENCY_BREAKS[] =
        {6.4463263225813347e-02};
    double PRIMARY_SPIN_FREQUENCY_BREAKS[] =
        {};
    double PRIMARY_TIDAL_FREQUENCY_POWERS[] =
        {3.0000000000000000e+00, 3.0000000000000000e+00};
    double PRIMARY_SPIN_FREQUENCY_POWERS[] =
        {0.0000000000000000e+00};
#endif
/******************************************************************************/

/******************************************************************************
 * SECONDARY PARAMETERS                                                       *
 ******************************************************************************/
const double SECONDARY_MASS = 5.0440022212988767e-01;

#if 1
    /* Wind */
    const double SECONDARY_WIND_STRENGTH = 1.7000000000000001e-01;
    const double SECONDARY_WIND_SATURATION_FREQUENCY =
        2.5400000000000000e+00;
    const double SECONDARY_DIFF_ROT_COUPLING_TIMESCALE =
        5.0000000000000001e-03;
#else
    const double SECONDARY_RADIUS = nan;
#endif

/* Dissipation */
#if 1
    const double SECONDARY_PHASE_LAG = 2.9841551829730380e-06;
    double SECONDARY_TIDAL_FREQUENCY_BREAKS[] =
        {6.4463263225813347e-02};
    double SECONDARY_SPIN_FREQUENCY_BREAKS[] =
        {};
    double SECONDARY_TIDAL_FREQUENCY_POWERS[] =
        {3.0000000000000000e+00, 3.0000000000000000e+00};
    double SECONDARY_SPIN_FREQUENCY_POWERS[] =
        {0.0000000000000000e+00};
#endif
/******************************************************************************/


/******************************************************************************
 * INITIAL ORBIT                                                              *
 ******************************************************************************/
const double INITIAL_SEMIMAJOR = 1.2438981391896482e+01;
const double DISK_FREQUENCY = 4.0999999999999996e+00;
const double DISK_DISSIPATION_AGE = 5.0000000000000001e-03;
const double INITIAL_INCLINATION = 0.0000000000000000e+00;
const double INITIAL_ECCENTRICITY = 7.3999999999999996e-02;
/******************************************************************************/

const double ZERO=0.0;

EvolvingStar *create_primary()
{
    MESAInterpolator *primary_interpolator = load_interpolator(
        "/home/ruskin/projects/poet/stellar_evolution_interpolators/f18f8304-fa82-4fe1-b082-43bbf8131cae"
    );

    EvolvingStar *primary = create_star(PRIMARY_MASS,
                                        FEH,
                                        PRIMARY_WIND_STRENGTH,
                                        PRIMARY_WIND_SATURATION_FREQUENCY,
                                        PRIMARY_DIFF_ROT_COUPLING_TIMESCALE,
                                        primary_interpolator);

    select_interpolation_region(primary, core_formation_age(primary));

#if 1
    set_star_dissipation(
        primary,
        0,          //zone index
        sizeof(PRIMARY_TIDAL_FREQUENCY_BREAKS) / sizeof(double),
        sizeof(PRIMARY_SPIN_FREQUENCY_BREAKS) / sizeof(double),
        (
            sizeof(PRIMARY_TIDAL_FREQUENCY_BREAKS)
            ? PRIMARY_TIDAL_FREQUENCY_BREAKS
            : NULL
        ),
        (
            sizeof(PRIMARY_SPIN_FREQUENCY_BREAKS)
            ? PRIMARY_SPIN_FREQUENCY_BREAKS
            : NULL
        ),
        PRIMARY_TIDAL_FREQUENCY_POWERS,
        PRIMARY_SPIN_FREQUENCY_POWERS,
        PRIMARY_PHASE_LAG
    );
#endif

    return primary;
}

#if 1
EvolvingStar *create_secondary()
{
    MESAInterpolator *secondary_interpolator = load_interpolator(
        "/home/ruskin/projects/poet/stellar_evolution_interpolators/f18f8304-fa82-4fe1-b082-43bbf8131cae"
    );

    double initial_secondary_angmom[] = {
        4.7025472297622317e-01,
        0.0000000000000000e+00
    };

    EvolvingStar *secondary = create_star(SECONDARY_MASS,
                                          FEH,
                                          0.0,
                                          1e10,
                                          SECONDARY_DIFF_ROT_COUPLING_TIMESCALE,
                                          secondary_interpolator);

    select_interpolation_region(secondary, DISK_DISSIPATION_AGE);
#if 1
    set_star_dissipation(
        secondary,
        0,          //zone index
        sizeof(SECONDARY_TIDAL_FREQUENCY_BREAKS) / sizeof(double),
        sizeof(SECONDARY_SPIN_FREQUENCY_BREAKS) / sizeof(double),
        (
            sizeof(SECONDARY_TIDAL_FREQUENCY_BREAKS)
            ? SECONDARY_TIDAL_FREQUENCY_BREAKS
            : NULL
        ),
        (
            sizeof(SECONDARY_SPIN_FREQUENCY_BREAKS)
            ? SECONDARY_SPIN_FREQUENCY_BREAKS
            : NULL
        ),
        SECONDARY_TIDAL_FREQUENCY_POWERS,
        SECONDARY_SPIN_FREQUENCY_POWERS,
        SECONDARY_PHASE_LAG
    );
#endif

    configure_star(secondary,
                   DISK_DISSIPATION_AGE,        //formation age
                   PRIMARY_MASS,                //companion mass
                   INITIAL_SEMIMAJOR,           //formation semimajor
                   INITIAL_ECCENTRICITY,        //formation eccentricity
                   initial_secondary_angmom,    //spin angular momentum
                   &ZERO,                       //inclination
                   &ZERO,                       //periapsis
                   false,                       //locked surface?
                   true,                        //zero outer inclination?
                   true);                       //zero outer periapsis?
    detect_stellar_wind_saturation(secondary);
    return secondary;
}
#else
CPlanet *create_secondary()
{
    CPlanet *secondary = create_planet(SECONDARY_MASS, SECONDARY_RADIUS);

#if 1
    set_planet_dissipation(
        secondary,
        sizeof(SECONDARY_TIDAL_FREQUENCY_BREAKS) / sizeof(double),
        sizeof(SECONDARY_SPIN_FREQUENCY_BREAKS) / sizeof(double),
        (
            sizeof(SECONDARY_TIDAL_FREQUENCY_BREAKS)
            ? SECONDARY_TIDAL_FREQUENCY_BREAKS
            : NULL
        ),
        (
            sizeof(SECONDARY_SPIN_FREQUENCY_BREAKS)
            ? SECONDARY_SPIN_FREQUENCY_BREAKS
            : NULL
        ),
        SECONDARY_TIDAL_FREQUENCY_POWERS,
        SECONDARY_SPIN_FREQUENCY_POWERS,
        SECONDARY_PHASE_LAG
    );
#endif

    configure_planet(secondary,
                     DISK_DISSIPATION_AGE,
                     PRIMARY_MASS,
                     INITIAL_SEMIMAJOR,
                     INITIAL_ECCENTRICITY,
                     &ZERO, //Spin angular momentum
                     NULL,  //Initial eccentricities
                     NULL,  //Initial periapses
                     false, //Locked to disk?
                     true,  //Outer inclination == 0?
                     true); //Outer periapsis == 0?
    return secondary;
}
#endif

DiskBinarySystem *create_system(
    EvolvingStar *primary,
#if 1
    EvolvingStar *
#else
    CPlanet *
#endif
    secondary
)
{

    return
#if 1
    create_star_star_system(
#else
    create_star_planet_system(
#endif
        primary,
        secondary,
        INITIAL_SEMIMAJOR,
        INITIAL_ECCENTRICITY,
        INITIAL_INCLINATION,
        DISK_FREQUENCY,
        DISK_DISSIPATION_AGE,
        DISK_DISSIPATION_AGE
    );
}

int main(int, char **)
{
    read_eccentricity_expansion_coefficients(
        "/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/eccentricity_expansion_coef.txt"
    );
    EvolvingStar *primary=create_primary();

#if 1
    EvolvingStar *
#else
    CPlanet *
#endif
    secondary = create_secondary();

    DiskBinarySystem *system=create_system(primary, secondary);

    configure_system(system,
                     core_formation_age(primary),
                     NaN,                       //semimajor
                     NaN,                       //eccentricity
                     &ZERO,                     //spin angmom
                     NULL,                      //inclination
                     NULL,                      //periapsis
                     LOCKED_SURFACE_SPIN_EVOL_MODE);
    detect_stellar_wind_saturation(primary);

    OrbitSolver *solver = evolve_system(
        system,
        4.0312213550728943e-01,
        1.0000000000000000e-03,
        9.9999999999999995e-07,
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
#if 1
    get_star_star_evolution(
#else
    get_star_planet_evolution(
#endif
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
#if 1
        NULL,           //secondary envelope inclination
        NULL,           //secondary core inclination
        NULL,           //secondary envelope periapsis
        NULL,           //secondary core periapsis
        secondary_lconv,//secondary envelope angmom
        secondary_lrad, //secondary core angmom
#endif
        NULL,           //evolution mode
        NULL           //primary wind saturation
#if 1
        , NULL            //secondary wind saturation
#endif
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



#if 1
        secondary_Iconv = (
            age[i] <= 2e-3
            ? Core::NaN
            : envelope_inertia(secondary, age[i])
        );
        if(
            age[i] < core_formation_age(secondary)
            ||
            age[i] > lifetime(secondary)
        )
            secondary_Irad = Core::NaN;
        else {
            secondary_Irad = core_inertia(secondary, age[i]);
        }

#else
        secondary_Iconv = secondary_Irad = Core::NaN;
#endif

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
//    destroy_interpolator(primary_interpolator);
#if 1
    destroy_star(secondary);
//    destroy_interpolator(secondary_interpolator);
#endif
    destroy_solver(solver);
    delete[] age;
    delete[] semimajor;
    delete[] primary_lconv;
    delete[] primary_lrad;
    delete[] secondary_lconv;
    delete[] secondary_lrad;
}
