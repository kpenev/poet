#include "%(poet_include_path)s/Planet/CInterface.h"
#include "%(poet_include_path)s/StellarEvolution/CInterface.h"
#include "%(poet_include_path)s/Star/CInterface.h"
#include "%(poet_include_path)s/Evolve/CInterface.h"
#include "%(poet_include_path)s/Core/OrbitalExpressions.h"
#include "%(poet_include_path)s/Core/AstronomicalConstants.h"
#include "dirent_hacked.h"
#include <iomanip>
#include <iostream>
#include <valarray>

/******************************************************************************
 * PRIMARY PARAMETERS                                                         *
 ******************************************************************************/
/* Physical properties */
const double PRIMARY_MASS = %(primary_mass).16e;
const double FEH = %(feh).16e;

/* Wind */
const double PRIMARY_WIND_STRENGTH = %(primary_wind_strength).16e;
const double PRIMARY_WIND_SATURATION_FREQUENCY =
    %(primary_wind_saturation_frequency).16e;
const double PRIMARY_DIFF_ROT_COUPLING_TIMESCALE =
    %(primary_diff_rot_coupling_timescale).16e;

/* Dissipation */
#if %(dissipative_primary)d
    const double PRIMARY_PHASE_LAG = %(primary_reference_phase_lag).16e;
    double PRIMARY_TIDAL_FREQUENCY_BREAKS[] =
        %(primary_tidal_frequency_breaks)s;
    double PRIMARY_SPIN_FREQUENCY_BREAKS[] =
        %(primary_spin_frequency_breaks)s;
    double PRIMARY_TIDAL_FREQUENCY_POWERS[] =
        %(primary_tidal_frequency_powers)s;
    double PRIMARY_SPIN_FREQUENCY_POWERS[] =
        %(primary_spin_frequency_powers)s;
#endif
/******************************************************************************/

/******************************************************************************
 * SECONDARY PARAMETERS                                                       *
 ******************************************************************************/
const double SECONDARY_MASS = %(secondary_mass).16e;

#if %(secondary_is_star)d
    /* Wind */
    const double SECONDARY_WIND_STRENGTH = %(secondary_wind_strength).16e;
    const double SECONDARY_WIND_SATURATION_FREQUENCY =
        %(secondary_wind_saturation_frequency).16e;
    const double SECONDARY_DIFF_ROT_COUPLING_TIMESCALE =
        %(secondary_diff_rot_coupling_timescale).16e;
#else
    const double SECONDARY_RADIUS = %(secondary_radius).16e;
#endif

/* Dissipation */
#if %(dissipative_secondary)d
    const double SECONDARY_PHASE_LAG = %(secondary_reference_phase_lag).16e;
    double SECONDARY_TIDAL_FREQUENCY_BREAKS[] =
        %(secondary_tidal_frequency_breaks)s;
    double SECONDARY_SPIN_FREQUENCY_BREAKS[] =
        %(secondary_spin_frequency_breaks)s;
    double SECONDARY_TIDAL_FREQUENCY_POWERS[] =
        %(secondary_tidal_frequency_powers)s;
    double SECONDARY_SPIN_FREQUENCY_POWERS[] =
        %(secondary_spin_frequency_powers)s;
#endif
/******************************************************************************/


/******************************************************************************
 * INITIAL ORBIT                                                              *
 ******************************************************************************/
const double INITIAL_SEMIMAJOR = %(initial_semimajor).16e;
const double DISK_FREQUENCY = %(disk_lock_frequency).16e;
const double DISK_DISSIPATION_AGE = %(disk_dissipation_age).16e;
const double INITIAL_INCLINATION = %(initial_inclination).16e;
const double INITIAL_ECCENTRICITY = %(initial_eccentricity).16e;
/******************************************************************************/

const double ZERO=0.0;

EvolvingStar *create_primary()
{
    MESAInterpolator *primary_interpolator = load_interpolator(
        "%(primary_interpolator_fname)s"
    );

    EvolvingStar *primary = create_star(PRIMARY_MASS,
                                        FEH,
                                        PRIMARY_WIND_STRENGTH,
                                        PRIMARY_WIND_SATURATION_FREQUENCY,
                                        PRIMARY_DIFF_ROT_COUPLING_TIMESCALE,
                                        primary_interpolator);

    select_interpolation_region(primary, core_formation_age(primary));

#if %(dissipative_primary)d
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

#if %(secondary_is_star)d
EvolvingStar *create_secondary()
{
    MESAInterpolator *secondary_interpolator = load_interpolator(
        "%(secondary_interpolator_fname)s"
    );

    double initial_secondary_angmom[] = {
        %(initial_secondary_envelope_angmom).16e
        %(initial_secondary_core_angmom).16e,
    };

    EvolvingStar *secondary = create_star(SECONDARY_MASS,
                                          FEH,
                                          0.0,
                                          1e10,
                                          DIFF_ROT_COUPLING_TIMESCALE,
                                          secondary_interpolator);

    select_interpolation_region(secondary, DISK_DISSIPATION_AGE);
#if %(dissipative_secondary)d
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

#if %(dissipative_secondary)d
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
#if %(secondary_is_star)d
    EvolvingStar *
#else
    CPlanet *
#endif
    secondary
)
{

    return
#if %(secondary_is_star)d
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
        "%(eccentricity_expansion_fname)s"
    );
    EvolvingStar *primary=create_primary();

#if %(secondary_is_star)d
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
        %(final_age).16e,
        %(max_time_step).16e,
        %(precision).16e,
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
#if %(secondary_is_star)d
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
#if %(secondary_is_star)d
        NULL,           //secondary envelope inclination
        NULL,           //secondary core inclination
        NULL,           //secondary envelope periapsis
        NULL,           //secondary core periapsis
        secondary_lconv,//secondary envelope angmom
        secondary_lrad, //secondary core angmom
#endif
        NULL,           //secondary envelope inclination
        NULL,           //secondary core inclination
        NULL,           //secondary core periapsis
        NULL,           //evolution mode
        NULL           //primary wind saturation
#if %(secondary_is_star)d
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



#if %(secondary_is_star)d
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
#if %(secondary_is_star)d
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
