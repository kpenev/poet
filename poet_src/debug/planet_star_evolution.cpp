#include "../Planet/CInterface.h"
#include "../StellarEvolution/CInterface.h"
#include "../Star/CInterface.h"
#include "../Evolve/CInterface.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/AstronomicalConstants.h"
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
    const double Mjup_to_Msun = (Core::AstroConst::jupiter_mass
                                 /
                                 Core::AstroConst::solar_mass),
                 Rjup_to_Msun = (Core::AstroConst::jupiter_radius
                                 /
                                 Core::AstroConst::solar_radius);

    const double STAR_MASS = 0.8021362093562747;
    const double PLANET_MASS = 0.006909308339572223;
    const double PLANET_RADIUS = 0.09464394270895214;
    const double FEH = 0.22600928;
//    const double INITIAL_PERIOD = ;
    const double LGQ_MIN = 2.0;
    const double LGQ_BREAK_PERIOD = 1.0;
    const double LGQ_POWERLAW = 0;
    const bool DISSIPATIVE_STAR = false;
    const bool DISSIPATIVE_PLANET = true;
    const double FINAL_AGE = 4.495718172918945;

//    double initial_secondary_angmom[] = {1.43525535, 0.43099626};
    double initial_secondary_angmom[] = {1.79981886e-05};


    const double DISK_PERIOD = 6.48175405473875;
    const double DISK_DISSIPATION_AGE = 0.02;
    const double WIND_SATURATION_FREQUENCY = 2.78;
    const double WIND_STRENGTH = 0.17;
    const double DIFF_ROT_COUPLING_TIMESCALE = 5e-2;
    const double INITIAL_ECCENTRICITY = 0.78;
    const double INCLINATION = 0.0;
    const double LOCK_PERIOD = 20.0;

    const double INITIAL_SEMIMAJOR = 7.929015451160677;
/*    Core::semimajor_from_period(
        PRIMARY_MASS,
        SECONDARY_MASS,
        INITIAL_PERIOD
    );*/
    const double DISK_FREQUENCY = 2.0 * M_PI / DISK_PERIOD;
    double ref_phase_lag = 15.0 / (16.0 * M_PI * std::pow(10.0, LGQ_MIN));


    MESAInterpolator *interpolator = get_interpolator(
        "../../stellar_evolution_interpolators/"
    );

    prepare_eccentricity_expansion(
        "../../eccentricity_expansion_coef_O400.sqlite",
        1e-4,
        true,
        true
    );

    CPlanet *planet = create_planet(PLANET_MASS, PLANET_RADIUS);
    configure_planet(planet,
                     DISK_DISSIPATION_AGE,
                     STAR_MASS,
                     INITIAL_SEMIMAJOR,
                     INITIAL_ECCENTRICITY,
                     initial_secondary_angmom,
                     NULL,
                     NULL,
                     false,
                     true,
                     true);

    EvolvingStar *star = create_star(STAR_MASS,
                                     FEH,
                                     WIND_STRENGTH,
                                     WIND_SATURATION_FREQUENCY,
                                     DIFF_ROT_COUPLING_TIMESCALE,
                                     interpolator);
    select_interpolation_region(star, core_formation_age(star));

    double zero = 0.0;
    unsigned num_breaks;
    if(LGQ_POWERLAW == 0)
        num_breaks = 0;
    else if(LGQ_POWERLAW > 0)
        num_breaks = 2;
    else
        num_breaks = 1;

    std::valarray<double> tidal_frequency_breaks(num_breaks),
                          tidal_frequency_powers(num_breaks + 1);

    tidal_frequency_powers[0] = 0.0;
    if(LGQ_POWERLAW > 0) {
        //ref_phase_lag *= std::pow(LGQ_BREAK_PERIOD / LOCK_PERIOD, LGQ_POWERLAW);
        tidal_frequency_breaks[0] = 2.0 * M_PI / LOCK_PERIOD;
        tidal_frequency_breaks[1] = 0.67344225;//2.0 * M_PI / LGQ_BREAK_PERIOD;
        tidal_frequency_powers[2] = 0.0;
    } else if(LGQ_POWERLAW < 0) {
        tidal_frequency_breaks[0] = 2.0 * M_PI / LGQ_BREAK_PERIOD;
    }
    if(LGQ_POWERLAW != 0)
        tidal_frequency_powers[1] = LGQ_POWERLAW;


    if(DISSIPATIVE_STAR)
        set_star_dissipation(
            star,
            0,                            //zone index
            tidal_frequency_breaks.size(),//# tidal frequency breaks
            0,                            //# spin frequency breaks
            0,                            //# age breaks
            &(tidal_frequency_breaks[0]), //tidal frequency breaks
            NULL,                         //spin frequency breaks
            &(tidal_frequency_powers[0]), //tidal frequency powers
            &zero,                        //spin frequency powers
            NULL,                         //age breaks
            &ref_phase_lag,
            1.0,
            0.0
        );

    if(DISSIPATIVE_PLANET)
        set_planet_dissipation(
            planet,
            tidal_frequency_breaks.size(),//# tidal frequency breaks
            0,                            //# spin frequency breaks
            0,                            //# age breaks
            &(tidal_frequency_breaks[0]), //tidal frequency breaks
            NULL,                         //spin frequency breaks
            &(tidal_frequency_powers[0]), //tidal frequency powers
            &zero,                        //spin frequency powers
            NULL,                         //age breaks
            &ref_phase_lag,
            1.0,
            0.0
        );

    DiskBinarySystem *system = create_star_planet_system(
        star,
        planet,
        INITIAL_SEMIMAJOR,
        INITIAL_ECCENTRICITY,
        0.0, //initial inclination
        DISK_FREQUENCY,
        DISK_DISSIPATION_AGE,
        DISK_DISSIPATION_AGE
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

    OrbitSolver *solver;
    try {
        solver = evolve_system(
            system,
            FINAL_AGE,
            1e-3,
            1e-6,
            NULL,
            0,
            true,
            0,
            0
        );
    } catch(std::exception &exception) {
        std::cerr << "Exception: "
                  << exception.what()
//                  << ": "
 //                 << exception.get_message()
                  << std::endl;
    }

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
                      STAR_MASS,
                      PLANET_MASS,
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
