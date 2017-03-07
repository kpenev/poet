/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h.
 *
 * \ingroup Evolve_group
 */

#include "CInterface.h"

const int LOCKED_SURFACE_SPIN_EVOL_MODE = Core::LOCKED_SURFACE_SPIN;
const int BINARY_EVOL_MODE = Core::BINARY;
const int SINGLE_EVOL_MODE = Core::SINGLE;
const int TABULATION_EVOL_MODE = Core::TABULATION;
const double NaN = Core::NaN;

DiskPlanetSystem *create_disk_planet_system(DissipatingBody *primary, 
                                            DissipatingBody *secondary,
                                            double initial_semimajor,
                                            double initial_eccentricity,
                                            double initial_inclination,
                                            double disk_lock_frequency,
                                            double disk_dissipation_age,
                                            double secondary_formation_age)
{
    return reinterpret_cast<DiskPlanetSystem*>(
        new Evolve::DiskPlanetSystem(
            *reinterpret_cast<Evolve::DissipatingBody*>(primary),
            *reinterpret_cast<Evolve::DissipatingBody*>(secondary),
            initial_semimajor,
            initial_eccentricity,
            initial_inclination,
            disk_lock_frequency,
            disk_dissipation_age,
            std::max(secondary_formation_age,
                     disk_dissipation_age)
        )
    );
}

void destroy_binary(DiskPlanetSystem *system)
{
    delete reinterpret_cast<Evolve::DiskPlanetSystem*>(system);
}

void configure_body(DissipatingBody *body,
                    double age,
                    double companion_mass,
                    double semimajor,
                    double eccentricity,
                    const double *spin_angmom,
                    const double *inclination,
                    const double *periapsis,
                    bool locked_surface,
                    bool zero_outer_inclination,
                    bool zero_outer_periapsis)
{
    reinterpret_cast<Evolve::DissipatingBody*>(body)->configure(
        age,
        companion_mass,
        semimajor,
        eccentricity,
        spin_angmom,
        inclination,
        periapsis,
        locked_surface,
        zero_outer_inclination,
        zero_outer_periapsis
    );
}

void configure_system(DiskPlanetSystem *system,
                      double age,
                      double semimajor,
                      double eccentricity,
                      const double *spin_angmom,
                      const double *inclination,
                      const double *periapsis,
                      int evolution_mode)
{
    assert(evolution_mode < TABULATION_EVOL_MODE);
    reinterpret_cast<Evolve::DiskPlanetSystem*>(system)->configure(
        age,
        semimajor,
        eccentricity,
        spin_angmom,
        inclination,
        periapsis,
        static_cast<Core::EvolModeType>(evolution_mode)
    );
}

OrbitSolver *evolve_system(DiskPlanetSystem *system,
                           double final_age,
                           double max_time_step,
                           double precision,
                           double *required_ages,
                           unsigned num_required_ages)
{
    Evolve::OrbitSolver *solver = new Evolve::OrbitSolver(final_age,
                                                          precision);
#ifndef DEBUG
    try {
#endif
        (*solver)(
            *reinterpret_cast<Evolve::DiskPlanetSystem*>(system),
            max_time_step,
            std::list<double>(required_ages,
                              required_ages + num_required_ages)
        );
#ifndef DEBUG
    } catch(std::exception)
    {
    }
#endif
    return reinterpret_cast<OrbitSolver*>(solver);
}
