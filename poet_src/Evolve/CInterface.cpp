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

void read_eccentricity_expansion_coefficients(const char *filename)
{
    Evolve::DissipatingZone::read_eccentricity_expansion(filename);
}

DiskBinarySystem *create_star_planet_system(EvolvingStar *star, 
                                            LockedPlanet *planet,
                                            double initial_semimajor,
                                            double initial_eccentricity,
                                            double initial_inclination,
                                            double disk_lock_frequency,
                                            double disk_dissipation_age,
                                            double secondary_formation_age)
{
    return reinterpret_cast<DiskBinarySystem*>(
        new Evolve::DiskBinarySystem(
            *reinterpret_cast<Star::InterpolatedEvolutionStar*>(star),
            *reinterpret_cast<Planet::LockedPlanet*>(planet),
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

DiskBinarySystem *create_star_star_system(EvolvingStar *primary, 
                                          EvolvingStar *secondary,
                                          double initial_semimajor,
                                          double initial_eccentricity,
                                          double initial_inclination,
                                          double disk_lock_frequency,
                                          double disk_dissipation_age,
                                          double secondary_formation_age)
{
    return reinterpret_cast<DiskBinarySystem*>(
        new Evolve::DiskBinarySystem(
            *reinterpret_cast<Star::InterpolatedEvolutionStar*>(primary),
            *reinterpret_cast<Star::InterpolatedEvolutionStar*>(secondary),
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

void destroy_binary(DiskBinarySystem *system)
{
    delete reinterpret_cast<Evolve::DiskBinarySystem*>(system);
}

void configure_star(EvolvingStar *star,
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
    reinterpret_cast<Star::InterpolatedEvolutionStar*>(star)->configure(
        true,
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

void configure_planet(LockedPlanet *planet,
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
    reinterpret_cast<Planet::LockedPlanet*>(planet)->configure(
        true,
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

void configure_system(DiskBinarySystem *system,
                      double age,
                      double semimajor,
                      double eccentricity,
                      const double *spin_angmom,
                      const double *inclination,
                      const double *periapsis,
                      int evolution_mode)
{
    assert(evolution_mode < TABULATION_EVOL_MODE);
    reinterpret_cast<Evolve::DiskBinarySystem*>(system)->configure(
        true,
        age,
        semimajor,
        eccentricity,
        spin_angmom,
        inclination,
        periapsis,
        static_cast<Core::EvolModeType>(evolution_mode)
    );
}

OrbitSolver *evolve_system(DiskBinarySystem *system,
                           double final_age,
                           double max_time_step,
                           double precision,
                           double *required_ages,
                           unsigned num_required_ages)
{
    Evolve::OrbitSolver *solver = new Evolve::OrbitSolver(final_age,
                                                          precision);
#ifdef NDEBUG
    try {
#endif
        (*solver)(
            *reinterpret_cast<Evolve::DiskBinarySystem*>(system),
            max_time_step,
            std::list<double>(required_ages,
                              required_ages + num_required_ages)
        );
#ifdef NDEBUG
    } catch(std::exception)
    {
    }
#endif
    return reinterpret_cast<OrbitSolver*>(solver);
}

void destroy_solver(OrbitSolver *solver)
{
    delete reinterpret_cast<Evolve::OrbitSolver*>(solver);
}

unsigned num_evolution_steps(OrbitSolver *solver)
{
    return reinterpret_cast<Evolve::OrbitSolver*>(
        solver
    )->evolution_ages().size();
}

template<typename T>
inline void list_to_array(const std::list<T> &source, T *destination)
{
    if(destination)
        std::copy(source.begin(), source.end(), destination);
}

///Fill the given arrays with the part of the evolution tracked by the star.
void get_star_evolution(const EvolvingStar *star_arg
                        double *envelope_inclination,
                        double *core_inclination,
                        double *envelope_periapsis,
                        double *core_periapsis,
                        double *envelope_angmom,
                        double *core_angmom,
                        bool *wind_saturation)
{
    const Star::InterpolatedEvolutionStar *star = 
        reinterpret_cast<const Star::InterpolatedEvolutionStar*>(star_arg);
    list_to_array(star->envelope().get_evolution_real(Evolve::INCLINATION),
                  envelope_inclination);

    list_to_array(star->core().get_evolution_real(Evolve::INCLINATION),
                  core_inclination);

    list_to_array(star->envelope().get_evolution_real(Evolve::PERIAPSIS),
                  envelope_periapsis);

    list_to_array(star->core().get_evolution_real(Evolve::PERIAPSIS),
                  core_periapsis);

    list_to_array(
        star->envelope().get_evolution_real(Evolve::ANGULAR_MOMENTUM),
        envelope_angmom
    );

    list_to_array(star->core().get_evolution_real(Evolve::ANGULAR_MOMENTUM),
                  core_angmom);

    list_to_array(star->wind_saturation_evolution(), wind_saturation);
}

///\brief Fill the given arrays with the part of the evolution (the orbital
///state) tracked by the binary system.
void get_binary_evolution(const DiskBinarySystem *system_arg,
                          double *semimajor,
                          double *eccentricity)
{
    const Evolve::DiskBinarySystem *system =
        reinterpret_cast<const Evolve::DiskBinarySystem*>(system_arg);

    list_to_array(system->semimajor_evolution(), semimajor);

    list_to_array(system->eccentricity_evolution(), eccentricity);
}

///\brief Fill the given arrays with the part of the evolution tracked by the
///orbit solver.
void get_solver_evolution(const OrbitSolver *solver_arg,
                          double *age,
                          int *evolution_mode)
{
    const Evolve::OrbitSolver *solver =
        reinterpret_cast<const Evolve::OrbitSolver*>(solver_arg);

    list_to_array(solver->evolution_ages(), age);

    if(evolution_mode)
        std::copy(solver->mode_evolution().begin(),
                  solver->mode_evolution().end(),
                  evolution_mode);
}

void get_star_planet_evolution(const OrbitSolver *solver,
                               const DiskBinarySystem *system,
                               const EvolvingStar *star,
                               double *age,
                               double *semimajor,
                               double *eccentricity,
                               double *envelope_inclination,
                               double *core_inclination,
                               double *envelope_periapsis,
                               double *core_periapsis,
                               double *envelope_angmom,
                               double *core_angmom,
                               int *evolution_mode,
                               bool *wind_saturation)
{

    get_solver_evolution(solver, age, evolution_mode);

    get_binary_evolution(system, semimajor, eccentricity);

    get_star_evolution(star,
                       envelope_inclination,
                       core_inclination,
                       envelope_periapsis,
                       core_periapsis,
                       envelope_angmom,
                       core_angmom,
                       wind_saturation);

}

void get_star_star_evolution(const OrbitSolver *solver,
                             const DiskBinarySystem *system,
                             const EvolvingStar *primary,
                             const EvolvingStar *secondary,
                             double *age,
                             double *semimajor,
                             double *eccentricity,
                             double *primary_envelope_inclination,
                             double *primary_core_inclination,
                             double *primary_envelope_periapsis,
                             double *primary_core_periapsis,
                             double *primary_envelope_angmom,
                             double *primary_core_angmom,
                             double *secondary_envelope_inclination,
                             double *secondary_core_inclination,
                             double *secondary_envelope_periapsis,
                             double *secondary_core_periapsis,
                             double *secondary_envelope_angmom,
                             double *secondary_core_angmom,
                             int *evolution_mode,
                             bool *primary_wind_saturation,
                             bool *secondary_wind_saturation)
{
    get_solver_evolution(solver, age, evolution_mode);

    get_binary_evolution(system, semimajor, eccentricity);

    get_star_evolution(primary,
                       primary_envelope_inclination,
                       primary_core_inclination,
                       primary_envelope_periapsis,
                       primary_core_periapsis,
                       primary_envelope_angmom,
                       primary_core_angmom,
                       primary_wind_saturation);

    get_star_evolution(secondary,
                       secondary_envelope_inclination,
                       secondary_core_inclination,
                       secondary_envelope_periapsis,
                       secondary_core_periapsis,
                       secondary_envelope_angmom,
                       secondary_core_angmom,
                       secondary_wind_saturation);
}

void get_final_state(const OrbitSolver *solver_arg,
                     const DiskBinarySystem *system_arg,
                     const EvolvingStar *star_arg,
                     double *age,
                     double *semimajor,
                     double *eccentricity,
                     double *envelope_inclination,
                     double *core_inclination,
                     double *envelope_periapsis,
                     double *core_periapsis,
                     double *envelope_angmom,
                     double *core_angmom,
                     int *evolution_mode,
                     bool *wind_saturation)
{
    const Evolve::OrbitSolver *solver =
        reinterpret_cast<const Evolve::OrbitSolver*>(solver_arg);
    const Evolve::DiskBinarySystem *system =
        reinterpret_cast<const Evolve::DiskBinarySystem*>(system_arg);
    const Star::InterpolatedEvolutionStar *star = 
        reinterpret_cast<const Star::InterpolatedEvolutionStar*>(star_arg);

    if(age)
        *age = solver->evolution_ages().back();

    if(semimajor)
        *semimajor = system->semimajor_evolution().back();

    if(eccentricity)
        *eccentricity = system->eccentricity_evolution().back();

    if(envelope_inclination)
        *envelope_inclination = (
            star->envelope().get_evolution_real(Evolve::INCLINATION).back()
        );

    if(core_inclination)
        *core_inclination = star->core().get_evolution_real(
            Evolve::INCLINATION
        ).back();

    if(envelope_periapsis)
        *envelope_periapsis = star->envelope().get_evolution_real(
            Evolve::PERIAPSIS
        ).back();

    if(core_periapsis)
        *core_periapsis = star->core().get_evolution_real(
            Evolve::PERIAPSIS
        ).back();

    if(envelope_angmom)
        *envelope_angmom = star->envelope().get_evolution_real(
            Evolve::ANGULAR_MOMENTUM
        ).back();

    if(core_angmom)
        *core_angmom = star->core().get_evolution_real(
            Evolve::ANGULAR_MOMENTUM
        ).back();

    if(evolution_mode)
        *evolution_mode = solver->mode_evolution().back();

    if(wind_saturation)
        *wind_saturation = star->wind_saturation_evolution().back();
}

