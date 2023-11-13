/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h.
 *
 * \ingroup Evolve_group
 */

#define BUILDING_LIBRARY
#include "CInterface.h"

const int LOCKED_SURFACE_SPIN_EVOL_MODE = Core::LOCKED_SURFACE_SPIN;
const int BINARY_EVOL_MODE = Core::BINARY;
const int SINGLE_EVOL_MODE = Core::SINGLE;
const int TABULATION_EVOL_MODE = Core::TABULATION;
const double NaN = Core::NaN;

void prepare_eccentricity_expansion(const char *filename,
                                    double precision,
                                    bool pre_load,
                                    bool disable_precision_fail)
{
    Evolve::TidalPotentialTerms::prepare(
        filename,
        precision,
        pre_load,
        disable_precision_fail
    );
}

void set_zone_dissipation(BrokenPowerlawPhaseLagZone *zone,
                          unsigned num_tidal_frequency_breaks,
                          unsigned num_spin_frequency_breaks,
                          unsigned num_age_breaks,
                          double *tidal_frequency_breaks,
                          double *spin_frequency_breaks,
                          double *tidal_frequency_powers,
                          double *spin_frequency_powers,
                          double *age_breaks,
                          double *reference_phase_lags,
                          double inertial_mode_enhancement,
                          double inertial_mode_sharpness)
{
    Evolve::BrokenPowerlawPhaseLagZone *real_zone =
        reinterpret_cast<Evolve::BrokenPowerlawPhaseLagZone*>(zone);
    std::cerr << "Defining zone dissipation" << std::endl;
    real_zone->setup(
        (
            num_tidal_frequency_breaks
            ? std::vector<double>(
                tidal_frequency_breaks,
                tidal_frequency_breaks + num_tidal_frequency_breaks
            )
            : std::vector<double>()
        ),
        (
            num_spin_frequency_breaks
            ? std::vector<double>(
                spin_frequency_breaks,
                spin_frequency_breaks + num_spin_frequency_breaks
            )
            : std::vector<double>()
        ),
        std::vector<double>(
            tidal_frequency_powers,
            tidal_frequency_powers + num_tidal_frequency_breaks + 1
        ),
        std::vector<double>(
            spin_frequency_powers,
            spin_frequency_powers + num_spin_frequency_breaks + 1
        ),
        std::vector<double>(
            reference_phase_lags,
            reference_phase_lags + num_age_breaks + 1
        ),
        inertial_mode_enhancement,
        inertial_mode_sharpness,
        (
            num_age_breaks
            ? std::vector<double>(
                age_breaks,
                age_breaks + num_age_breaks
            )
            : std::vector<double>()
        )
    );

}

void set_single_term_zone_dissipation(SingleTermZone *zone,
                                      int orbital_frequency_multiplier,
                                      int spin_frequency_multiplier,
                                      double phase_lag)
{
    Evolve::SingleTermZone *real_zone =
        reinterpret_cast<Evolve::SingleTermZone*>(zone);
    real_zone->setup(orbital_frequency_multiplier,
                     spin_frequency_multiplier,
                     phase_lag);
}

DiskBinarySystem *create_star_planet_system(EvolvingStar *star,
                                            CPlanet *planet,
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
            *reinterpret_cast<Planet::Planet*>(planet),
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

DiskBinarySystem *create_planet_planet_system(CPlanet *primary,
                                              CPlanet *secondary,
                                              double initial_semimajor,
                                              double initial_eccentricity,
                                              double initial_inclination,
                                              double disk_lock_frequency,
                                              double disk_dissipation_age)
{
    return reinterpret_cast<DiskBinarySystem*>(
        new Evolve::DiskBinarySystem(
            *reinterpret_cast<Planet::Planet*>(primary),
            *reinterpret_cast<Planet::Planet*>(secondary),
            initial_semimajor,
            initial_eccentricity,
            initial_inclination,
            disk_lock_frequency,
            disk_dissipation_age,
            disk_dissipation_age
        )
    );
}

DiskBinarySystem *create_single_term_system(
    CSingleTermNonEvolvingBody *primary,
    CSingleTermNonEvolvingBody *secondary,
    double initial_semimajor,
    double initial_eccentricity,
    double initial_inclination,
    double disk_lock_frequency,
    double disk_dissipation_age
)
{
    return reinterpret_cast<DiskBinarySystem*>(
        new Evolve::DiskBinarySystem(
            *reinterpret_cast<SingleTermNonEvolvingBody::SingleTermNonEvolvingBody*>(primary),
            *reinterpret_cast<SingleTermNonEvolvingBody::SingleTermNonEvolvingBody*>(secondary),
            initial_semimajor,
            initial_eccentricity,
            initial_inclination,
            disk_lock_frequency,
            disk_dissipation_age,
            disk_dissipation_age
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
#ifndef NDEBUG
    std::cerr << "Finished c-interface configuring star."<< std::endl;
#endif

}

void configure_planet(CPlanet *planet,
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
    reinterpret_cast<Planet::Planet*>(planet)->configure(
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

void configure_single_term_body(CSingleTermNonEvolvingBody *body,
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
    reinterpret_cast<SingleTermNonEvolvingBody::SingleTermNonEvolvingBody*>(
        body
    )->configure(
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


void set_expansion_order(DiskBinarySystem *system,
                         unsigned order)
{
    reinterpret_cast<Evolve::DiskBinarySystem*>(system)->change_expansion_order(
        order
    );
}


void differential_equations(DiskBinarySystem *system,
                            double age,
                            const double *parameters,
                            Core::EvolModeType evolution_mode,
                            double *differential_equations)
{
    Evolve::DiskBinarySystem* real_system =  (
        reinterpret_cast<Evolve::DiskBinarySystem*>(system)
    );
    std::valarray<double> diff_eq_parameters;
    if(!parameters) {
        real_system->fill_orbit(diff_eq_parameters);
        parameters = &diff_eq_parameters[0];
    }
    real_system->differential_equations(
        age,
        parameters,
        evolution_mode,
        differential_equations
    );
}

OrbitSolver *evolve_system(DiskBinarySystem *system,
                           double final_age,
                           double max_time_step,
                           double precision,
                           double *required_ages,
                           unsigned num_required_ages,
                           bool print_progress,
                           double max_runtime,
                           unsigned max_time_steps)
{
	std::cerr.setf(std::ios_base::scientific);
	std::cerr.precision(16);
    Evolve::OrbitSolver *solver = new Evolve::OrbitSolver(final_age,
                                                          precision,
                                                          print_progress);
#ifdef NDEBUG
    try {
#endif
        (*solver)(
            *reinterpret_cast<Evolve::DiskBinarySystem*>(system),
            max_time_step,
            std::list<double>(required_ages,
                              required_ages + num_required_ages),
            max_runtime,
            max_time_steps
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
void get_star_evolution(const EvolvingStar *star_arg,
                        double *envelope_inclination,
                        double *core_inclination,
                        double *envelope_periapsis,
                        double *core_periapsis,
                        double *envelope_angmom,
                        double *core_angmom,
                        bool *wind_saturation,
                        double *envelope_inclination_rate,
                        double *core_inclination_rate,
                        double *envelope_periapsis_rate,
                        double *core_periapsis_rate,
                        double *envelope_angmom_rate,
                        double *core_angmom_rate)
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

    list_to_array(
        star->envelope().get_evolution_real(Evolve::INCLINATION_DERIV),
        envelope_inclination_rate
    );

    list_to_array(
        star->core().get_evolution_real(Evolve::INCLINATION_DERIV),
        core_inclination_rate
    );

    list_to_array(star->envelope().get_evolution_real(Evolve::PERIAPSIS_DERIV),
                  envelope_periapsis_rate);

    list_to_array(star->core().get_evolution_real(Evolve::PERIAPSIS_DERIV),
                  core_periapsis_rate);

    list_to_array(
        star->envelope().get_evolution_real(Evolve::ANGULAR_MOMENTUM_DERIV),
        envelope_angmom_rate
    );

    list_to_array(
        star->core().get_evolution_real(Evolve::ANGULAR_MOMENTUM_DERIV),
        core_angmom_rate
    );
}

///Fill the given array with the part of the evolution tracked by the planet.
void get_planet_evolution(const CPlanet *planet_arg,
                          double *inclination,
                          double *periapsis,
                          double *angmom,
                          double *inclination_rate,
                          double *periapsis_rate,
                          double *angmom_rate)
{
    const Planet::Planet *planet = reinterpret_cast<const Planet::Planet*>(
        planet_arg
    );

    list_to_array(
        planet->zone().get_evolution_real(Evolve::INCLINATION),
        inclination
    );
    list_to_array(
        planet->zone().get_evolution_real(Evolve::PERIAPSIS),
        periapsis
    );
    list_to_array(
        planet->zone().get_evolution_real(Evolve::ANGULAR_MOMENTUM),
        angmom
    );

    list_to_array(
        planet->zone().get_evolution_real(Evolve::INCLINATION_DERIV),
        inclination_rate
    );
    list_to_array(
        planet->zone().get_evolution_real(Evolve::PERIAPSIS_DERIV),
        periapsis_rate
    );
    list_to_array(
        planet->zone().get_evolution_real(Evolve::ANGULAR_MOMENTUM_DERIV),
        angmom_rate
    );
}

///\brief Fill the given arrays with the part of the evolution (the orbital
///state) tracked by the binary system.
void get_binary_evolution(const DiskBinarySystem *system_arg,
                          double *semimajor,
                          double *eccentricity,
                          double *semimajor_rate,
                          double *eccentricity_rate)
{
    const Evolve::DiskBinarySystem *system =
        reinterpret_cast<const Evolve::DiskBinarySystem*>(system_arg);

    list_to_array(system->semimajor_evolution(), semimajor);
    list_to_array(system->eccentricity_evolution(), eccentricity);

    list_to_array(system->semimajor_evolution_rate(), semimajor_rate);
    list_to_array(system->eccentricity_evolution_rate(), eccentricity_rate);
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
                               const CPlanet *planet,
                               double *age,
                               double *semimajor,
                               double *eccentricity,
                               double *envelope_inclination,
                               double *core_inclination,
                               double *envelope_periapsis,
                               double *core_periapsis,
                               double *envelope_angmom,
                               double *core_angmom,
                               double *planet_inclination,
                               double *planet_periapsis,
                               double *planet_angmom,
                               int *evolution_mode,
                               bool *wind_saturation,
                               double *semimajor_rate,
                               double *eccentricity_rate,
                               double *envelope_inclination_rate,
                               double *core_inclination_rate,
                               double *envelope_periapsis_rate,
                               double *core_periapsis_rate,
                               double *envelope_angmom_rate,
                               double *core_angmom_rate,
                               double *planet_inclination_rate,
                               double *planet_periapsis_rate,
                               double *planet_angmom_rate)
{

    get_solver_evolution(solver, age, evolution_mode);

    get_binary_evolution(system,
                         semimajor,
                         eccentricity,
                         semimajor_rate,
                         eccentricity_rate);

    get_star_evolution(star,
                       envelope_inclination,
                       core_inclination,
                       envelope_periapsis,
                       core_periapsis,
                       envelope_angmom,
                       core_angmom,
                       wind_saturation,
                       envelope_inclination_rate,
                       core_inclination_rate,
                       envelope_periapsis_rate,
                       core_periapsis_rate,
                       envelope_angmom_rate,
                       core_angmom_rate);

    get_planet_evolution(planet,
                         planet_inclination,
                         planet_periapsis,
                         planet_angmom,
                         planet_inclination_rate,
                         planet_periapsis_rate,
                         planet_angmom_rate);
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
                             bool *secondary_wind_saturation,
                             double *semimajor_rate,
                             double *eccentricity_rate,
                             double *primary_envelope_inclination_rate,
                             double *primary_core_inclination_rate,
                             double *primary_envelope_periapsis_rate,
                             double *primary_core_periapsis_rate,
                             double *primary_envelope_angmom_rate,
                             double *primary_core_angmom_rate,
                             double *secondary_envelope_inclination_rate,
                             double *secondary_core_inclination_rate,
                             double *secondary_envelope_periapsis_rate,
                             double *secondary_core_periapsis_rate,
                             double *secondary_envelope_angmom_rate,
                             double *secondary_core_angmom_rate)
{
    get_solver_evolution(solver, age, evolution_mode);

    get_binary_evolution(system,
                         semimajor,
                         eccentricity,
                         semimajor_rate,
                         eccentricity_rate);

    get_star_evolution(primary,
                       primary_envelope_inclination,
                       primary_core_inclination,
                       primary_envelope_periapsis,
                       primary_core_periapsis,
                       primary_envelope_angmom,
                       primary_core_angmom,
                       primary_wind_saturation,
                       primary_envelope_inclination_rate,
                       primary_core_inclination_rate,
                       primary_envelope_periapsis_rate,
                       primary_core_periapsis_rate,
                       primary_envelope_angmom_rate,
                       primary_core_angmom_rate);

    get_star_evolution(secondary,
                       secondary_envelope_inclination,
                       secondary_core_inclination,
                       secondary_envelope_periapsis,
                       secondary_core_periapsis,
                       secondary_envelope_angmom,
                       secondary_core_angmom,
                       secondary_wind_saturation,
                       secondary_envelope_inclination_rate,
                       secondary_core_inclination_rate,
                       secondary_envelope_periapsis_rate,
                       secondary_core_periapsis_rate,
                       secondary_envelope_angmom_rate,
                       secondary_core_angmom_rate);
}

void get_planet_planet_evolution(const OrbitSolver *solver,
                                 const DiskBinarySystem *system,
                                 const CPlanet *primary,
                                 const CPlanet *secondary,
                                 double *age,
                                 double *semimajor,
                                 double *eccentricity,
                                 double *primary_inclination,
                                 double *primary_periapsis,
                                 double *primary_angmom,
                                 double *secondary_inclination,
                                 double *secondary_periapsis,
                                 double *secondary_angmom,
                                 int *evolution_mode,
                                 double *semimajor_rate,
                                 double *eccentricity_rate,
                                 double *primary_inclination_rate,
                                 double *primary_periapsis_rate,
                                 double *primary_angmom_rate,
                                 double *secondary_inclination_rate,
                                 double *secondary_periapsis_rate,
                                 double *secondary_angmom_rate)
{

    get_solver_evolution(solver, age, evolution_mode);

    get_binary_evolution(system,
                         semimajor,
                         eccentricity,
                         semimajor_rate,
                         eccentricity_rate);

    get_planet_evolution(primary,
                         primary_inclination,
                         primary_periapsis,
                         primary_angmom,
                         primary_inclination_rate,
                         primary_periapsis_rate,
                         primary_angmom_rate);

    get_planet_evolution(secondary,
                         secondary_inclination,
                         secondary_periapsis,
                         secondary_angmom,
                         secondary_inclination_rate,
                         secondary_periapsis_rate,
                         secondary_angmom_rate);
}

///\brief Fill the given pointers with the state of the given star at the end
///of the evolution.
void get_star_final_state(const EvolvingStar *star_arg,
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

    if(wind_saturation)
        *wind_saturation = star->wind_saturation_evolution().back();
}

///\brief Fill the given pointers with the state of the given planet at the end
///of the evolution.
void get_planet_final_state(const CPlanet *planet_arg,
                          double *inclination,
                          double *periapsis,
                          double *angmom)
{
    const Planet::Planet *planet = reinterpret_cast<const Planet::Planet*>(
        planet_arg
    );

    if(inclination)
        *inclination = planet->zone().get_evolution_real(
            Evolve::INCLINATION
        ).back();

    if(periapsis)
        *periapsis = planet->zone().get_evolution_real(
            Evolve::PERIAPSIS
        ).back();

    if(angmom)
        *angmom = planet->zone().get_evolution_real(
            Evolve::ANGULAR_MOMENTUM
        ).back();
}

///\brief Fill the given pointers with the final orbital state of a
///previously evolved system.
void get_binary_final_state(const DiskBinarySystem *system_arg,
                            double *semimajor,
                            double *eccentricity)
{
    const Evolve::DiskBinarySystem *system =
        reinterpret_cast<const Evolve::DiskBinarySystem*>(system_arg);

    if(semimajor)
        *semimajor = system->semimajor_evolution().back();

    if(eccentricity)
        *eccentricity = system->eccentricity_evolution().back();
}

///\brief Fill the given pointers with the final state of an orbit solver
///used to calculate an evolution.
void get_solver_final_state(const OrbitSolver *solver_arg,
                            double *age,
                            int *evolution_mode)
{
    const Evolve::OrbitSolver *solver =
        reinterpret_cast<const Evolve::OrbitSolver*>(solver_arg);

    if(age)
        *age = solver->evolution_ages().back();

    if(evolution_mode)
        *evolution_mode = solver->mode_evolution().back();
}

void get_star_planet_final_state(const OrbitSolver *solver,
                                 const DiskBinarySystem *system,
                                 const EvolvingStar *star,
                                 const CPlanet *planet,
                                 double *age,
                                 double *semimajor,
                                 double *eccentricity,
                                 double *envelope_inclination,
                                 double *core_inclination,
                                 double *envelope_periapsis,
                                 double *core_periapsis,
                                 double *envelope_angmom,
                                 double *core_angmom,
                                 double *planet_inclination,
                                 double *planet_periapsis,
                                 double *planet_angmom,
                                 int *evolution_mode,
                                 bool *wind_saturation)
{
    get_solver_final_state(solver, age, evolution_mode);

    get_binary_final_state(system, semimajor, eccentricity);

    get_star_final_state(star,
                         envelope_inclination,
                         core_inclination,
                         envelope_periapsis,
                         core_periapsis,
                         envelope_angmom,
                         core_angmom,
                         wind_saturation);

    get_planet_final_state(planet,
                           planet_inclination,
                           planet_periapsis,
                           planet_angmom);
}

void get_star_star_final_state(const OrbitSolver *solver,
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
    get_solver_final_state(solver, age, evolution_mode);

    get_binary_final_state(system, semimajor, eccentricity);

    get_star_final_state(primary,
                         primary_envelope_inclination,
                         primary_core_inclination,
                         primary_envelope_periapsis,
                         primary_core_periapsis,
                         primary_envelope_angmom,
                         primary_core_angmom,
                         primary_wind_saturation);

    get_star_final_state(secondary,
                         secondary_envelope_inclination,
                         secondary_core_inclination,
                         secondary_envelope_periapsis,
                         secondary_core_periapsis,
                         secondary_envelope_angmom,
                         secondary_core_angmom,
                         secondary_wind_saturation);
}

double get_expansion_coeff_precision(int m, int s)
{
    return Evolve::TidalPotentialTerms::expansion_coefficient_evaluator(
    ).interp_precision(
        m,
        s
    );
}

double evaluate_expansion_coeff(int m,
                                int s,
                                double e,
                                bool deriv)
{
    return Evolve::TidalPotentialTerms::expansion_coefficient_evaluator()(
        m,
        s,
        e,
        deriv
    );
}

void destroy_expansion_coef(
    const EccentricityExpansionCoefficients *expansion_arg
)
{
    delete reinterpret_cast<const Evolve::EccentricityExpansionCoefficients*>(
        expansion_arg
    );
}
