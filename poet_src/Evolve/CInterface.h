/**\file
 * 
 * \brief Declare C-style functions for accessing the functionality of the
 * Evolve library.
 *
 * \ingroup Evolve_group
 */

#include "DiskPlanetSystem.h"
#include "OrbitSolver.h"
#include "../Star/CInterface.h"
#include "../Planet/CInterface.h"

extern "C" {
    ///\brief Evolution mode ID for when the surface rotation of one of the
    ///bodies is locked to a prescribed value.
    extern const int LOCKED_SURFACE_SPIN_EVOL_MODE;

    ///\brief Evolution mode ID for when the two bodies orbit each other.
    extern const int BINARY_EVOL_MODE;

    ///\brief Evolution mode ID for when there is only one body in the system
    ///(only its rotation evolves).
    extern const int SINGLE_EVOL_MODE;

    ///\brief Evolution mode ID used as the mode to transform to from all
    ///other modes when storing the computed evolution.
    extern const int TABULATION_EVOL_MODE;

    ///Not a number.
    extern const double NaN;

    ///Opaque struct to cast to/from Evolve::DiskPlanetSystem.
    struct StarPlanetSystem;

    ///Opaque struct to cast to/from Evolve::OrbitSolver.
    struct OrbitSolver;

    ///Read eccentricity expansion coefficients from a file.
    void read_eccentricity_expansion_coefficients(const char *filename);

    ///Create a binary system out of two bodies.
    StarPlanetSystem *create_star_planet_system(
        ///The first body in the system. Assumed to always be there, so
        ///for a star-planet system this should be the star.
        EvolvingStar *star, 

        ///The second body in the system, initially may not be there and
        ///later may be engulfed by the first body.
        LockedPlanet *planet,

        ///The semimajor axis of the orbit at which the secondary forms in
        /// \f$R_\odot\f$.
        double initial_semimajor,

        ///The eccentricity of the orbit at which the secondary forms.
        double initial_eccentricity,

        ///Inclination between surface zone of primary and initial orbit in
        ///radians.
        double initial_inclination,

        ///Frequency of the surface spin of the primary when disk is present
        ///in rad/day.
        double disk_lock_frequency,

        ///Age when disk dissipates in Gyrs.
        double disk_dissipation_age,

        ///Age when the secondary forms.
        double secondary_formation_age
    );

    ///Destroy a previously created binary system.
    void destroy_binary(
        ///The system to destroy.
        StarPlanetSystem *system
    );

    ///\brief Defines the orbit a planet is in.
    ///
    ///The inclinations and arguments of periapsis must be already set for
    ///all zones.
    void configure_planet(
        ///The body to configure.
        LockedPlanet *planet,

        ///The age to set the body to.
        double age,

        ///The mass of the second body in the system.
        double companion_mass,

        ///The semimajor axis of the orbit in \f$R_\odot\f$.
        double semimajor,

        ///The eccentricity of the orbit
        double eccentricity,

        ///The spin angular momenta of the non-locked zones of the body
        ///(outermost zone to innermost).
        const double *spin_angmom,

        ///The inclinations of the zones of the body (same order as 
        ///spin_angmom). If NULL, all inclinations are assumed zero.
        const double *inclination,

        ///The arguments of periapsis of the zones of the bodies (same
        ///order as spin_angmom). If NULL, all periapses are assumed
        ///zero.
        const double *periapsis,

        ///If true, the outermost zone's spin is assumed locked to a 
        ///disk and spin_angmom is assumed to start from the next zone.
        bool locked_surface,

        ///If true, the outermost zone's inclination is assumed to be
        ///zero and the inclination argument is assumed to start from the
        ///next zone.
        bool zero_outer_inclination,

        ///If true, the outermost zone's periapsis is assumed to be
        ///zero and the inclination argument is assumed to start from the
        ///next zone.
        bool zero_outer_periapsis
    );

    ///\brief Defines the orbit a star is in.
    ///
    ///The inclinations and arguments of periapsis must be already set for
    ///all zones.
    void configure_star(
        ///The body to configure.
        EvolvingStar *star,

        ///The age to set the body to.
        double age,

        ///The mass of the second body in the system.
        double companion_mass,

        ///The semimajor axis of the orbit in \f$R_\odot\f$.
        double semimajor,

        ///The eccentricity of the orbit
        double eccentricity,

        ///The spin angular momenta of the non-locked zones of the body
        ///(outermost zone to innermost).
        const double *spin_angmom,

        ///The inclinations of the zones of the body (same order as 
        ///spin_angmom). If NULL, all inclinations are assumed zero.
        const double *inclination,

        ///The arguments of periapsis of the zones of the bodies (same
        ///order as spin_angmom). If NULL, all periapses are assumed
        ///zero.
        const double *periapsis,

        ///If true, the outermost zone's spin is assumed locked to a 
        ///disk and spin_angmom is assumed to start from the next zone.
        bool locked_surface,

        ///If true, the outermost zone's inclination is assumed to be
        ///zero and the inclination argument is assumed to start from the
        ///next zone.
        bool zero_outer_inclination,

        ///If true, the outermost zone's periapsis is assumed to be
        ///zero and the inclination argument is assumed to start from the
        ///next zone.
        bool zero_outer_periapsis
    );

    ///Sets the current state of a system.
    void configure_system(
        ///The system to set the state of.
        StarPlanetSystem *system,

        ///The age to set the system to.
        double age,

        ///The semimajor axis of the orbit.
        double semimajor,

        ///The eccentricity of the orbit.
        double eccentricity,

        ///The spin angular momenta of the zones of the bodies (body 1 
        ///first, outermost zone to innermost, followed by body 2).
        const double *spin_angmom,

        ///The inclinations of the zones of the bodies (same order as
        ///spin_angmom). The surface zone inclination must be omitted for
        ///single body systems.
        const double *inclination,

        ///The arguments of periapsis of the zones of the bodies (same
        ///order as spin_angmom, but not including the surface zone of
        ///the first body).
        const double *periapsis,

        ///The evolution mode to assume. Must be one of the constants
        ///defined.
        int evolution_mode
    );

    ///Calculate the evolution of a previously configured binary system.
    OrbitSolver *evolve_system(
        ///The system to evolve.
        StarPlanetSystem *system,

        ///The age at which to stop the evolution in Gyrs. The starting age
        ///must be already set for the system.
        double final_age,

        ///The maximum size of the time step allowed in Gyrs.
        double max_time_step,

        ///The precision to require of the solution.
        double precision,

        ///Ages at which the evolution must stop precisely.
        double *required_ages,

        ///The number of required ages.
        unsigned num_required_ages
    );

    ///Destroy a solver created by evolve_system.
    void destroy_solver(
        ///The solver to destroy.
        OrbitSolver *solver
    );

    ///At how many points was the evolution saved.
    unsigned num_evolution_steps(
        ///The solver used to follow the evolution.
        OrbitSolver *solver
    );

    ///\brief Fill C-style arrays with the calculated evolution of a system.
    ///
    ///All return arrays must either be allocated to the correct size or be
    ///NULL. In the latter case, the corresponding quantity is not returned.
    void get_evolution(
        ///The solver which was used to calculate the orbital evolution.
        const OrbitSolver *solver,

        ///The system which was evolved.
        const StarPlanetSystem *system,

        ///The star which was evolved.
        const EvolvingStar *star,

        ///An array to fill with the ages at which the evolution was
        ///calculated.
        double *age,

        ///An array to fill with the semimajor axis at the saved evolution
        ///steps.
        double *semimajor,

        ///An array to fill with the orbital eccentricity at the saved 
        ///evolution steps.
        double *eccentricity,

        ///An array to fill with the angle between the stellar convective
        ///envelope and the orbital angular momenta.
        double *envelope_inclination,

        ///An array to fill with the angle between the stellar radiative core
        ///and the orbital angular momenta.
        double *core_inclination,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar convective envelope.
        double *envelope_periapsis,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar radiative core.
        double *core_periapsis,

        ///An array to fill with the angular momentum of the stellar
        ///convective envelope.
        double *envelope_angmom,

        ///An array to fill with the angular momentum of the stellar
        ///radiative core.
        double *core_angmom,

        ///An array to fill with the evolution mode of the system.
        int *evolution_mode,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *wind_saturation
    );

}//End Extern "C"
