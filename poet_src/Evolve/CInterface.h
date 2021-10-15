/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * Evolve library.
 *
 * \ingroup Evolve_group
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "DiskBinarySystem.h"
#include "OrbitSolver.h"
#include "../Star/CInterface.h"
#include "../Planet/CInterface.h"

extern "C" {
    ///\brief Evolution mode ID for when the surface rotation of one of the
    ///bodies is locked to a prescribed value.
    LIB_PUBLIC extern const int LOCKED_SURFACE_SPIN_EVOL_MODE;

    ///\brief Evolution mode ID for when the two bodies orbit each other.
    LIB_PUBLIC extern const int BINARY_EVOL_MODE;

    ///\brief Evolution mode ID for when there is only one body in the system
    ///(only its rotation evolves).
    LIB_PUBLIC extern const int SINGLE_EVOL_MODE;

    ///\brief Evolution mode ID used as the mode to transform to from all
    ///other modes when storing the computed evolution.
    LIB_PUBLIC extern const int TABULATION_EVOL_MODE;

    ///Not a number.
    LIB_PUBLIC extern const double NaN;

    ///Opaque struct to cast to/from Evolve::DiskBinarySystem.
    struct LIB_PUBLIC DiskBinarySystem;

    ///Opaque struct to cast to/from Evolve::OrbitSolver.
    struct LIB_PUBLIC OrbitSolver;

    ///Opaque struct to cast to/from Evolve::BrokenPowerlawPhasLagZone
    struct LIB_PUBLIC BrokenPowerlawPhaseLagZone;
    
    ///Opaque struct to cast to/from Evolve::EccentricityExpansionCoefficients
    struct LIB_PUBLIC EccentricityExpansionCoefficients;

    ///Read eccentricity expansion coefficients from a file.
    LIB_PUBLIC void read_eccentricity_expansion_coefficients(
        const char *filename
    );

    LIB_PUBLIC void set_zone_dissipation(
        ///The zone to set the dissipation of.
        BrokenPowerlawPhaseLagZone *zone,

        ///The number of breaks in the tidal frequency dependence.
        unsigned num_tidal_frequency_breaks,

        ///The number of breaks in the spin frequency dependence.
        unsigned num_spin_frequency_breaks,

        ///The locations of the breaks in tidal frequency in rad/day.
        ///Entries should be sorted.
        double *tidal_frequency_breaks,

        ///The locations of the breaks in spin frequency in rad/day.
        ///Entries should be sorted.
        double *spin_frequency_breaks,

        ///The powerlaw indices for the tidal frequency dependence.
        ///Should be indexed in the same order as tidal_frequency_breaks,
        ///but must contain an additional starting entry for the powerlaw
        ///index before the first break.
        double *tidal_frequency_powers,

        ///The powerlaw indices for the spin frequency dependence.
        ///Should be indexed in the same order as spin_frequency_breaks,
        ///but must contain an additional starting entry for the powerlaw
        ///index before the first break.
        double *spin_frequency_powers,

        ///The phase lag at the first tidal and first spin frequency break.
        ///The rest are calculated by imposing continuity.
        double reference_phase_lag
    );

    ///Create a binary system out of a star and a planet.
    LIB_PUBLIC DiskBinarySystem *create_star_planet_system(
        ///The first body in the system. Assumed to always be there, so
        ///for a star-planet system this should be the star.
        EvolvingStar *star,

        ///The second body in the system, initially may not be there and
        ///later may be engulfed by the first body.
        CPlanet *planet,

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

    ///Create a binary system out of two stars.
    LIB_PUBLIC DiskBinarySystem *create_star_star_system(
        ///The first body in the system. Assumed to always be there, so
        ///for a star-planet system this should be the star.
        EvolvingStar *primary,

        ///The second body in the system, initially may not be there and
        ///later may be engulfed by the first body.
        EvolvingStar *secondary,

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
    LIB_PUBLIC void destroy_binary(
        ///The system to destroy.
        DiskBinarySystem *system
    );

    ///\brief Defines the orbit a planet is in.
    ///
    ///The inclinations and arguments of periapsis must be already set for
    ///all zones.
    LIB_PUBLIC void configure_planet(
        ///The body to configure.
        CPlanet *planet,

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
    LIB_PUBLIC void configure_star(
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
    LIB_PUBLIC void configure_system(
        ///The system to set the state of.
        DiskBinarySystem *system,

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
    LIB_PUBLIC OrbitSolver *evolve_system(
        ///The system to evolve.
        DiskBinarySystem *system,

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
        unsigned num_required_ages,

        ///See print_progress argument to Evolve::OrbitSolver::OrbitSolver()
        bool print_progress
    );

    ///Destroy a solver created by evolve_system.
    LIB_PUBLIC void destroy_solver(
        ///The solver to destroy.
        OrbitSolver *solver
    );

    ///At how many points was the evolution saved.
    LIB_PUBLIC unsigned num_evolution_steps(
        ///The solver used to follow the evolution.
        OrbitSolver *solver
    );

    ///\brief Fill C-style arrays with the calculated evolution of a
    ///star-planet system.
    ///
    ///All return arrays must either be allocated to the correct size or be
    ///NULL. In the latter case, the corresponding quantity is not returned.
    LIB_PUBLIC void get_star_planet_evolution(
        ///The solver which was used to calculate the orbital evolution.
        const OrbitSolver *solver,

        ///The system which was evolved.
        const DiskBinarySystem *system,

        ///The star which was evolved.
        const EvolvingStar *star,

        ///The planet which was evolved,
        const CPlanet *planet,

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

        ///An array to fill wit the angle between the planet's and the
        ///orbital angular momentum (pass NULL if the planet is not
        ///dissipative).
        double *planet_inclination,

        ///An array to fill wit the  periapsis of the orbit in the equatorial
        ///plane of the planet (pass NULL if the planet is not dissipative).
        double *planet_periapsis,

        ///An array to fill wit the angular momentum of the planet (pass NULL if
        ///the planet is not dissipative).
        double *planet_angmom,

        ///An array to fill with the evolution mode of the system.
        int *evolution_mode,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *wind_saturation
    );

    ///\brief Fill C-style arrays with the calculated evolution of a
    ///binary star system.
    ///
    ///All return arrays must either be allocated to the correct size or be
    ///NULL. In the latter case, the corresponding quantity is not returned.
    LIB_PUBLIC void get_star_star_evolution(
        ///The solver which was used to calculate the orbital evolution.
        const OrbitSolver *solver,

        ///The system which was evolved.
        const DiskBinarySystem *system,

        ///The primary star which was evolved.
        const EvolvingStar *primary,

        ///The secondary star which was evolved.
        const EvolvingStar *secondary,

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
        double *primary_envelope_inclination,

        ///An array to fill with the angle between the stellar radiative core
        ///and the orbital angular momenta.
        double *primary_core_inclination,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar convective envelope.
        double *primary_envelope_periapsis,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar radiative core.
        double *primary_core_periapsis,

        ///An array to fill with the angular momentum of the stellar
        ///convective envelope.
        double *primary_envelope_angmom,

        ///An array to fill with the angular momentum of the stellar
        ///radiative core.
        double *primary_core_angmom,

        ///An array to fill with the angle between the stellar convective
        ///envelope and the orbital angular momenta.
        double *secondary_envelope_inclination,

        ///An array to fill with the angle between the stellar radiative core
        ///and the orbital angular momenta.
        double *secondary_core_inclination,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar convective envelope.
        double *secondary_envelope_periapsis,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar radiative core.
        double *secondary_core_periapsis,

        ///An array to fill with the angular momentum of the stellar
        ///convective envelope.
        double *secondary_envelope_angmom,

        ///An array to fill with the angular momentum of the stellar
        ///radiative core.
        double *secondary_core_angmom,

        ///An array to fill with the evolution mode of the system.
        int *evolution_mode,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *primary_wind_saturation,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *secondary_wind_saturation
    );

    ///\brief Fill destiantions with the calculated final state of a
    ///star-planet system.
    ///
    ///All return variables must either be allocated or be NULL. In the
    ///latter case, the corresponding quantity is not returned.
    LIB_PUBLIC void get_star_planet_final_state(
        ///The solver which was used to calculate the orbital evolution.
        const OrbitSolver *solver,

        ///The system which was evolved.
        const DiskBinarySystem *system,

        ///The star which was evolved.
        const EvolvingStar *star,

        ///The planet which was evolved.
        const CPlanet *planet,

        ///An array to fill with the ages at which the evolution was
        ///calculated.
        double *age,

        ///An array to fill with the semimajor axis at the saved evolution
        ///steps.
        double *semimajor,

        ///An array to fill with the orbital eccentricity at the saved
        ///evolution steps.
        double *eccentricity,

        ///To overwrite with the angle between the stellar convective
        ///envelope and the orbital angular momenta.
        double *envelope_inclination,

        ///To overwrite with the angle between the stellar radiative core
        ///and the orbital angular momenta.
        double *core_inclination,

        ///To overwrite with the periapsis of the orbit in the equatorial
        ///plane of the stellar convective envelope.
        double *envelope_periapsis,

        ///To overwrite with the periapsis of the orbit in the equatorial
        ///plane of the stellar radiative core.
        double *core_periapsis,

        ///To overwrite with the angular momentum of the stellar
        ///convective envelope.
        double *envelope_angmom,

        ///To overwrite with the angular momentum of the stellar
        ///radiative core.
        double *core_angmom,

        ///To overwrite with the angle between the planet's and the
        ///orbital angular momentum (pass NULL if the planet is not
        ///dissipative).
        double *planet_inclination,

        ///To overwrite with the  periapsis of the orbit in the equatorial
        ///plane of the planet (pass NULL if the planet is not dissipative).
        double *planet_periapsis,

        ///To overwrite with the angular momentum of the planet (pass NULL if
        ///the planet is not dissipative).
        double *planet_angmom,

        ///An array to fill with the evolution mode of the system.
        int *evolution_mode,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *wind_saturation
    );

    ///\brief Fill destiantions with the calculated final state of a
    ///binary star system.
    ///
    ///All return variables must either be allocated or be NULL. In the
    ///latter case, the corresponding quantity is not returned.
    LIB_PUBLIC void get_star_star_final_state(
        ///The solver which was used to calculate the orbital evolution.
        const OrbitSolver *solver,

        ///The system which was evolved.
        const DiskBinarySystem *system,

        ///The primary star in the system.
        const EvolvingStar *primary,

        ///The secondary star in the system.
        const EvolvingStar *secondary,

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
        double *primary_envelope_inclination,

        ///An array to fill with the angle between the stellar radiative core
        ///and the orbital angular momenta.
        double *primary_core_inclination,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar convective envelope.
        double *primary_envelope_periapsis,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar radiative core.
        double *primary_core_periapsis,

        ///An array to fill with the angular momentum of the stellar
        ///convective envelope.
        double *primary_envelope_angmom,

        ///An array to fill with the angular momentum of the stellar
        ///radiative core.
        double *primary_core_angmom,

        ///An array to fill with the angle between the stellar convective
        ///envelope and the orbital angular momenta.
        double *secondary_envelope_inclination,

        ///An array to fill with the angle between the stellar radiative core
        ///and the orbital angular momenta.
        double *secondary_core_inclination,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar convective envelope.
        double *secondary_envelope_periapsis,

        ///An array to fill with the periapsis of the orbit in the equatorial
        ///plane of the stellar radiative core.
        double *secondary_core_periapsis,

        ///An array to fill with the angular momentum of the stellar
        ///convective envelope.
        double *secondary_envelope_angmom,

        ///An array to fill with the angular momentum of the stellar
        ///radiative core.
        double *secondary_core_angmom,

        ///An array to fill with the evolution mode of the system.
        int *evolution_mode,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *primary_wind_saturation,

        ///An array to fill with a flag indicating whether the angular
        ///momentum loss due to stellar wind is in the satured state (true)
        ///or not (false).
        bool *secondary_wind_saturation
    );
    
    ///\brief doc
    LIB_PUBLIC EccentricityExpansionCoefficients *coeff_new();
    
    ///\brief doc
    LIB_PUBLIC void coeff_read(
        ///doc
        const EccentricityExpansionCoefficients *expansion_arg,
        
        ///The name of the file to read pre-tabulated coefficients from.
        const std::string &tabulated_pms_fname,

        ///Only reads the coefficients of the expansion which
        ///most closely achieves the indicated precision, erring on
        ///the side of being more precise.
        double *precision,
        
        ///Which loading style to use. Setting to true means we keep
        ///a 9 GB database in memory.
        bool *load_style
    );
    
    ///\brief doc
    LIB_PUBLIC unsigned coeff_max_e(
        ///doc
        const EccentricityExpansionCoefficients *expansion_arg
    );
    
    ///\brief doc
    LIB_PUBLIC double coeff_max_precision(
        ///doc
        const EccentricityExpansionCoefficients *expansion_arg,
        
        ///doc
        int *m,
        
        ///doc
        int *s
    );
    
    ///\brief doc
    LIB_PUBLIC double coeff_operator(
        ///doc
        const EccentricityExpansionCoefficients *expansion_arg,
        
        ///The first index (0 or +-2).
        int *m, 

        ///The second index.
        int *s,

        ///The value of the eccentricity to use.
        double *e,

        ///Prevision, the maximum eccentricity order to include in the Taylor
        ///series. Currently does nothing.
        unsigned *max_e_power,
        
        ///Previously, if true the result was differentiated w.r.t. to the
        ///eccentricity. Currently does nothing.
        bool *deriv
    );

}//End Extern "C"
