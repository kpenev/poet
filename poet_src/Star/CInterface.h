/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * Star library.
 *
 * \ingroup Star_group
 */

#include "EvolvingStar.h"
#include "../StellarEvolution/CInterface.h"
#include "../Evolve/DissipationQuantities.h"

extern "C" {
    ///Identifier for not differentiating the phase lag.
    extern const int NO_DERIV;

    ///Identifier for differentiating the phase lag w.r.t. age.
    extern const int AGE_DERIV;

    ///Identifier for differentiating the phase lag w.r.t. spin freuqency.
    extern const int SPIN_FREQUENCY_DERIV;

    ///Identifier for differentiating the phase lag w.r.t. orbital freuqency.
    extern const int ORBITAL_FREQUENCY_DERIV;

    ///Opaque struct to cast to/from Star::InterpolatedEvolutionStar.
    struct EvolvingStar;

    ///Create a star to participate in the orbital evolution calculation.
    EvolvingStar *create_star(
        ///Mass of the star
        double mass,

        ///The metallicity ([Fe/H]) of the star
        double metallicity,

        ///The strength of the wind.
        double wind_strength,

        ///The frequency at which the wind loss saturates in rad/day.
        double wind_saturation_frequency,

        ///The timescale for differential rotation coupling.
        double diff_rot_coupling_timescale,

        ///A StellarEvolution interpolator.
        const MESAInterpolator *interpolator
    );

    ///Destroy a previously created star.
    void destroy_star(
        ///The star to destroy. Must have previously been created using
        ///create_star.
        EvolvingStar *star
    );

    ///Set the dissipative properties of one of the zones of a star.
    void set_dissipation(
        ///The star to set the dissipation for.
        EvolvingStar *star,

        ///Which zone to set the dissiaption for (0 - envelope, 1 - core).
        unsigned zone_index,

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

    ///Tell the star to detect its wind saturation state per its current
    ///configuration.
    void detect_stellar_wind_saturation(
        ///The star to detect the saturation state of.
        EvolvingStar *star
    );

    ///\brief Prepare the zone quantities for interpolation around the
    ///given age.
    ///
    ///See EvolvingStar::select_interpolation_region() for details.
    void select_interpolation_region(
        ///The star to select the interpolation region for.
        const EvolvingStar *star,

        ///The age around which interpolation should work.
        double age
    );

    ///\brief See Evolve::BrokenPowerlawPhaseLagZone::modified_phase_lag for
    ///details.
    double modified_phase_lag(
        ///The star to get the modified phase lag for.
        const EvolvingStar *star,

        ///The index of the zone whose modified phase lag to return.
        unsigned zone_index,

        ///The multiplier of the orbital frequency in the
        ///expression for the forcing frequency.
        int orbital_frequency_multiplier,

        ///The multiplier of the spin frequency in the
        ///expression for the forcing frequency.
        int spin_frequency_multiplier,

        ///The current forcing frequency in rad/day.
        double forcing_frequency,

        ///The return value should be either the phase lag itself
        ///(NO_DERIV) or its derivative w.r.t. the specified quantity.
        int deriv,

        ///If the lag of a locked term is calculated this should be set
        ///to the lag assuming the spin frequency is just above the lock.
        ///Otherwise, leave untouched.
        double *above_lock_value
    );

    ///The age at which the core of a star forms.
    double core_formation_age(
        ///The star for which to return the core formation age.
        const EvolvingStar *star
    );

    ///The lifetime of a star (the maximum age at which it can be queried)
    double lifetime(
        ///The star whose lifetime to return.
        const EvolvingStar *star
    );

    ///The luminosity of a star at a given age.
    double luminosity(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        double age
    );

    ///The luminosity of a star at a given array of ages.
    void luminosity_array(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The ages at which to return the luminosity.
        const double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );

    ///The moment of inertia of the stellar core at a given age.
    double core_inertia(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        double age
    );

    ///The moment of inertia of the stellar core at a given array of ages.
    void core_inertia_array(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The ages at which to return the moment of inertia.
        const double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );

    ///The moment of inertia of the stellar envelope at a given age.
    double envelope_inertia(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        double age
    );

    ///\brief The moment of inertia of the stellar envelope at a given array
    ///of ages.
    void envelope_inertia_array(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        const double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );

    ///The radius of the star at a given age.
    double star_radius(
        ///The star whose radius to return.
        EvolvingStar *star,

        ///The age at which to return the stellar radius.
        double age
    );

    ///The radius of the star at an array of ages.
    void star_radius_array(
        ///The star whose radius to return.
        EvolvingStar *star,

        ///The ages at which to return the stellar radius.
        const double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );

} //End extern "C"
