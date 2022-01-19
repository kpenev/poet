/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * Star library.
 *
 * \ingroup Star_group
 */


#include "../StellarEvolution/CInterface.h"
#include "../Core/SharedLibraryExportMacros.h"
#include "EvolvingStar.h"
#include "../Evolve/DissipationQuantities.h"


extern "C" {
    ///Identifier for not differentiating the phase lag.
    LIB_PUBLIC extern const int NO_DERIV;

    ///Identifier for differentiating the phase lag w.r.t. age.
    LIB_PUBLIC extern const int AGE_DERIV;

    ///Identifier for differentiating the phase lag w.r.t. spin freuqency.
    LIB_PUBLIC extern const int SPIN_FREQUENCY_DERIV;

    ///Identifier for differentiating the phase lag w.r.t. orbital freuqency.
    LIB_PUBLIC extern const int ORBITAL_FREQUENCY_DERIV;

    ///Opaque struct to cast to/from Star::InterpolatedEvolutionStar.
    struct LIB_PUBLIC EvolvingStar;

    ///Create a star to participate in the orbital evolution calculation.
    LIB_PUBLIC EvolvingStar *create_star(
        ///Mass of the star
        double mass,

        ///The [Fe/H] of the star
        double feh,

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
    LIB_PUBLIC void destroy_star(
        ///The star to destroy. Must have previously been created using
        ///create_star.
        EvolvingStar *star
    );

    ///Set the dissipative properties of one of the zones of a star.
    LIB_PUBLIC void set_star_dissipation(
        ///The star to set the dissipation for.
        EvolvingStar *star,

        ///Which zone to set the dissiaption for (0 - envelope, 1 - core).
        unsigned zone_index,

        ///See same name argument to set_zone_dissipation()
        unsigned num_tidal_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        unsigned num_spin_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        double *tidal_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        double *spin_frequency_breaks,

        ///See same name argument to set_zone_dissipation()
        double *tidal_frequency_powers,

        ///See same name argument to set_zone_dissipation()
        double *spin_frequency_powers,

        ///See same name argument to set_zone_dissipation()
        double reference_phase_lag
    );

    ///Tell the star to detect its wind saturation state per its current
    ///configuration.
    LIB_PUBLIC void detect_stellar_wind_saturation(
        ///The star to detect the saturation state of.
        EvolvingStar *star
    );

    ///\brief Prepare the zone quantities for interpolation around the
    ///given age.
    ///
    ///See EvolvingStar::select_interpolation_region() for details.
    LIB_PUBLIC void select_interpolation_region(
        ///The star to select the interpolation region for.
        const EvolvingStar *star,

        ///The age around which interpolation should work.
        double age
    );

    ///\brief See Evolve::BrokenPowerlawPhaseLagZone::modified_phase_lag for
    ///details.
    LIB_PUBLIC double modified_phase_lag(
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
    LIB_PUBLIC double core_formation_age(
        ///The star for which to return the core formation age.
        const EvolvingStar *star
    );

    ///The lifetime of a star (the maximum age at which it can be queried)
    LIB_PUBLIC double lifetime(
        ///The star whose lifetime to return.
        const EvolvingStar *star
    );

    ///The luminosity of a star at a given age.
    LIB_PUBLIC double luminosity(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        double age
    );

    ///The luminosity of a star at a given array of ages.
    LIB_PUBLIC void luminosity_array(
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
    LIB_PUBLIC double core_inertia(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        double age
    );

    ///The moment of inertia of the stellar core at a given array of ages.
    LIB_PUBLIC void core_inertia_array(
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
    LIB_PUBLIC double envelope_inertia(
        ///The star whose luminosity to return.
        EvolvingStar *star,

        ///The age at which to return the luminosity.
        double age
    );

    ///\brief The moment of inertia of the stellar envelope at a given array
    ///of ages.
    LIB_PUBLIC void envelope_inertia_array(
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
    LIB_PUBLIC double star_radius(
        ///The star whose radius to return.
        EvolvingStar *star,

        ///The age at which to return the stellar radius.
        double age
    );

    ///The radius of the star at an array of ages.
    LIB_PUBLIC void star_radius_array(
        ///The star whose radius to return.
        EvolvingStar *star,

        ///The ages at which to return the stellar radius.
        const double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );

    ///The radius of the star at a given age.
    LIB_PUBLIC double core_radius(
        ///The star whose radius to return.
        EvolvingStar *star,

        ///The age at which to return the stellar radius.
        double age
    );

    ///The radius of the star at an array of ages.
    LIB_PUBLIC void core_radius_array(
        ///The star whose radius to return.
        EvolvingStar *star,

        ///The ages at which to return the stellar radius.
        const double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );



    ///Converts lg(Q) to a tidal phase lag.
    LIB_PUBLIC double lag_from_lgQ(double lgQ);

} //End extern "C"
