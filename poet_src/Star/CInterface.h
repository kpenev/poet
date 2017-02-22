/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * Star library.
 *
 * \ingroup Star_group
 */

#include "EvolvingStar.h"
#include "../StellarEvolution/CInterface.h"

extern "C" {
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
} //End extern "C"
