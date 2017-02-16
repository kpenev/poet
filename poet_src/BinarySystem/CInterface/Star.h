#include "../EvolvingStar.h"

extern "C" {
    ///Create a star to participate in the orbital evolution calculation.
    InterpolatedEvolutionStar *create_star(
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
        const StellarEvolution::Interpolator *interpolator
    );

    ///Destroy a previously created star.
    void destroy_star(
        ///The star to destroy. Must have previously been created using
        ///create_star.
        InterpolatedEvolutionStar *star
    );
} //End extern "C"
