/**\file
 * 
 * \brief Implement the functions defined in MakeStar.h.
 *
 * \ingroup UnitTests_group
 */

#include "MakeStar.h"

Star::InterpolatedEvolutionStar *make_const_lag_star(
    const StellarEvolution::Interpolator &evolution,
    double wind_strength,
    double wind_sat_freq,
    double coupling_timescale,
    double phase_lag
)
{
    Star::InterpolatedEvolutionStar *star = (
        new Star::InterpolatedEvolutionStar(1.0,//mass
                                            0.0,//feh
                                            wind_strength,
                                            wind_sat_freq,
                                            coupling_timescale,
                                            evolution)
    );
    Evolve::BrokenPowerlawPhaseLagZone *zone = &(star->envelope());
    for(int i = 0; i < 2; ++i) {
        zone->setup(std::vector<double>(),//Wtide breaks
                    std::vector<double>(),//W* breaks
                    std::vector<double>(1, 0.0),//Wtide pow.
                    std::vector<double>(1, 0.0),//W* pow.
                    (i==0 ? phase_lag : 0.0));
        zone = &(star->core());
    }
    return star;
}
