/**\file
 * 
 * \brief Declares functions for creating stars used by unit tests.
 *
 * \ingroup UnitTests_group
 */

#ifndef __MAKE_STAR_H
#define __MAKE_STAR_H

#include "../../Star/EvolvingStar.h"

///Create a star with the given parameters with a constant phase lag.
Star::InterpolatedEvolutionStar *make_const_lag_star(
    const StellarEvolution::Interpolator &evolution,
    double wind_strength,
    double wind_sat_freq,
    double coupling_timescale,
    double phase_lag = 0
);

#endif
