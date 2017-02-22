/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h.
 *
 * \ingroup Star_group
 */

#include "CInterface.h"

///Converts lg(Q) to a tidal phase lag.
double lag_from_lgQ(double lgQ)
{
	return 15.0 / (16.0 * M_PI * std::pow(10.0, lgQ));
}

EvolvingStar *create_star(double mass,
                          double metallicity,
                          double wind_strength,
                          double wind_saturation_frequency,
                          double diff_rot_coupling_timescale,
                          const MESAInterpolator *interpolator)
{
    return reinterpret_cast<EvolvingStar*>(
        new Star::InterpolatedEvolutionStar(
            mass,
            metallicity,
            wind_strength,
            wind_saturation_frequency,
            diff_rot_coupling_timescale,
            *reinterpret_cast<const StellarEvolution::MESA::Interpolator*>(
                interpolator
            )
        )
    );
}

///Destroy a previously created star.
void destroy_star(EvolvingStar *star)
{
    delete reinterpret_cast<Star::InterpolatedEvolutionStar*>(star);
}

void set_dissipation(EvolvingStar *star,
                     unsigned zone_index,
                     unsigned num_tidal_frequency_breaks,
                     unsigned num_spin_frequency_breaks,
                     double *tidal_frequency_breaks,
                     double *spin_frequency_breaks,
                     double *tidal_frequency_powers,
                     double *spin_frequency_powers,
                     double reference_phase_lag)
{
    Star::InterpolatedEvolutionStar* real_star = 
        reinterpret_cast<Star::InterpolatedEvolutionStar*>(star);
    Evolve::BrokenPowerlawPhaseLagZone *zone;
    if(zone_index == 0) zone = &real_star->envelope();
    else zone = &real_star->core();
    
    zone->setup(
        std::vector<double>(
            tidal_frequency_breaks, 
            tidal_frequency_breaks + num_tidal_frequency_breaks
        ),
        std::vector<double>(
            spin_frequency_breaks,
            spin_frequency_breaks + num_spin_frequency_breaks
        ),
        std::vector<double>(
            tidal_frequency_powers,
            tidal_frequency_powers + num_tidal_frequency_breaks + 1
        ),
        std::vector<double>(
            spin_frequency_powers,
            spin_frequency_powers + num_spin_frequency_breaks + 1
        ),
        reference_phase_lag
    );
}

void detect_stellar_wind_saturation(EvolvingStar *star)
{
    reinterpret_cast<Star::InterpolatedEvolutionStar*>(
        star
    )->detect_saturation();
}
