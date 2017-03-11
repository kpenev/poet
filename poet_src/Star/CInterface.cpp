/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h.
 *
 * \ingroup Star_group
 */

#include "CInterface.h"

const int NO_DERIV = Evolve::Dissipation::NO_DERIV;
const int AGE_DERIV = Evolve::Dissipation::AGE;
const int SPIN_FREQUENCY_DERIV = Evolve::Dissipation::SPIN_FREQUENCY;
const int ORBITAL_FREQUENCY_DERIV = Evolve::Dissipation::ORBITAL_FREQUENCY;

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
    if(zone_index == 0) zone = &(real_star->envelope());
    else zone = &(real_star->core());
    
    zone->setup(
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
        reference_phase_lag
    );
}

void detect_stellar_wind_saturation(EvolvingStar *star)
{
    reinterpret_cast<Star::InterpolatedEvolutionStar*>(
        star
    )->detect_saturation();
}

void select_interpolation_region(const EvolvingStar *star,
                                 double age)
{
    return reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
        star
    )->select_interpolation_region(age);
}

double modified_phase_lag(const EvolvingStar *star,
                          unsigned zone_index,
                          int orbital_frequency_multiplier,
                          int spin_frequency_multiplier,
                          double forcing_frequency,
                          int deriv,
                          double *above_lock_value)
{
    return reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
        star
    )->zone(zone_index).modified_phase_lag(
        orbital_frequency_multiplier,
        spin_frequency_multiplier,
        forcing_frequency,
        static_cast<Evolve::Dissipation::Derivative>(deriv),
        *above_lock_value
    );
}

double core_formation_age(const EvolvingStar *star)
{
    return reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
        star
    )->core().formation_age();
}

double lifetime(const EvolvingStar *star)
{
    return reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
        star
    )->lifetime();
}

double luminosity(EvolvingStar *star, double age)
{
    const Star::InterpolatedEvolutionStar *actual_star = 
        reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
            star
        );
    actual_star->select_interpolation_region(age);
    return actual_star->luminosity(age);
}

void luminosity_array(EvolvingStar *star,
                      const double *age,
                      unsigned nvalues,
                      double *result)
{
    Star::InterpolatedEvolutionStar *actual_star = 
        reinterpret_cast<Star::InterpolatedEvolutionStar*>(
            star
        );
    for(unsigned i = 0; i < nvalues; ++i) {
        if(
            i == 0
            ||
            age[i] < age[i - 1]
        )
            actual_star->select_interpolation_region(age[i]);
        else 
            actual_star->reached_critical_age(age[i]);
        result[i] = actual_star->luminosity(age[i]);
    }
}

double core_inertia(EvolvingStar *star, double age)
{
    const Star::EvolvingStellarCore &core =
        reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
            star
        )->core();
    core.select_interpolation_region(age);
    return core.moment_of_inertia(age);
}

void zone_inertia_array(Star::EvolvingStellarZone &zone,
                        const double *age,
                        unsigned nvalues,
                        double *result)
{
    for(unsigned i = 0; i < nvalues; ++i) {
        if(
            i == 0
            ||
            age[i] < age[i - 1]
        )
            zone.select_interpolation_region(age[i]);
        else
            zone.reached_critical_age(age[i]);
        result[i] = zone.moment_of_inertia(age[i]);
    }
}

void core_inertia_array(EvolvingStar *star,
                        const double *age,
                        unsigned nvalues,
                        double *result)
{
    zone_inertia_array(
        reinterpret_cast<Star::InterpolatedEvolutionStar*>(
            star
        )->core(),
        age,
        nvalues,
        result
    );
}

double envelope_inertia(EvolvingStar *star, double age)
{
    const Star::EvolvingStellarEnvelope &envelope =
        reinterpret_cast<const Star::InterpolatedEvolutionStar*>(
            star
        )->envelope();
    envelope.select_interpolation_region(age);
    return envelope.moment_of_inertia(age);
}

void envelope_inertia_array(EvolvingStar *star,
                            const double *age,
                            unsigned nvalues,
                            double *result)
{
    zone_inertia_array(
        reinterpret_cast<Star::InterpolatedEvolutionStar*>(
            star
        )->envelope(),
        age,
        nvalues,
        result
    );
}