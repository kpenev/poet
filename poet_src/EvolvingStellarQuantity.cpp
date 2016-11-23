/**\file
 *
 * \brief Defines some of the methods of the classes for interpolating among
 * stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#include "EvolvingStellarQuantity.h"
#include "Error.h"
#include <math.h>
#include <memory>
#include <sstream>

double mass_metallicity_interp(alglib::real_1d_array interp_masses,
                               alglib::real_1d_array interp_metallicities,
                               alglib::real_1d_array interp_values,
                               double stellar_mass,
                               double stellar_metallicity)
{
    spline2dinterpolant spline;
    spline2dbuildbicubicv(interp_masses, 
                          interp_masses.length(), 
                          interp_metallicities,
                          interp_metallicities.length(),
                          interp_values,
                          1,
                          spline);
    return spline2dcalc(spline, stellar_mass, stellar_metallicity);
}

double LogArgDerivatives::transform_log_arg_deriv(unsigned order) const
{
	if(order==1) return underlying_deriv_values[0]/x;
	else if(order==2)
		return (underlying_deriv_values[1]-underlying_deriv_values[0])/(x*x);
	else throw Error::BadFunctionArguments(
			"Transforming log-derivatives of order higher than 2 is not "
			"implemented."
    );
}

double LogArgDerivatives::order(unsigned deriv_order) const
{
	if(deriv_order==0) {
		if(std::isnan(value)) value=calc_deriv(0);
		return value;
	}
	if(deriv_values.size()<deriv_order)
		for(unsigned i=deriv_values.size(); i<deriv_order; ++i) {
			if(correct_log_arg) {
				underlying_deriv_values.push_back(calc_deriv(i+1));
				deriv_values.push_back(transform_log_arg_deriv(i+1));
			} else deriv_values.push_back(calc_deriv(i+1));
		}
	return deriv_values[deriv_order-1];
}

double InterpolatedDerivatives::calc_deriv(unsigned deriv_order) const
{
	if(deriv_order > 2) return 0.0;
	alglib::real_1d_array interp_values;
    interp_values.setlength(interp_deriv->size());

	std::valarray<double> interp_values(interp_deriv->size());
	for(unsigned i = 0; i < interp_deriv->size(); ++i)
		interp_values[i] = (*interp_deriv)[i]->order(deriv_order);

    return mass_metallicity_interp(__interp_masses,
                                   __interp_metallicities,
                                   interp_values,
                                   stellar_mass,
                                   stellar_metallicity);
}

InterpolatedDerivatives::InterpolatedDerivatives(
    double mass,
    double metallicity,
    std::vector<const FunctionDerivatives*> *derivatives,
    const alglib::real_1d_array &interp_masses,
    const alglib::real_1d_array &interp_metallicities,
    double age,
    bool delete_derivatives,
) :
    LogArgDerivatives(age),
    __stellar_mass(mass),
    __stellar_metallicity(metallicity),
    __interp_deriv(derivatives),
    __interp_masses(interp_masses),
    __interp_metallicities(interp_metallicities),
    __delete_derivatives(delete_derivatives)
{
    assert(derivatives->size() == masses.size() * metallicities.size());
}

double ScaledDerivatives::calc_deriv(unsigned deriv_order=1) const
{
	return (correct_log_arg ? 1.0 : 
			std::pow(scaling, static_cast<int>(deriv_order)))*
		underlying_deriv->order(deriv_order);
}

void EvolvingStellarQuantity::check_grid_range()
{
    if(
        __track_masses.size() != 1
        && 
        __track_metallicities.size() != 1
        &&
        (
            __mass < __track_masses[0]
            ||
            __mass > __track_masses[__track_masess.size() - 1]
            ||
            __metallicity < __track_metallicities[0]
            ||
            __metallicity > __track_metallicities[
                __track_metallicities.size() - 1
            ]
        )
    ) {
        std::ostringstream msg;
		msg << "Stellar mass: "
            << __mass
            << " and metallicity: "
            << __metallicity
			<< " are outside the range of masses: "
            << __track_masses[0] 
            << " - "
            << __track_masses[__track_masses.size() - 1]
            << " and metallicities: "
            << __track_metallicities[0]
            << " - "
            << __track_metallicities[__track_metallicities.size() - 1]
            << " covered by the stellar evolution grid.";
		throw Error::BadFunctionArguments(msg.str());
    }
}

void EvolvingStellarQuantity::find_cell(
    const std::valarray<double> &boundaries
    const double &value,
    size_t &below_index,
    size_t &above_index,
)
{
    const double *first_ptr = &boundaries[0],
                 *above_ptr = std::lower_bound(first_ptr,
                                               first_ptr + boundaries.size(),
                                               value),
                 *below_ptr = above_ptr;
    if( value < *above_ptr) --below_ptr;
    below_index = below_ptr - first_ptr;
    above_index = above_ptr - first_ptr;
}

void EvolvingStellarQuantity::set_interp_age_ranges()
{
    size_t track_i = 0;
    for(
        size_t metallicity_i = 0;
        metallicity_i <= __track_metallicities.size();
        ++metallicity_i
    ) {
        double track_metallicity = __track_metallicities[metallicity_i];
        for(
            size_t mass_i = 0; 
            mass_i <= __track_masses.size();
            ++mass_i
        ) {
            double track_mass = __track_masses[mass_i];
            double min_track_age = evolution_tracks[track_i]->range_low(),
                   max_track_age = evolution_tracks[track_i]->range_high();
            if(__use_log_age) {
                min_track_age = std::exp(min_track_age);
                max_track_age = std::exp(max_track_age);
            }
            __min_interp_ages[track_i] = (
                __initially_zero 
                ? -Inf 
                : interp_param_to_age(
                    age_to_interp_param(min_track_age,
                                        track_mass,
                                        track_metallicity)
                )
            );
            __max_interp_ages[track_i] = interp_param_to_age(
                age_to_interp_param(max_track_age,
                                    track_mass,
                                    track_metallicity)
            );
            ++track_i;
        }
    }
}

double EvolvingStellarQuantity::evaluate_track(double track_age,
		const OneArgumentDiffFunction &track,
		const FunctionDerivatives **derivatives) const
{
	bool too_young = (
        (__use_log_age ? std::log(age) : age) < track.range_low()
    );

    if(derivatives==NULL) return (too_young ? 0.0 : track(age));
    else {
        if(too_young) *derivatives=new ZeroDerivatives;
        else if(use_log_age) 
            *derivatives=new RemoveLogDeriv(age,
                                            track.deriv(std::log(age)),
                                            true);
        else *derivatives=track.deriv(age);
        return (*derivatives)->order(0);
    }
}

void EvolvingStellarQuantity::check_grid_expansion_directions(
    AllowedGridGrowth &grow
)
{
    if(__min_interp_mass_index == 0)
        grow.block_lighter();
    if(__max_interp_mass_index == __track_masses.size())
        grow.block_heavier();
    if(__min_interp_metallicity_index == 0)
        grow.block_poorer();
    if(__max_interp_metallicity_index == __track_metallicities.size())
        grow.block_richer()

    for(
        size_t metallicity_index = __min_interp_metallicity_index;
        (
            (grow.lighter() || grow.heavier())
            &&
            metallicity_index < __max_interp_metallicity_index
        );
        ++metallicity_index
    ) {
        if(
            grow.lighter()
            &&
            !track_in_range(__min_interp_mass_index - 1, metallicity_index)
        )
            grow.block_lighter();

        if(
            grow.heavier()
            &&
            !track_in_range(__max_interp_mass_index, feh_index)
        )
            grow.block_heavier();
    }

    for(
        size_t mass_index = __min_interp_mass_index;
        (
            (grow.poorer() || grow.richer())
            &&
            mass_index < __max_interp_mass_index
        );
        ++mass_index
    ) {
        if(
            grow.poorer()
            &&
            !track_in_range(mass_index, __min_interp_metallicity_index - 1)
        )
            grow.disable_poorer();
        if(
            grow.richer()
            &&
            !track_in_range(mass_index, __max_interp_metallicity_index)
        )
            grow.disable_richer();
    }
}

void EvolvingStellarQuantity::expand_grid(const AllowedGridGrowth &grow)
{
    if(
        grow.lighter()
        &&
        track_in_range(__min_interp_mass_index - 1,
                       std::max(0, __min_interp_metallicity_index - 1))
        &&
        track_in_range(__min_interp_mass_index - 1,
                       std::min(__track_metallicities.size() - 1,
                                __max_interp_metallicity_index))
    ) {
        --__min_interp_mass_index;
        return;
    }

    if(
        grow.heavier()
        &&
        track_in_range(__max_interp_mass_index,
                       std::max(0, __min_interp_metallicity_index - 1))
        &&
        track_in_range(__max_interp_mass_index,
                       std::min(__track_metallicities.size() - 1,
                                __max_interp_metallicity_index))
    ) {
        ++__max_interp_mass_index;
        return;
    }

    if(
        grow.poorer()
        &&
        track_in_range(std::max(0, __min_interp_mass_index - 1),
                       __min_interp_metallicity_index - 1)
        &&
        track_in_range(std::min(__max_interp_mass_index,
                                __track_masses.size()),
                       __min_interp_metallicity_index - 1)
    ) {
        --__min_interp_metallicity_index;
        return;
    }

    if(
        grow.richer()
        &&
        track_in_range(std::max(0, __min_interp_mass_index - 1),
                       __max_interp_metallicity_index)
        &&
        track_in_range(std::min(__max_interp_mass_index,
                                __track_masses.size()),
                       __max_interp_metallicity_index)
    ) {
        ++__max_interp_metallicity_index;
        return;
    }

    std::vallarray<size_t> padding(4);
    padding[0] = __min_interp_mass_index - __mass_index_below;
    padding[1] = __min_interp_metallicity_index - __metallicity_index_below;
    padding[2] =
        __max_interp_metallicity_index - __metallicity_index_above - 1
    padding[3] = __max_interp_mass_index - __mass_index_above - 1,

    if(grow.lighter() && padding[0] == padding.min()) {
        --__min_interp_mass_index;
        return;
    }
    padding[0] = std::numeric_limits<size_t>::max();
    if(grow_poorer() && padding[1] == padding.min()) {
        --__min_interp_metallicity_index;
        return;
    }
    padding[1] = std::numeric_limits<size_t>::max();
    if(grow_richer() && padding[2] == padding.min()) {
        ++__max_interp_metallicity_index;
        return;
    }
    ++__max_interp_mass_index;
}

void EvolvingStellarQuantity::update_interpolation_grid()
{
    std::vector<double>::const_iterator lower_age = __next_grid_change_age;
    --lower_age;
    double age = 0.5 * (*lower_age + *__next_grid_change_age);
    if(min_interp_ages.max() < age && max_interp_ages.min() > age) {
        __min_interp_mass_index = 0;
        __max_interp_mass_index = __track_masses.size();
        __min_interp_metallicity_index = 0;
        __max_interp_metallicity_index = __track_metallicities.size();
        return;
    }

    min_mass_index = __mass_index_below;
    max_mass_index = __mass_index_above + 1;
    min_feh_index = __metallicity_index_below;
    max_feh_index = __metallicity_index_above + 1;

    AllowedGridGrowth grow;

    if(__mass_index_above == __mass_index_below) 
        grow.block_lighter().block_heavier();
    if(__metallicity_index_above == __metallicity_index_below)
        grow.block_poorer().block_richer();

    while(grow) {
        check_grid_expansion_directions(grow);
        if(grow) expand_grid(grow);
    }
    __interp_masses.setcontent(
        __max_interp_mass_index - __min_interp_mass_index,
        &__track_masses[__min_interp_mass_index]
    );
    __interp_metallicities.setcontent(
        __max_interp_metallicity_index - __min_interp_metallicity_index,
        &__track_metallicities[__min_interp_metallicity_index]
    );
}

double EvolvingStellarQuantity::interpolate(
    double age,
    const FunctionDerivatives **derivatives
) const
{
    assert(__min_age <= age && age <= __max_age);
    assert(__next_grid_change_age != __interp_grid_change_ages.begin());
#ifndef NDEBUG
    std::vector<double>::const_iterator lower_age = __next_grid_change_age;
    --lower_age;
    assert(*lower_age <= age && age <= *__next_grid_change_age);
#endif

    double interp_param = age_to_interp_param(age);
    size_t num_interp_tracks = (
        (__max_interp_metallicity_index - __min_interp_metallicity_index)
        *
        (__max_interp_mass_index - __min_interp_mass_index)
    )
    std::valarray<double> track_values(num_interp_tracks);
    std::vector<const FunctionDerivatives *> 
        *track_derivatives 
        = 
        new std::vector<const FunctionDerivatives *>(num_interp_tracks);
    size_t value_index = 0;
    for(
        size_t metallicity_index = __min_interp_metallicity_index;
        metallicity_index < __max_interp_metallicity_index;
        ++metallicity_index
    ) {
        double track_metallicity = __track_metallicities[metallicity_index];
        for(
            size_t mass_index = __min_interp_mass_index;
            mass_index < __max_interp_mass_index;
            ++mass_index
        ) {
            double track_mass = __track_masses[mass_index];
            double track_age = interp_param_to_age(interp_param,
                                                   track_mass,
                                                   track_metallicity);
            track_values[value_index] = evaluate_track(
                track_age,
                *__evolution_tracks[track_index(mass_index,
                                                metallicity_index)],
                (derivatives ? &((*track_derivatives)[good_ind]) : NULL)
            );
            ++value_index;
        }
    }
    if(derivatives == NULL)
        return mass_metallicity_interp(__interp_masses,
                                       __interp_metallicities,
                                       track_values,
                                       __mass,
                                       __metallicity);
    else {
		*derivatives=new InterpolatedDerivatives(__mass,
                                                 __metallicity
                                                 track_derivatives,
                                                 __interp_masses,
                                                 __interp_metallicities,
                                                 NaN,
                                                 true);
        return (*derivatives)->order(0);
    }
}

double EvolvingStellarQuantity::age_to_interp_param(double age,
                                                    double mass,
                                                    double metallicity)
{
    double metallicity_factor = std::pow(10.0, -metallicity / 5.0)
    return std::log(
        age 
        *
        (
            1.0
            +
            age / 5.0 * std::pow(mass, 5) * metallicity_factor
        )
        *
        std::pow(mass, 2.3)
        *
        std::pow(metallicity_factor, 2)
    );
}

double EvolvingStellarQuantity::interp_param_to_age(double interp_param,
                                                    double mass,
                                                    double metallicity)
{
    double metallicity_factor = std::pow(10.0, -metallicity / 5.0),
           c = (std::exp(interp_param)
                *
                std::pow(mass, 2.3)
                *
                std::pow(metallicity_factor, 2)),
           a = 0.2 * std::pow(mass, 5) * metallicity_factor,
           discr = 1.0 + 4.0 * a * c;
    return (std::sqrt(discr) - 1.0) / (2.0 * a);
}

EvolvingStellarQuantity::EvolvingStellarQuantity(
    double mass,
    double metallicity,
    const std::valarray<double> &track_masses,
    const std::valarray<double> &track_metallicities,
    const std::list<const OneArgumentDiffFunction *> &evolution_tracks,
    bool log_age,
    double extrapolate,
    bool starts_zero
) :
    __mass(mass), 
    __metallicity(metallicity),
    __use_log_age(log_age),
    __initially_zero(starts_zero),
    __track_masses(masses_of_tracks),
    __track_metallicities(metallicities_of_tracks),
    __min_interp_ages(evolution_tracks.size()),
    __max_interp_ages(evolution_tracks.size()),
    __evolution_tracks(evolution_tracks.begin(), evolution_tracks.end()),
	__age_scaling_mass(low_mass_age_scaling),
	__extrapolate(extrapolate),
{
    assert(evolution_tracks.size() == (masses_of_tracks.size()
                                       *
                                       metallicities_of_tracks.size()));
    check_grid_range();
    set_interp_age_ranges();

    if(starts_zero) {
        __interp_grid_change_ages.reserve(evolution_tracks.size() + 1);
        __interp_grid_change_ages.insert(__interp_grid_change_ages.end(),
                                         -Inf);

    } else {
        __interp_grid_change_ages.reserve(2 * evolution_tracks.size());
        __interp_grid_change_ages.insert(
            __interp_grid_change_ages.end(),
            &__min_interp_ages[0],
            &__min_interp_ages[0] + __min_interp_ages.size()]
        );
    }
    __interp_grid_change_ages.insert(
        __interp_grid_change_ages.end(),
        &__max_interp_ages[0],
        &__max_interp_ages[0] + __max_interp_ages.size()]
    );

    std::sort(__interp_grid_change_ages.begin(),
              __interp_grid_change_ages.end());
    std::unique(__interp_grid_change_ages.begin(),
                __interp_grid_change_ages.end());
    __next_grid_change_age = __interp_grid_change_ages.begin();

    find_cell(__track_masses, mass, __mass_index_above, __mass_index_below);
    find_cell(__track_metallicities,
              metallicity,
              __metallicity_index_above,
              __metallicity_index_below);

    __min_age = __min_interp_ages[track_index(__mass_index_below,
                                              __metallicity_index_below)];
    __max_age = __max_interp_ages[track_index(__mass_index_below,
                                              __metallicity_index_below)];
    __min_age = std::max(
        __min_age,
        __min_interp_ages[track_index(__mass_index_above,
                                      __metallicity_index_below)]
    );
    __min_age = std::min(
        __max_age,
        __max_interp_ages[track_index(__mass_index_above,
                                      __metallicity_index_below)]
    );

    __min_age = std::max(
        __min_age,
        __min_interp_ages[track_index(__mass_index_below,
                                      __metallicity_index_above)]
    );
    __max_age = std::min(
        __max_age,
        __max_interp_ages[track_index(__mass_index_below,
                                      __metallicity_index_above)]
    );

    __min_age = std::max(
        __min_age,
        __min_interp_ages[track_index(__mass_index_above,
                                      __metallicity_index_above)]
    );
    __max_age = std::min(
        __max_age,
        __max_interp_ages[track_index(__mass_index_above,
                                      __metallicity_index_above)]
    );
}

void EvolvingStellarQuantity::select_interpolation_region(double age)
{
    __next_grid_change_age = std::upper_bound(
        __interp_grid_change_ages.begin(),
        __interp_grid_change_ages.end()
    );
    update_interpolation_grid();
}

const FunctionDerivatives *EvolvingStellarQuantity::deriv(double age) const
{
    const FunctionDerivatives *deriv;
    interpolate(age, &deriv);
    return deriv;
}

double EvolvingStellarQuantity::enable_next_interpolation_region()
{
    assert(__next_grid_change_age != __interp_grid_change_ages.end());
    ++__next_grid_change_age;
    update_interpolation_grid();
}
