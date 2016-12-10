/**\file
 *
 * \brief Defines some of the methods of the EvolvingStellarQuantity class
 * used for interpolating among stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#include "EvolvingStellarQuantity.h"
#include "mass_metallicity_interp.h"
#include "InterpolatedDerivatives.h"
#include <math.h>
#include <memory>
#include <sstream>

namespace StellarEvolution {

    void EvolvingStellarQuantity::check_grid_range() const
    {
        if(
            __track_masses.size() != 1
            && 
            __track_metallicities.size() != 1
            &&
            (
                __mass < __track_masses[0]
                ||
                __mass > __track_masses[__track_masses.size() - 1]
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
            throw Core::Error::BadFunctionArguments(msg.str());
        }
    }

    void EvolvingStellarQuantity::find_cell(
        const std::valarray<double> &boundaries,
        double value,
        size_t &below_index,
        size_t &above_index
    ) const
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
            metallicity_i < __track_metallicities.size();
            ++metallicity_i
        ) {
            double track_metallicity = __track_metallicities[metallicity_i];
            for(
                size_t mass_i = 0; 
                mass_i < __track_masses.size();
                ++mass_i
            ) {
                double track_mass = __track_masses[mass_i];
                double min_track_age = 
                           __evolution_tracks[track_i]->range_low(),
                       max_track_age = 
                           __evolution_tracks[track_i]->range_high();
                if(__use_log_age) {
                    min_track_age = std::exp(min_track_age);
                    max_track_age = std::exp(max_track_age);
                }
                __min_interp_ages[track_i] = (
                    __initially_zero 
                    ? -Core::Inf 
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

    double EvolvingStellarQuantity::evaluate_track(
        double age,
        const OneArgumentDiffFunction &track,
        const FunctionDerivatives **derivatives
    ) const
    {
        bool too_young = (
            (__use_log_age ? std::log(age) : age) < track.range_low()
        );

        if(derivatives==NULL) return (too_young ? 0.0 : track(age));
        else {
            if(too_young) *derivatives=new Core::ZeroDerivatives;
            else if(__use_log_age) 
                *derivatives=new RemoveLogDeriv(age,
                                                track.deriv(std::log(age)),
                                                true);
            else *derivatives=track.deriv(age);
            return (*derivatives)->order(0);
        }
    }

    void EvolvingStellarQuantity::check_grid_expansion_directions(
        AllowedGridGrowth &grow,
        double age
    ) const
    {
        if(__min_interp_mass_index == 0)
            grow.block_lighter();
        if(__max_interp_mass_index == __track_masses.size())
            grow.block_heavier();
        if(__min_interp_metallicity_index == 0)
            grow.block_poorer();
        if(__max_interp_metallicity_index == __track_metallicities.size())
            grow.block_richer();

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
                !track_in_range(__min_interp_mass_index - 1, 
                                metallicity_index,
                                age)
            )
                grow.block_lighter();

            if(
                grow.heavier()
                &&
                !track_in_range(__max_interp_mass_index,
                                metallicity_index,
                                age)
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
                !track_in_range(mass_index,
                                __min_interp_metallicity_index - 1,
                                age)
            )
                grow.block_poorer();
            if(
                grow.richer()
                &&
                !track_in_range(mass_index,
                                __max_interp_metallicity_index,
                                age)
            )
                grow.block_richer();
        }
    }

    void EvolvingStellarQuantity::expand_grid(const AllowedGridGrowth &grow,
                                              double age) const
    {
        if(
            grow.lighter()
            &&
            track_in_range(__min_interp_mass_index - 1,
                           std::max(size_t(0),
                                    __min_interp_metallicity_index - 1),
                           age)
            &&
            track_in_range(__min_interp_mass_index - 1,
                           std::min(__track_metallicities.size() - 1,
                                    __max_interp_metallicity_index),
                           age)
        ) {
            --__min_interp_mass_index;
            return;
        }

        if(
            grow.heavier()
            &&
            track_in_range(__max_interp_mass_index,
                           std::max(size_t(0),
                                    __min_interp_metallicity_index - 1),
                           age)
            &&
            track_in_range(__max_interp_mass_index,
                           std::min(__track_metallicities.size() - 1,
                                    __max_interp_metallicity_index),
                           age)
        ) {
            ++__max_interp_mass_index;
            return;
        }

        if(
            grow.poorer()
            &&
            track_in_range(std::max(size_t(0), __min_interp_mass_index - 1),
                           __min_interp_metallicity_index - 1,
                           age)
            &&
            track_in_range(std::min(__max_interp_mass_index,
                                    __track_masses.size()),
                           __min_interp_metallicity_index - 1,
                           age)
        ) {
            --__min_interp_metallicity_index;
            return;
        }

        if(
            grow.richer()
            &&
            track_in_range(std::max(size_t(0), __min_interp_mass_index - 1),
                           __max_interp_metallicity_index,
                           age)
            &&
            track_in_range(std::min(__max_interp_mass_index,
                                    __track_masses.size()),
                           __max_interp_metallicity_index,
                           age)
        ) {
            ++__max_interp_metallicity_index;
            return;
        }

        std::valarray<size_t> padding(4);
        padding[0] = __min_interp_mass_index - __mass_index_below;
        padding[1] = (__min_interp_metallicity_index
                      -
                      __metallicity_index_below);
        padding[2] =
            __max_interp_metallicity_index - __metallicity_index_above - 1;
        padding[3] = __max_interp_mass_index - __mass_index_above - 1;

        if(grow.lighter() && padding[0] == padding.min()) {
            --__min_interp_mass_index;
            return;
        }
        padding[0] = std::numeric_limits<size_t>::max();
        if(grow.poorer() && padding[1] == padding.min()) {
            --__min_interp_metallicity_index;
            return;
        }
        padding[1] = std::numeric_limits<size_t>::max();
        if(grow.richer() && padding[2] == padding.min()) {
            ++__max_interp_metallicity_index;
            return;
        }
        ++__max_interp_mass_index;
    }

    void EvolvingStellarQuantity::update_interpolation_grid() const
    {
        std::vector<double>::const_iterator
            lower_age = __next_grid_change_age;
        --lower_age;
        double age = 0.5 * (*lower_age + *__next_grid_change_age);
        if(__min_interp_ages.max() < age && __max_interp_ages.min() > age) {
            __min_interp_mass_index = 0;
            __max_interp_mass_index = __track_masses.size();
            __min_interp_metallicity_index = 0;
            __max_interp_metallicity_index = __track_metallicities.size();
            return;
        }

        AllowedGridGrowth grow;

        if(__mass_index_above == __mass_index_below) 
            grow.block_lighter().block_heavier();
        if(__metallicity_index_above == __metallicity_index_below)
            grow.block_poorer().block_richer();

        __max_interp_mass_index = __mass_index_above;
        __min_interp_mass_index = __mass_index_below;
        __max_interp_metallicity_index = __metallicity_index_above;
        __min_interp_metallicity_index = __metallicity_index_below;
        while(grow) {
            check_grid_expansion_directions(grow, age);
            if(grow) expand_grid(grow, age);
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
        std::vector<double>::const_iterator
            lower_age = __next_grid_change_age;
        --lower_age;
        assert(*lower_age <= age && age <= *__next_grid_change_age);
#endif

        double interp_param = age_to_interp_param(age);
        size_t num_interp_tracks = (
            (
                __max_interp_metallicity_index
                -
                __min_interp_metallicity_index
                +
                1
            )
            *
            (__max_interp_mass_index - __min_interp_mass_index + 1)
        );
        alglib::real_1d_array track_values;
        track_values.setlength(num_interp_tracks);
        std::vector<const FunctionDerivatives *> 
            *track_derivatives 
            = 
            new std::vector<const FunctionDerivatives *>(num_interp_tracks);
        size_t value_index = 0;
        for(
            size_t metallicity_index = __min_interp_metallicity_index;
            metallicity_index <= __max_interp_metallicity_index;
            ++metallicity_index
        ) {
            double
                track_metallicity = __track_metallicities[metallicity_index];
            for(
                size_t mass_index = __min_interp_mass_index;
                mass_index <= __max_interp_mass_index;
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
                    (derivatives
                     ? &((*track_derivatives)[value_index])
                     : NULL)
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
            *derivatives = new InterpolatedDerivatives(__mass,
                                                       __metallicity,
                                                       track_derivatives,
                                                       __interp_masses,
                                                       __interp_metallicities,
                                                       NaN,
                                                       true);
            return (*derivatives)->order(0);
        }
    }

    double EvolvingStellarQuantity::age_to_interp_param(
        double age,
        double mass,
        double metallicity
    ) const
    {
        double metallicity_factor = std::pow(10.0, -metallicity / 5.0);
        return std::log(
            age 
            *
            (1.0 + age / 5.0 * std::pow(mass, 5) * metallicity_factor)
            *
            std::pow(mass, 2.3)
            *
            std::pow(metallicity_factor, 2)
        );
    }

    double EvolvingStellarQuantity::interp_param_to_age(
        double interp_param,
        double mass,
        double metallicity
    ) const
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
        const std::vector<const OneArgumentDiffFunction *> &evolution_tracks,
        bool log_age,
        bool starts_zero
    ) :
        __mass(mass), 
        __metallicity(metallicity),
        __use_log_age(log_age),
        __initially_zero(starts_zero),
        __track_masses(track_masses),
        __track_metallicities(track_metallicities),
        __min_interp_ages(evolution_tracks.size()),
        __max_interp_ages(evolution_tracks.size()),
        __evolution_tracks(evolution_tracks)
    {
        assert(evolution_tracks.size() == (track_masses.size()
                                           *
                                           track_metallicities.size()));
        check_grid_range();
        set_interp_age_ranges();

        if(starts_zero) {
            __interp_grid_change_ages.reserve(evolution_tracks.size() + 1);
            __interp_grid_change_ages.insert(__interp_grid_change_ages.end(),
                                             -Core::Inf);

        } else {
            __interp_grid_change_ages.reserve(2 * evolution_tracks.size());
            __interp_grid_change_ages.insert(
                __interp_grid_change_ages.end(),
                &__min_interp_ages[0],
                &__min_interp_ages[0] + __min_interp_ages.size()
            );
        }
        __interp_grid_change_ages.insert(
            __interp_grid_change_ages.end(),
            &__max_interp_ages[0],
            &__max_interp_ages[0] + __max_interp_ages.size()
        );

        std::sort(__interp_grid_change_ages.begin(),
                  __interp_grid_change_ages.end());
        std::unique(__interp_grid_change_ages.begin(),
                    __interp_grid_change_ages.end());
        __next_grid_change_age = __interp_grid_change_ages.begin();

        find_cell(__track_masses,
                  mass,
                  __mass_index_above,
                  __mass_index_below);
        find_cell(__track_metallicities,
                  metallicity,
                  __metallicity_index_above,
                  __metallicity_index_below);

        __min_age = __min_interp_ages[
            track_index(__mass_index_below, __metallicity_index_below)
        ];
        __max_age = __max_interp_ages[
            track_index(__mass_index_below, __metallicity_index_below)
        ];
        __min_age = std::max(
            __min_age,
            __min_interp_ages[track_index(__mass_index_above,
                                          __metallicity_index_below)]
        );
        __max_age = std::min(
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
        const
    {
        std::vector<double>::const_iterator new_grid_change_age =
            std::upper_bound(
                __interp_grid_change_ages.begin(),
                __interp_grid_change_ages.end(),
                age
            );
        if(__next_grid_change_age != new_grid_change_age) {
            __next_grid_change_age = new_grid_change_age;
            assert(__next_grid_change_age
                   !=
                   __interp_grid_change_ages.end());
            update_interpolation_grid();
        }
    }

    const FunctionDerivatives *EvolvingStellarQuantity::deriv(double age) 
        const
    {
        const FunctionDerivatives *deriv;
        interpolate(age, &deriv);
        return deriv;
    }

    void EvolvingStellarQuantity::enable_next_interpolation_region() const
    {
        assert(__next_grid_change_age != __interp_grid_change_ages.end());
        ++__next_grid_change_age;
        update_interpolation_grid();
    }

} //End StellarEvolution namespace
