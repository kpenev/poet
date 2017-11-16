/**\file
 *
 * \brief Defines some of the methods of the EvolvingStellarQuantity class
 * used for interpolating among stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#define BUILDING_LIBRARY
#include "EvolvingStellarQuantity.h"
#include <memory>

namespace StellarEvolution {

    void EvolvingStellarQuantity::check_grid_range() const
    {
        if(
            __track_masses.size() != 1
            && 
            __track_feh.size() != 1
            &&
            (
                __mass < __track_masses[0]
                ||
                __mass > __track_masses[__track_masses.size() - 1]
                ||
                __feh < __track_feh[0]
                ||
                __feh > __track_feh[__track_feh.size() - 1]
            )
        ) {
            std::ostringstream msg;
            msg << "Stellar mass: "
                << __mass
                << " and [Fe/H]: "
                << __feh
                << " are outside the range of masses: "
                << __track_masses[0] 
                << " - "
                << __track_masses[__track_masses.size() - 1]
                << " and [Fe/H]: "
                << __track_feh[0]
                << " - "
                << __track_feh[__track_feh.size() - 1]
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
                     *last_ptr = first_ptr + boundaries.size(),
                     *above_ptr = std::lower_bound(first_ptr,
                                                   last_ptr,
                                                   value),
                     *below_ptr = above_ptr;
#ifndef NDEBUG
        if(above_ptr >= last_ptr )
            std::cerr << "Value: "
                      << value
                      << " out of bounds: ["
                      << boundaries[0]
                      << ", "
                      << boundaries[boundaries.size() - 1]
                      << "]"
                      << std::endl;
#endif

        assert( above_ptr < last_ptr );

        if( value < *below_ptr) --below_ptr;
#ifndef NDEBUG
        if(*above_ptr < value)
            std::cerr << "Value: " << value
                      << " not in [" << *below_ptr << ", " << *above_ptr << "]"
                      << std::endl;
#endif
        assert(*above_ptr >= value);
        assert(value >= *below_ptr);
        below_index = below_ptr - first_ptr;
        above_index = above_ptr - first_ptr;
    }

    void EvolvingStellarQuantity::set_interp_age_ranges()
    {
        size_t track_i = 0;
        for(
            size_t feh_i = 0;
            feh_i < __track_feh.size();
            ++feh_i
        ) {
            double track_feh = __track_feh[feh_i];
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
                if(__log_age) {
                    min_track_age = std::exp(min_track_age);
                    max_track_age = std::exp(max_track_age);
                }
                __min_interp_ages[track_i] = interp_param_to_age(
                    age_to_interp_param(min_track_age,
                                        track_mass,
                                        track_feh)
                );
                __max_interp_ages[track_i] = interp_param_to_age(
                    age_to_interp_param(max_track_age,
                                        track_mass,
                                        track_feh)
                );
#ifndef NDEBUG
                std::cerr << "Track "
                          << track_i
                          << " age range: ["
                          << __min_interp_ages[track_i]
                          << ", "
                          << __max_interp_ages[track_i]
                          << "]"
                          << std::endl;
#endif
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
        double track_argument = (__log_age ? std::log(age) : age);
        if(track_argument < track.range_low()) {
            assert(std::abs(track_argument - track.range_low()) 
                   <
                   1e-8 * std::abs(track_argument)); 
            track_argument = track.range_low();
        }

        if(derivatives) {
            *derivatives = new RemoveLogDeriv(
                (__log_age ? age : NaN),
                false,
                track.deriv(track_argument),
                true
            );
            return (*derivatives)->order(0);
        } else return track(track_argument);
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
        if(__min_interp_feh_index == 0)
            grow.block_poorer();
        if(__max_interp_feh_index == __track_feh.size())
            grow.block_richer();

        for(
            size_t feh_index = __min_interp_feh_index;
            (
                (grow.lighter() || grow.heavier())
                &&
                feh_index < __max_interp_feh_index
            );
            ++feh_index
        ) {
            if(
                grow.lighter()
                &&
                !track_in_range(__min_interp_mass_index - 1,
                                feh_index,
                                age)
            )
                grow.block_lighter();

            if(
                grow.heavier()
                &&
                !track_in_range(__max_interp_mass_index,
                                feh_index,
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
                                __min_interp_feh_index - 1,
                                age)
            )
                grow.block_poorer();
            if(
                grow.richer()
                &&
                !track_in_range(mass_index,
                                __max_interp_feh_index,
                                age)
            )
                grow.block_richer();
        }
    }

    bool EvolvingStellarQuantity::expand_grid(const AllowedGridGrowth &grow,
                                              double age) const
    {
        if(
            grow.lighter()
            &&
            track_in_range(__min_interp_mass_index - 1,
                           std::max(size_t(1),
                                    __min_interp_feh_index) - 1,
                           age)
            &&
            track_in_range(__min_interp_mass_index - 1,
                           std::min(__track_feh.size() - 1,
                                    __max_interp_feh_index),
                           age)
        ) {
            --__min_interp_mass_index;
            return true;
        }

        if(
            grow.heavier()
            &&
            track_in_range(__max_interp_mass_index,
                           std::max(size_t(1),
                                    __min_interp_feh_index) - 1,
                           age)
            &&
            track_in_range(__max_interp_mass_index,
                           std::min(__track_feh.size() - 1,
                                    __max_interp_feh_index),
                           age)
        ) {
            ++__max_interp_mass_index;
            return true;
        }

        if(
            grow.poorer()
            &&
            track_in_range(std::max(size_t(1), __min_interp_mass_index) - 1,
                           __min_interp_feh_index - 1,
                           age)
            &&
            track_in_range(std::min(__max_interp_mass_index,
                                    __track_masses.size() - 1),
                           __min_interp_feh_index - 1,
                           age)
        ) {
            --__min_interp_feh_index;
            return true;
        }

        if(
            grow.richer()
            &&
            track_in_range(std::max(size_t(1), __min_interp_mass_index) - 1,
                           __max_interp_feh_index,
                           age)
            &&
            track_in_range(std::min(__max_interp_mass_index,
                                    __track_masses.size() - 1),
                           __max_interp_feh_index,
                           age)
        ) {
            ++__max_interp_feh_index;
            return true;
        }

        std::valarray<size_t> padding(4);
        padding[0] = __min_interp_mass_index - __mass_index_below;
        padding[1] = (__min_interp_feh_index
                      -
                      __feh_index_below);
        padding[2] =
            __max_interp_feh_index - __feh_index_above - 1;
        padding[3] = __max_interp_mass_index - __mass_index_above - 1;

        if(grow.lighter() && padding[0] == padding.min()) {
            --__min_interp_mass_index;
            return true;
        }
        padding[0] = std::numeric_limits<size_t>::max();
        if(grow.poorer() && padding[1] == padding.min()) {
            --__min_interp_feh_index;
            return true;
        }
        padding[1] = std::numeric_limits<size_t>::max();
        if(grow.richer() && padding[2] == padding.min()) {
            ++__max_interp_feh_index;
            return true;
        }
        if(grow.heavier()) {
            ++__max_interp_mass_index;
            return true;
        }
        return false;
    }

    void EvolvingStellarQuantity::update_interpolation_grid() const
    {
        std::vector<double>::const_iterator
            lower_age = __next_grid_change_age;
        if(__next_grid_change_age == __interp_grid_change_ages.begin()) {
            assert(__initially_zero);
            __min_interp_mass_index = 0;
            __max_interp_mass_index = 0;
            __min_interp_feh_index = 0;
            __max_interp_feh_index = 0;
            __interp_masses.setlength(0);
            __interp_feh.setlength(0);
            return;
        }
        --lower_age;
        double age = 0.5 * (*lower_age + *__next_grid_change_age);

        if(__min_interp_ages.max() < age && __max_interp_ages.min() > age) {
            __min_interp_mass_index = 0;
            __max_interp_mass_index = __track_masses.size();
            __min_interp_feh_index = 0;
            __max_interp_feh_index = __track_feh.size();
        } else {

            AllowedGridGrowth grow;

            if(__mass_index_above == __mass_index_below) 
                grow.block_lighter().block_heavier();
            if(__feh_index_above == __feh_index_below)
                grow.block_poorer().block_richer();

            __max_interp_mass_index = __mass_index_above + 1;
            __min_interp_mass_index = __mass_index_below;
            __max_interp_feh_index = __feh_index_above + 1;
            __min_interp_feh_index = __feh_index_below;
            while(grow) {
                check_grid_expansion_directions(grow, age);
                if(grow && !expand_grid(grow, age)) break;
            }
        }

        __interp_masses.setcontent(
            __max_interp_mass_index - __min_interp_mass_index,
            &__track_masses[__min_interp_mass_index]
        );
        __interp_feh.setcontent(
            __max_interp_feh_index - __min_interp_feh_index,
            &__track_feh[__min_interp_feh_index]
        );
    }

    double EvolvingStellarQuantity::interpolate(
        double age,
        const FunctionDerivatives **derivatives
    ) const
    {
        if(age < __min_age && __initially_zero) {
            if(derivatives) *derivatives=new Core::ZeroDerivatives;
            return 0.0;
        }
#ifndef NDEBUG
        if(age < __min_age || age > __max_age)
            std::cerr << "Age: " << age
                      << " not in [" << __min_age << ", " << __max_age << "]"
                      << std::endl;
#endif
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
            (__max_interp_feh_index - __min_interp_feh_index)
            *
            (__max_interp_mass_index - __min_interp_mass_index)
        );
        alglib::real_1d_array track_values;
        track_values.setlength(num_interp_tracks);
        std::vector<const FunctionDerivatives *> 
            *track_derivatives 
            = 
            new std::vector<const FunctionDerivatives *>(num_interp_tracks);
        size_t value_index = 0;
        for(
            size_t feh_index = __min_interp_feh_index;
            feh_index < __max_interp_feh_index;
            ++feh_index
        ) {
            double
                track_feh = __track_feh[feh_index];
            for(
                size_t mass_index = __min_interp_mass_index;
                mass_index < __max_interp_mass_index;
                ++mass_index
            ) {
                double track_mass = __track_masses[mass_index];
                double track_age = interp_param_to_age(interp_param,
                                                       track_mass,
                                                       track_feh);
                track_values[value_index] = evaluate_track(
                    track_age,
                    *__evolution_tracks[track_index(mass_index,
                                                    feh_index)],
                    (derivatives
                     ? &((*track_derivatives)[value_index])
                     : NULL)
                );
                ++value_index;
            }
        }
        if(derivatives == NULL) {
            double result = mass_feh_interp(__interp_masses,
                                            __interp_feh,
                                            track_values,
                                            __mass,
                                            __feh);
            delete track_derivatives;
            return (__log_quantity ? std::exp(result) : result);
        } else {
            *derivatives = new InterpolatedDerivatives(
                __mass,
                __feh,
                track_derivatives,
                __interp_masses,
                __interp_feh,
                NaN,
                __log_quantity,
                true
            );
            return (*derivatives)->order(0);
        }
    }

    double EvolvingStellarQuantity::age_to_interp_param(
        double age,
        double mass,
        double feh
    ) const
    {
        //t * (1.0 + (t / 5.0) * m**5 * 10.0**(-0.2*feh))* m**2.3 * 10.0**(-0.4*feh)
        double feh_factor = std::pow(10.0, -feh / 5.0);
        return std::log(
            age 
            *
            (1.0 + age / 5.0 * std::pow(mass, 5) * feh_factor)
            *
            std::pow(mass, 2.3)
            *
            std::pow(feh_factor, 2)
        );
    }

    double EvolvingStellarQuantity::interp_param_to_age(
        double interp_param,
        double mass,
        double feh
    ) const
    {
        double feh_factor = std::pow(10.0, -feh / 5.0),
               c = (std::exp(interp_param)
                    *
                    std::pow(mass, -2.3)
                    /
                    std::pow(feh_factor, 2)),
               a = 0.2 * std::pow(mass, 5) * feh_factor,
               discr = 1.0 + 4.0 * a * c;
        return (std::sqrt(discr) - 1.0) / (2.0 * a);
    }

    EvolvingStellarQuantity::EvolvingStellarQuantity(
        double mass,
        double feh,
        const std::valarray<double> &track_masses,
        const std::valarray<double> &track_feh,
        const std::vector<const OneArgumentDiffFunction *> &evolution_tracks,
        bool log_age,
        bool log_quantity,
        bool starts_zero
    ) :
        __mass(mass), 
        __feh(feh),
        __log_age(log_age),
        __log_quantity(log_quantity),
        __initially_zero(starts_zero),
        __track_masses(track_masses),
        __track_feh(track_feh),
        __min_interp_ages(evolution_tracks.size()),
        __max_interp_ages(evolution_tracks.size()),
        __evolution_tracks(evolution_tracks)
    {
#ifndef NDEBUG
        std::cerr << "Creating quantity with mass: " << mass
                  << ", [Fe/H]: " << feh << std::endl
                  << ", with track masses: " << track_masses << std::endl
                  << " and track [Fe/H]: " << track_feh << std::endl
                  << " with " << evolution_tracks.size() << " tracks"
                  << std::endl;
#endif

        assert(evolution_tracks.size() == (track_masses.size()
                                           *
                                           track_feh.size()));
        check_grid_range();
        set_interp_age_ranges();

        if(starts_zero) {
            __interp_grid_change_ages.reserve(2 * evolution_tracks.size()
                                              +
                                              1);
            __interp_grid_change_ages.insert(__interp_grid_change_ages.end(),
                                             -Core::Inf);

        } else {
            __interp_grid_change_ages.reserve(2 * evolution_tracks.size());
        }
        __interp_grid_change_ages.insert(
            __interp_grid_change_ages.end(),
            &__min_interp_ages[0],
            &__min_interp_ages[0] + __min_interp_ages.size()
        );
        __interp_grid_change_ages.insert(
            __interp_grid_change_ages.end(),
            &__max_interp_ages[0],
            &__max_interp_ages[0] + __max_interp_ages.size()
        );

        std::sort(__interp_grid_change_ages.begin(),
                  __interp_grid_change_ages.end());
        std::unique(__interp_grid_change_ages.begin(),
                    __interp_grid_change_ages.end());
        __next_grid_change_age = __interp_grid_change_ages.end();

        find_cell(__track_masses,
                  mass,
                  __mass_index_below,
                  __mass_index_above);
        find_cell(__track_feh,
                  feh,
                  __feh_index_below,
                  __feh_index_above);

        __min_age = __min_interp_ages[
            track_index(__mass_index_below, __feh_index_below)
        ];
        __max_age = __max_interp_ages[
            track_index(__mass_index_below, __feh_index_below)
        ];
        __min_age = std::max(
            __min_age,
            __min_interp_ages[track_index(__mass_index_above,
                                          __feh_index_below)]
        );
        __max_age = std::min(
            __max_age,
            __max_interp_ages[track_index(__mass_index_above,
                                          __feh_index_below)]
        );

        __min_age = std::max(
            __min_age,
            __min_interp_ages[track_index(__mass_index_below,
                                          __feh_index_above)]
        );
        __max_age = std::min(
            __max_age,
            __max_interp_ages[track_index(__mass_index_below,
                                          __feh_index_above)]
        );

        __min_age = std::max(
            __min_age,
            __min_interp_ages[track_index(__mass_index_above,
                                          __feh_index_above)]
        );
        __max_age = std::min(
            __max_age,
            __max_interp_ages[track_index(__mass_index_above,
                                          __feh_index_above)]
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

    double EvolvingStellarQuantity::previous_discontinuity() const
    {
        if(__next_grid_change_age == __interp_grid_change_ages.begin()) {
            assert(__initially_zero);
            return -Core::Inf;
        }
        std::vector<double>::const_iterator result = __next_grid_change_age;
        return *(--result);
    }

    void EvolvingStellarQuantity::enable_next_interpolation_region() const
    {
        assert(__next_grid_change_age != __interp_grid_change_ages.end());
        ++__next_grid_change_age;
        update_interpolation_grid();
    }

} //End StellarEvolution namespace
