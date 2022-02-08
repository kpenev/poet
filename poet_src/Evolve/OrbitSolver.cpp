/**\file
 *
 * \brief Implements some of the members of the OrbitSolver class, the
 * various stopping conditions and a number of other classes used while
 * calculating the orbital evolution.
 *
 * \ingroup OrbitSolver_group
 */

#define BUILDING_LIBRARY
#include "OrbitSolver.h"
#include <iostream>
#include <iomanip>

namespace Evolve {

    const double MIN_RELATIVE_STEP = (
        1.0
        +
        10.0 * std::numeric_limits<double>::epsilon()
    );

    int stellar_system_diff_eq(double age,
                               const double *parameters,
                               double *derivatives,
                               void *system_mode)
    {
        void **input_params = static_cast<void **>(system_mode);
        BinarySystem &system = *static_cast< BinarySystem* >(input_params[0]);
        Core::EvolModeType evol_mode = *static_cast< Core::EvolModeType* >(
            input_params[1]
        );
        return system.differential_equations(age,
                                             parameters,
                                             evol_mode,
                                             derivatives);
    }

#ifdef ENABLE_DERIVATIVES
    int stellar_system_jacobian(double age, const double *orbital_parameters,
            double *param_derivs, double *age_derivs,void *system_mode)
    {
        void **input_params = static_cast<void **>(system_mode);
        BinarySystem &system = *static_cast< BinarySystem* >(input_params[0]);
        Core::EvolModeType evol_mode = *static_cast< Core::EvolModeType* >(
            input_params[1]
        );
        return system.jacobian(age,
                               orbital_parameters,
                               evol_mode,
                               param_derivs,
                               age_derivs);
    }
#endif

#ifndef NDEBUG
    void OrbitSolver::output_history_and_discarded(std::ostream &os)
    {
        assert(
            !__stop_history_ages.empty()
            ||
            !__discarded_stop_ages.empty()
        );
    //	return;
        std::streamsize orig_precision=os.precision();
        os.precision(16);
        std::ios_base::fmtflags orig_flags=os.flags();
        os.setf(std::ios_base::scientific);

        os << "Stored stop condition information:" << std::endl
            << std::setw(20) << "Age:";
        bool skip_check = true;
        std::list<double>::const_iterator age_i;
        if(__stop_history_ages.size() > 5) {
            age_i = __stop_history_ages.end();
            std::advance(age_i, -5);
        } else
            age_i = __stop_history_ages.begin();
        for(
                ;
                skip_check || age_i != __discarded_stop_ages.end();
                ++age_i
        ) {
            if(skip_check && age_i == __stop_history_ages.end()) {
                os << "|";
                age_i = __discarded_stop_ages.begin();
                skip_check=false;
            }
            os << std::setw(28) << *age_i;
        }
        std::string hline;
        hline.assign(
            20
            +
            25 * (std::min(static_cast<int>(__stop_history_ages.size()), 5)
                  +
                  __discarded_stop_ages.size()),
            '_'
        );
        os << std::endl << hline << std::endl;

        for(size_t i = 0; i < __stop_cond_discarded.front().size(); i++) {
            assert(
                __stopping_conditions->expected_crossing_deriv_sign(i) > 0
                ||
                __stopping_conditions->expected_crossing_deriv_sign(i) < 0
            );
            os << std::setw(13) << "Condition["
                << std::setw(5) << i
                << "]:"
                << (__stopping_conditions->expected_crossing_deriv_sign(i) > 0
                    ? " / "
                    : " \\ ");
            size_t point_num=0;
            std::list<double>::const_iterator age_i=__stop_history_ages.begin();
            bool marked_skip_extremum=false;
            skip_check = true;
            std::list< std::valarray<double> >::const_iterator cond_i;
            if(__stop_cond_history.size() > 5) {
                cond_i = __stop_cond_history.end();
                std::advance(cond_i, -5);
            } else
                cond_i = __stop_cond_history.begin();
            for(
                    ;
                    skip_check || cond_i != __stop_cond_discarded.end();
                    cond_i++
            ) {
                if(skip_check && cond_i == __stop_cond_history.end()) {
                    os << "|";
                    cond_i=__stop_cond_discarded.begin();
                    age_i = __discarded_stop_ages.begin();
                    skip_check = false;
                }
                bool marked=false;
                if(point_num==__skip_history_zerocrossing[i]) {
                    os << "z"; marked=true;
                } else os << " ";
                if((*age_i)>__skip_history_extremum[i]
                        && !marked_skip_extremum) {
                    os << "e"; marked_skip_extremum=true; marked=true;
                } else os << (marked ? "-" : " ");
                os << (marked ? ">" : " ");
                os << std::setw(25) << (*cond_i)[i];
                point_num++;
                age_i++;
            }
            os << std::endl;
            os << std::setw(13) << "Derivative[" << std::setw(5) << i
                << "]:";
            skip_check=true;
            std::list< std::valarray<double> >::const_iterator deriv_i;
            if(__stop_deriv_history.size() > 5) {
                deriv_i = __stop_deriv_history.end();
                std::advance(deriv_i, -5);
            } else
                deriv_i = __stop_deriv_history.begin();
            for(
                    ;
                    skip_check || deriv_i != __stop_deriv_discarded.end();
                    deriv_i++
            ) {
                if(skip_check && deriv_i == __stop_deriv_history.end()) {
                    os << "|";
                    deriv_i=__stop_deriv_discarded.begin();
                    skip_check=false;
                }
                os << std::setw(28) << (*deriv_i)[i];
            }
            os << std::endl;
        }

        os.precision(orig_precision);
        os.flags(orig_flags);
    }
#endif

    void OrbitSolver::clear_discarded()
    {
        __discarded_stop_ages.clear();
        __stop_cond_discarded.clear();
        __stop_deriv_discarded.clear();
    }

    void OrbitSolver::insert_discarded(double age,
            const std::valarray<double> &current_stop_cond,
            const std::valarray<double> &current_stop_deriv)
    {
        std::list<double>::iterator age_i = __discarded_stop_ages.begin();
        std::list< std::valarray<double> >::iterator
            cond_i = __stop_cond_discarded.begin(),
            deriv_i = __stop_deriv_discarded.begin();
        while(age_i != __discarded_stop_ages.end() && age > (*age_i)) {
            age_i++; cond_i++; deriv_i++;
        }
        __discarded_stop_ages.insert(age_i, age);
        __stop_cond_discarded.insert(cond_i, current_stop_cond);
        __stop_deriv_discarded.insert(deriv_i, current_stop_deriv);
    }

    void OrbitSolver::add_to_evolution(double age,
                                       Core::EvolModeType evolution_mode,
                                       BinarySystem &system)
    {
        clear_discarded();
        system.add_to_evolution();
        __tabulated_ages.push_back(age);
        __tabulated_evolution_modes.push_back(evolution_mode);
    }

    double OrbitSolver::go_back(double max_age,
                                BinarySystem &system,
                                std::valarray<double> &orbit)
    {
        unsigned nsteps = 0;
        while(!__tabulated_ages.empty() && max_age < __tabulated_ages.back()) {
            __tabulated_ages.pop_back();
            __tabulated_evolution_modes.pop_back();
            ++nsteps;
        }
        system.rewind_evolution(nsteps);

        if(max_age < __stop_history_ages.back()) clear_discarded();
        while(
            !__stop_history_ages.empty()
            &&
            max_age < __stop_history_ages.back()
        ) {
            __stop_history_ages.pop_back();
            __orbit_history.pop_back();
            __orbit_deriv_history.pop_back();
            __stop_cond_history.pop_back();
            __stop_deriv_history.pop_back();
        }
        for(size_t i = 0; i < __skip_history_zerocrossing.size(); ++i) {
            __skip_history_zerocrossing[i] = std::min(
                __skip_history_zerocrossing[i],
                __stop_history_ages.size()
            );
            __skip_history_extremum[i] = std::min(
                __stop_history_ages.back(),
                __skip_history_extremum[i]
            );
        }
        orbit = __orbit_history.back();
        return __stop_history_ages.back();
    }

    void OrbitSolver::clear_history()
    {
        __orbit_history.clear();
        __orbit_deriv_history.clear();
        __stop_history_ages.clear();
        __stop_cond_history.clear();
        __stop_deriv_history.clear();
        clear_discarded();
    }

    StopHistoryInterval OrbitSolver::select_stop_condition_interval(
        bool crossing,
        size_t cond_ind,
        size_t max_points
    ) const
    {
        if(crossing)
            assert(
                __stop_cond_history.size()
                -
                __skip_history_zerocrossing[cond_ind]
                >=
                1
            );
        else
            assert(
                __stop_cond_history.size()
                -
                __skip_history_zerocrossing[cond_ind]
                >=
                2
            );
        size_t num_points = std::min(
            max_points,
            (
                __stop_history_ages.size()
                +
                __discarded_stop_ages.size()
                -
                (crossing ? __skip_history_zerocrossing[cond_ind] : 0)
            )
        );
        std::list<double>::const_iterator first_age = __stop_history_ages.end();
        std::list< std::valarray<double> >::const_iterator
            first_stop_cond = __stop_cond_history.end(),
            first_stop_deriv = __stop_deriv_history.end();
        int go_back = (static_cast<int>(num_points)
                       -
                       static_cast<int>(__discarded_stop_ages.size()));
        go_back = std::max(go_back, (crossing ? 1 : 2));
        size_t failed_back = 0;
        for(int i = 0; i < go_back; ++i) {
            --first_age;
            --first_stop_cond;
            --first_stop_deriv;
            if(!crossing && *first_age < __skip_history_extremum[cond_ind]) {
                ++failed_back;
                ++first_age;
                ++first_stop_cond;
                ++first_stop_deriv;
                --num_points;
            }
        }
        if(go_back - failed_back < (crossing ? 1 : 2))
            return StopHistoryInterval();
        StopHistoryInterval interval(num_points,
                                     first_age,
                                     __stop_history_ages.end(),
                                     __discarded_stop_ages.begin(),
                                     first_stop_cond,
                                     __stop_cond_history.end(),
                                     __stop_cond_discarded.begin(),
                                     first_stop_deriv,
                                     __stop_deriv_history.end(),
                                     __stop_deriv_discarded.begin()),
                            result = interval;
        int max_left_shift = std::min(__discarded_stop_ages.size() - 1,
                                      num_points - (crossing ? 2 : 3));
        int history_limit = 0;
        if(crossing)
            history_limit = (__stop_history_ages.size()
                             -
                             go_back
                             +
                             failed_back
                             -
                             __skip_history_zerocrossing[cond_ind]);
        else
            while(
                first_age != __stop_history_ages.begin()
                &&
                *(--first_age) > __skip_history_extremum[cond_ind]
            )
                ++history_limit;
        max_left_shift = std::min(history_limit, max_left_shift);
        for(int i = 0; i < max_left_shift; i++) {
            interval << 1;
            if(
                (interval.last_age() - interval.first_age())
                <
                result.last_age() - result.first_age()
            )
                result = interval;
        }
        return result;
    }

    ExtremumInformation OrbitSolver::extremum_from_history_no_deriv(
        size_t condition_index
    ) const
    {
        ExtremumInformation result;
        if(
            __stop_history_ages.size()
            -
            __skip_history_zerocrossing[condition_index] < 2
        )
            return result;
        std::list< std::valarray<double> >::const_iterator stop_cond_i =
            __stop_cond_history.end();
        double pre1_cond = (*(--stop_cond_i))[condition_index],
               pre2_cond = (*(--stop_cond_i))[condition_index],
               post_cond = __stop_cond_discarded.front()[condition_index];
        if(
            std::abs(pre1_cond) > std::abs(pre2_cond)
            ||
            std::abs(pre1_cond) > std::abs(post_cond)
            ||
            (
                !std::isfinite(pre1_cond)
                &&
                !std::isfinite(pre2_cond)
                &&
                !std::isfinite(post_cond)
            )
        )
            return result;

        StopHistoryInterval stop_interval = select_stop_condition_interval(
            false,
            condition_index,
            4
        );
        if(stop_interval.num_points() < 3)
            return result;
        double t0 = stop_interval.age(),
               c0 = stop_interval.stop_condition_value(condition_index),
               t1 = (++stop_interval).age(),
               c1 = stop_interval.stop_condition_value(condition_index),
               t2 = (++stop_interval).age(),
               c2 = stop_interval.stop_condition_value(condition_index),
               abs_c0 = std::abs(c0),
               abs_c1 = std::abs(c1),
               abs_c2 = std::abs(c2);

        const double min_fractional_diff = (
            1000.0
            *
            std::numeric_limits<double>::epsilon()
        );

        bool ignore_10_diff = (std::abs(c1 - c0) / std::max(abs_c0, abs_c1)
                               <
                               min_fractional_diff),
             ignore_21_diff = (std::abs(c2 - c1) / std::max(abs_c1, abs_c2)
                               <
                               min_fractional_diff);


        if(stop_interval.num_points() == 3) {
            if((c1 - c0) * (c2 - c1) > 0 || ignore_10_diff || ignore_21_diff)
            {
                return result;
            }
            result.x() = Core::quadratic_extremum(t0, c0, t1, c1, t2, c2,
                                                  &(result.y()));
        } else {
            double t3 = (++stop_interval).age(),
                   c3 = stop_interval.stop_condition_value(condition_index),
                   abs_c3 = std::abs(c3);
            bool ignore_32_diff = (
                std::abs(c3 - c2) / std::max(abs_c3, abs_c2)
                <
                min_fractional_diff
            );
            if(
                (
                    (c1 - c0) * (c2 - c1) > 0
                    ||
                    abs_c1 >= std::max(abs_c0, abs_c2)
                    ||
                    ignore_10_diff
                    ||
                    ignore_21_diff
                )
                &&
                (
                    (c2 - c1) * (c3 - c2) > 0
                    ||
                    abs_c2 >= std::max(abs_c1, abs_c3)
                    ||
                    ignore_21_diff
                    ||
                    ignore_32_diff
                )
            ){
                return result;
            }
            double range_low,
                   range_high;
            if(
                !ignore_10_diff
                &&
                std::abs(c1) <= std::abs(c0)
                &&
                std::abs(c1) <= std::abs(c2)
            ) {
                range_low = t0;
                range_high = t2;
            } else if(
                !ignore_32_diff
                &&
                abs_c2 <= abs_c1
                &&
                abs_c2 <= abs_c3
            ) {
                range_low = t1;
                range_high = t3;
            } else throw Core::Error::BadFunctionArguments(
                "Searching for extremum among monotonic stopping condition "
                "values in OrbitSolver::extremum_from_history_no_deriv."
            );
            result.x() = Core::cubic_extremum(t0, c0, t1, c1, t2, c2, t3, c3,
                                              &(result.y()),
                                              range_low, range_high);
        }
        return result;
    }

    ExtremumInformation OrbitSolver::extremum_from_history(
        size_t condition_index,
        double min_extremum_x
    ) const
    {
        ExtremumInformation result;
        if(
            __stop_history_ages.size() == 0
            ||
            (
                __stop_history_ages.back()
                <
                __skip_history_extremum[condition_index]
            )
        )
            return result;
        double prev_stop_cond = __stop_cond_history.back()[condition_index],
               next_stop_cond = __stop_cond_discarded.front()[condition_index];
        if(next_stop_cond * prev_stop_cond <= 0)
            return result;
        if(std::isnan(__stop_deriv_history.back()[condition_index])) {
            result = extremum_from_history_no_deriv(condition_index);
        } else {
            double
                prev_stop_deriv = __stop_deriv_history.back()[
                    condition_index
                ],
                next_stop_deriv = __stop_deriv_discarded.front()[
                    condition_index
                ];
            if(
                next_stop_cond * next_stop_deriv < 0
                ||
                next_stop_deriv * prev_stop_deriv >= 0
            )
                return result;
            result.x() = Core::estimate_extremum(__stop_history_ages.back(),
                                                 prev_stop_cond,
                                                 __discarded_stop_ages.front(),
                                                 next_stop_cond,
                                                 prev_stop_deriv,
                                                 next_stop_deriv,
                                                 &(result.y()));
        }

        if(result.x() < min_extremum_x) {
            result.x() = min_extremum_x;
            result.y() = Core::NaN;
        }

        return result;
    }

    double OrbitSolver::crossing_from_history_no_deriv(size_t condition_index)
        const
    {
        StopHistoryInterval stop_interval = select_stop_condition_interval(
            true,
            condition_index,
            4
        );
        if(stop_interval.num_points() < 2) return Core::Inf;
        double t0 = stop_interval.age(),
               c0 = stop_interval.stop_condition_value(condition_index),
               t1 = (++stop_interval).age(),
               c1 = stop_interval.stop_condition_value(condition_index);
        if(stop_interval.num_points() == 2)
            return Core::estimate_zerocrossing(t0, c0, t1, c1);
        double t2 = (++stop_interval).age(),
               c2 = stop_interval.stop_condition_value(condition_index);
        double range_low = Core::NaN,
               range_high = Core::NaN;
        short crossing_sign =
            __stopping_conditions->expected_crossing_deriv_sign(
                condition_index
            );
        if(c0 * c1 <= 0 && c1 * crossing_sign > 0) {
            range_low = t0;
            range_high = t1;
        } else if(c1 * c2 <= 0 && c2 * crossing_sign > 0) {
            range_low = t1;
            range_high = t2;
        }
        if(stop_interval.num_points() == 3) {
            assert(!std::isnan(range_low) && !std::isnan(range_high));
            return Core::quadratic_zerocrossing(
                t0, c0, t1, c1, t2, c2, range_low, range_high
            );
        }
        double t3 = (++stop_interval).age(),
               c3 = stop_interval.stop_condition_value(condition_index);
        if(std::isnan(range_low)) {
            range_low = t2;
            range_high = t3;
            assert(c2 * c3 < 0);
            assert(c3 * crossing_sign > 0);
        }
        return Core::cubic_zerocrossing(
            t0, c0, t1, c1, t2, c2, t3, c3, range_low, range_high
        );
    }

    double OrbitSolver::crossing_from_history(size_t condition_index) const
    {
        if(
            __stop_history_ages.size() == 0
            ||
            (
                __stop_history_ages.size()
                ==
                __skip_history_zerocrossing[condition_index]
            )
        )
            return Core::Inf;
        double prev_stop_cond = __stop_cond_history.back()[condition_index];
        double next_stop_cond = __stop_cond_discarded.front()[condition_index];
        if(
            next_stop_cond * prev_stop_cond > 0
            ||
            (
                next_stop_cond
                *
                __stopping_conditions->expected_crossing_deriv_sign(
                    condition_index
                )
                <
                0
            )
            ||
            (
                next_stop_cond == 0
                &&
                prev_stop_cond
                *
                __stopping_conditions->expected_crossing_deriv_sign(
                    condition_index
                )
                >=
                0
            )
        )
            return Core::Inf;
        if(std::isnan(__stop_deriv_history.back()[condition_index]))
            return crossing_from_history_no_deriv(condition_index);
        double
            prev_stop_deriv = __stop_deriv_history.back()[condition_index],
            prev_age = __stop_history_ages.back(),
            next_stop_deriv = __stop_deriv_discarded.front()[condition_index],
            next_age = __discarded_stop_ages.front();
        return Core::estimate_zerocrossing(prev_age,
                                           prev_stop_cond,
                                           next_age,
                                           next_stop_cond,
                                           prev_stop_deriv,
                                           next_stop_deriv);
    }

    void OrbitSolver::initialize_skip_history(
        const StoppingCondition &stop_cond,
        const std::valarray<double> &stop_cond_values,
        StoppingConditionType stop_reason
    )
    {
#ifndef NDEBUG
        std::cerr << "Initializing skip history with stop reason: "
                  << stop_reason
                  << std::endl;
#endif
        __skip_history_zerocrossing.resize(stop_cond.num_subconditions(), 0);
        __skip_history_extremum.resize(stop_cond.num_subconditions(), 0);
        for(
            size_t cond_ind = 0;
            cond_ind < stop_cond.num_subconditions();
            cond_ind++
        ) {
            StoppingConditionType stop_cond_type=stop_cond.type(cond_ind);
            if(
                (
                    (
                        stop_reason == BREAK_LOCK
                        &&
                        stop_cond_type == SYNCHRONIZED
                    )
                    ||
                    (
                        stop_reason != NO_STOP
                        &&
                        stop_reason != BREAK_LOCK
                        &&
                        stop_cond_type == stop_reason
                    )
                )
                &&
                std::abs(stop_cond_values[cond_ind]) <= __precision
            ) {
#ifndef NDEBUG
                std::cerr << "Skipping first step of condition "
                          << cond_ind
                          << "(" << stop_cond_type << ")"
                          << std::endl;
#endif
                __skip_history_zerocrossing[cond_ind]=1;
                __skip_history_extremum[cond_ind]=(
                    __stop_history_ages.front()
                    *
                    (1.0+std::numeric_limits<double>::epsilon())
                );
            }
#ifdef VERBOSE_DEBUG
            else
                std::cerr << "Not skipping condition "
                          << cond_ind
                          << "(" << stop_cond_type << ")"
                          << std::endl;
#endif
        }
    }

    void OrbitSolver::update_skip_history(
            const std::valarray<double> &current_stop_cond,
            const StopInformation &stop_info)
    {
        size_t history_size = __stop_history_ages.size();
        for(size_t i = 0; i < current_stop_cond.size(); i++) {
            if(
                    __skip_history_zerocrossing[i] == history_size
                    &&
                    std::abs(current_stop_cond[i]) <= __precision
            )
                ++__skip_history_zerocrossing[i];
            if(
                    stop_info.stop_reason() == __stopping_conditions->type(i)
                    &&
                    !stop_info.is_crossing()
                    &&
                    (
                        std::abs(stop_info.stop_condition_precision())
                        <=
                        __precision
                    )
            )
                __skip_history_extremum[i]=stop_info.stop_age();
        }
    }

    bool OrbitSolver::at_exact_condition(double previous_age,
                                         const StopInformation &stop_info)
    {
        return (
            (
                std::abs(stop_info.stop_condition_precision())
                <=
                __precision
            )
            ||
            (
                stop_info.stop_age()
                <
                previous_age * MIN_RELATIVE_STEP
            )
        );
    }

    bool OrbitSolver::acceptable_step(double current_age,
                                      double previous_age,
                                      const StopInformation &stop_info)
    {
#ifdef VERBOSE_DEBUG
        std::cerr << "From t = " << previous_age
                  << ", stepped to t = " << current_age
                  << ", stop at t = " << stop_info.stop_age()
                  << ", must be at least: "  << previous_age * MIN_RELATIVE_STEP
                  << std::endl;
        if(
            at_exact_condition(previous_age, stop_info)
            &&
            std::abs(stop_info.stop_condition_precision()) > __precision
        ) {
            std::cerr << "Failed to meet precision for "
                      << stop_info
                      << std::endl;
        }
#endif
        return (
            stop_info.stop_age() >= current_age
            ||
            (
                at_exact_condition(previous_age, stop_info)
                &&
                (stop_info.crossed_zero() || !stop_info.is_crossing())
            )
        );
    }

    StopInformation OrbitSolver::update_stop_condition_history(
        double age,
        const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        Core::EvolModeType evolution_mode,
        unsigned current_expansion_order,
        StoppingConditionType stop_reason
    )
    {
        for(unsigned i = 0; i < orbit.size(); ++i)
            if(
                std::isnan(orbit[i])
                ||
                (evolution_mode == Core::BINARY && orbit[0] <= 0)
            ) {
#ifndef NDEBUG
                std::cerr << "Bad orbit: " << orbit << std::endl;
#endif
                return StopInformation(
                    0.5 * (age + __stop_history_ages.back()),
                    Core::Inf
                );
            }

        std::pair<double, double> expansion_range =
            TidalPotentialTerms::get_expansion_range(current_expansion_order);

#ifdef VERBOSE_DEBUG
        if(evolution_mode == Core::BINARY)
            std::cerr << "Updating stop condition history. Current e = "
                      << orbit[1]
                      << " current expansion order: "
                      << current_expansion_order
                      << " current expansion range: "
                      << expansion_range.first
                      << " < e < "
                      << expansion_range.second
                      << std::endl;
#endif


        if(
                evolution_mode == Core::BINARY
                &&
                orbit[1] > expansion_range.second
        ) {
#ifdef VERBOSE_DEBUG
            std::cerr << "Eccentricity ("
                      << orbit[1]
                      << ") exceeds current expansion error limit of "
                      << expansion_range.second
                      << ". Choosing to stop half way between t = "
                      << age
                      << " and t = "
                      << __stop_history_ages.back()
                      << std::endl;
#endif
            return StopInformation(
                0.5 * (age + __stop_history_ages.back()),
                Core::Inf,
                LARGE_EXPANSION_ERROR,
                true,
                true
            );
        }

        std::valarray<double> current_stop_cond(
            __stopping_conditions->num_subconditions()
        );
        std::valarray<double> current_stop_deriv;
        current_stop_cond = (*__stopping_conditions)(evolution_mode,
                                                     orbit,
                                                     derivatives,
                                                     current_stop_deriv);

        if(__stop_history_ages.empty())
            initialize_skip_history(*__stopping_conditions,
                                    current_stop_cond,
                                    stop_reason);
        StopInformation result;
        insert_discarded(age, current_stop_cond, current_stop_deriv);
#ifdef VERBOSE_DEBUG
    	std::cerr << std::string(77, '@') << std::endl;
    	output_history_and_discarded(std::cerr);
        std::cerr << "Decreasing expansion order is "
                  << (current_expansion_order == 0 ? "not" : "")
                  << " allowed"
                  <<std::endl;
#endif
        for(
            size_t cond_ind = 0;
            (
                __stop_cond_history.size() > 0
                &&
                cond_ind < current_stop_cond.size()
            );
            cond_ind++
        ) {
            double stop_cond_value = current_stop_cond[cond_ind],
                   crossing_age = crossing_from_history(cond_ind),
                   crossing_precision =
                       __stop_cond_history.back()[cond_ind];
            bool crossed_zero = false;
            if(std::abs(crossing_precision) >= std::abs(stop_cond_value)) {
                crossing_precision = stop_cond_value;
                crossed_zero = true;
            }
            ExtremumInformation extremum = extremum_from_history(
                cond_ind,
                __stop_history_ages.back() + __min_extremum_search_step
            );
            double extremum_precision;
            if(std::isnan(extremum.y())) extremum_precision = Core::NaN;
            else
                extremum_precision = (
                    std::min(
                        std::abs(extremum.y() - stop_cond_value),
                        std::abs(extremum.y()
                                 -
                                 __stop_cond_history.back()[cond_ind])
                    )
                    /
                    std::abs(extremum.y())
                );
            bool is_crossing = crossing_age <= extremum.x();
            short deriv_sign = 0;
            if(is_crossing) deriv_sign = (stop_cond_value > 0 ? 1 : -1);
            StopInformation &stop_info = __stop_info[cond_ind];
            stop_info.stop_age() = std::min(crossing_age, extremum.x());
            stop_info.stop_condition_precision() = (is_crossing
                                                    ? crossing_precision
                                                    : extremum_precision);
            stop_info.is_crossing() = is_crossing;
            stop_info.crossed_zero() = crossed_zero;
            stop_info.deriv_sign_at_crossing() = (is_crossing
                                                  ? deriv_sign
                                                  : 0.0);
#ifdef VERBOSE_DEBUG
            std::cerr << "Condition " << cond_ind << " "
                      << __stopping_conditions->describe(cond_ind)
                      << ": " << stop_info
                      << std::endl;
#endif
            if(
                !__stop_history_ages.empty()
                &&
                (
                    !acceptable_step(age, __stop_history_ages.back(), stop_info)
                    ||
                    is_crossing
                )
                &&
                stop_info.stop_age() < result.stop_age()
            ) {
                result = stop_info;
#ifdef VERBOSE_DEBUG
                std::cerr << "SELECTED" << std::endl;
            } else {
                std::cerr << "NOT SELECTED!" << std::endl;
#endif
            }
        }
        if(
            __stop_history_ages.empty()
            ||
            acceptable_step(age, __stop_history_ages.back(), result)
        ) {
            if(!__stop_history_ages.empty())
                update_skip_history(current_stop_cond, result);
            __stop_history_ages.push_back(age);
            __stop_cond_history.push_back(current_stop_cond);
            __stop_deriv_history.push_back(current_stop_deriv);
            __orbit_history.push_back(orbit);
            __orbit_deriv_history.push_back(derivatives);

            if(
                result.stop_reason() == NO_STOP
                &&
                evolution_mode == Core::BINARY
                &&
                orbit[1] < expansion_range.first
                &&
                age > (0.99 * __last_order_upgrade_age
                       +
                       0.01 * __end_age)
                &&
                current_expansion_order > 0
            ) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Eccentricity ("
                          << orbit[1]
                          << " is sufficiently small (< "
                          << expansion_range.first
                          << ") to decrease expansion order."
                          << std::endl;
#endif
                return StopInformation(Core::Inf,
                                       Core::Inf,
                                       SMALL_EXPANSION_ERROR,
                                       true,
                                       true);
            }

        } else {
#ifndef NDEBUG
    		std::cerr << "Step to age = "
                      << age
                      << " deemed unacceptable: "
                      << result
                      << std::endl;
#endif
        }
        return result;
    }

    void OrbitSolver::reject_step(
        double &t,
        StopInformation &stop,
        BinarySystem &system,
        std::valarray<double> &orbit,
        double &max_next_t,
        double &step_size
#ifndef NDEBUG
        , std::string reason
#endif
    )
    {
        double last_good_t = go_back(stop.stop_age(),
                                     system,
                                     orbit);
#ifndef NDEBUG
        std::cerr
            << "Reverting step from t = "
            << t
            << " to "
            << last_good_t
            << " due to "
            << reason
            << "Stop info: "
            << stop
            << std::endl;
#endif

        if(
            t
            <
            last_good_t * MIN_RELATIVE_STEP
        ) {
#ifndef NDEBUG
            std::cerr << "Stepped only "
                << t - last_good_t
                << "Gyr, aborting!"
                << std::endl;
#endif
            throw Core::Error::NonGSLZeroStep();
        }
        t = last_good_t;
        if(stop.is_crossing())
            stop.stop_condition_precision() = (
                __stop_cond_history.back()[
                stop.stop_condition_index()
                ]
            );
        max_next_t = stop.stop_age();
        step_size = 0.1 * (max_next_t - t);
    }

    StopInformation OrbitSolver::evolve_until(
        BinarySystem &system,
        double &max_age,
        std::valarray<double> &orbit,
        StoppingConditionType &stop_reason,
        double max_step,
        Core::EvolModeType evolution_mode
    )
    {
        size_t nargs = orbit.size();
#ifndef NDEBUG
        std::cerr << "Starting evolution leg in " << evolution_mode
            << " from t=" << system.age() << " with initial orbit:\n";
        for(size_t i = 0; i < nargs; ++i) {
            if(i) std::cerr << ", ";
            std::cerr << "\t" << orbit[i] << std::endl;
        }
        std::cerr << std::endl;
        std::cerr << "Stopping conditions:" << std::endl
                  << __stopping_conditions->describe() << std::endl;
#endif

        const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;

        gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(step_type, nargs);
        gsl_odeiv2_control *step_control = gsl_odeiv2_control_standard_new(
            __precision,
            __precision,
            1,
            0
        );
        gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(nargs);

        void *sys_mode[2]={&system, &evolution_mode};
        gsl_odeiv2_system ode_system={stellar_system_diff_eq,
#ifdef ENABLE_DERIVATIVES
                                      stellar_system_jacobian,
#else
                                      NULL,
#endif
                                      nargs,
                                      sys_mode};
        double t=system.age();
        std::valarray<double> derivatives(nargs),
                              param_derivatives(nargs),
                              age_derivatives(nargs);

        add_to_evolution(t, evolution_mode, system);

        clear_discarded();
        double step_size = std::min(0.1 * (max_age - t),
                                    max_step);

        StopInformation stop;
        bool first_step = true;
        double from_t;
        while(t<max_age) {
            double max_next_t = std::min(t + max_step, max_age);
            int status=GSL_SUCCESS;
            bool step_rejected=false;
            do {
#ifndef NDEBUG
                std::cerr << "Attempting step from t = " << t
                          << " not to miss t = " << max_next_t
                          << ", suggested step = " << step_size
                          << ", orbit:\n";
                for(size_t i=0; i<nargs; ++i) {
                    if(i) std::cerr << ", ";
                    std::cerr << "\t" << orbit[i] << std::endl;
                }
                std::cerr << std::endl;
#endif
                from_t = t;
                if(!first_step) {
                    step_size = std::max(step_size,
                                         3.0 * (MIN_RELATIVE_STEP * t - t));
                    status = gsl_odeiv2_evolve_apply(evolve,
                                                     step_control,
                                                     step,
                                                     &ode_system,
                                                     &t,
                                                     max_next_t,
                                                     &step_size,
                                                     &(orbit[0]));
                }
                if (status == GSL_FAILURE) {
#ifndef NDEBUG
                    std::cerr << "Failed, (presume zero step size)!"
                              << std::endl;
#endif
                    throw Core::Error::GSLZeroStep("rkf45");
                } else if (status != GSL_SUCCESS && status != GSL_EDOM) {
                    std::ostringstream msg;
                    msg << "GSL signaled failure while evolving (error code " <<
                        status << ")";
                    throw Core::Error::Runtime(msg.str());
                } else if (
                    (
                        __runtime_limit > 0
                        &&
                        difftime(time(NULL), __evolution_start_time)
                        >
                        __runtime_limit
                    )
                    ||
                    (
                        __num_step_limit > 0
                        &&
                        (
                            __tabulated_ages.size()
                            +
                            __discarded_stop_ages.size()
                            >
                            __num_step_limit
                        )
                    )
                ) {
                    std::ostringstream msg;
                    msg << "After "
                        << __tabulated_ages.size()
                        << "steps (+"
                        << __discarded_stop_ages.size()
                        << " discarded), exceeded evolution time limit of "
                        << __runtime_limit
                        << " seconds or step limit of "
                        << __num_step_limit
                        << " steps!";
                    throw Core::Error::Runtime(msg.str());
                }
                if(status == GSL_SUCCESS) {
#ifdef NDEBUG
                    if(__print_progress)
#endif
                        std::cerr << "Succeeded! Now t = " << t << std::endl;
#ifdef VERBOSE_DEBUG
                    std::cerr << "GSL suggested new step size:"
                              << step_size
                              << std::endl;
#endif

                    stellar_system_diff_eq(t,
                                           &(orbit[0]),
                                           &(derivatives[0]),
                                           sys_mode);

                    stop = update_stop_condition_history(
                        t,
                        orbit,
                        derivatives,
                        evolution_mode,
                        system.expansion_order(),
                        stop_reason
                    );
                    stop_reason = NO_STOP;
                }
                if(
                    (status == GSL_EDOM || !acceptable_step(t, from_t, stop))
                    &&
                    !first_step
                ) {
                    reject_step(
                        t,
                        stop,
                        system,
                        orbit,
                        max_next_t,
                        step_size
#ifndef NDEBUG
                        , (status == GSL_EDOM ? "EDOM error" : "bad step")
#endif
                    );
                    step_rejected = true;
                    gsl_odeiv2_evolve_reset(evolve);
                } else {
                    if(!first_step && t < from_t * MIN_RELATIVE_STEP) {
#ifndef NDEBUG
                        std::cerr << "Stepped only "
                                  << t - from_t
                                  << "Gyr, aborting!"
                                  << std::endl;
#endif
                        throw Core::Error::GSLZeroStep("rkf45");
                    }
                    step_rejected=false;
                }

            } while(
                step_rejected
                &&
                !first_step
                &&
                (
                    !at_exact_condition(from_t, stop)
                    ||
                    __stop_history_ages.size() == 1
                )
                &&
                stop.stop_reason() != LARGE_EXPANSION_ERROR
            );

            if(!step_rejected) {
#ifndef NDEBUG
                std::cerr << "Stepped to t = " << t << std::endl;
#endif
                add_to_evolution(t, evolution_mode, system);
            }
#ifndef NDEBUG
            std::cerr << "Stop: " << stop
                      << "expansion order: " << system.expansion_order()
                      << "last order upgrade at t=" << __last_order_upgrade_age
                      << "max age: " << max_age
                      << std::endl;
#endif
            if(
                (stop.is_crossing() && stop.stop_reason() != NO_STOP)
                ||
                stop.stop_reason() == LARGE_EXPANSION_ERROR
                ||
                (
                    stop.stop_reason() == SMALL_EXPANSION_ERROR
                    &&
                    system.expansion_order() > 0
                    &&
                    !first_step
                )
            ) {
                stop_reason = stop.stop_reason();
#ifndef NDEBUG
                std::cerr << "Breaking for = " << stop << std::endl;
#endif
                break;
            }

            first_step = false;
        }
        max_age=t;
        clear_history();
        stellar_system_diff_eq(t,
                               &(orbit[0]),
                               &(derivatives[0]),
                               sys_mode);

        gsl_odeiv2_evolve_free(evolve);
        gsl_odeiv2_control_free(step_control);
        gsl_odeiv2_step_free(step);
        return stop;
    }

    CombinedStoppingCondition *OrbitSolver::get_stopping_condition(
            BinarySystem &system
    )
    {
        CombinedStoppingCondition *result = system.stopping_conditions();
#ifdef EXTERNAL_CONDITION
        (*result) |= new EXTERNAL_CONDITION;
#endif
        __stop_info.clear();
        __stop_info.resize(result->num_subconditions());
        for(size_t cond_ind = 0; cond_ind < __stop_info.size(); ++cond_ind)
        {
            __stop_info[cond_ind].stop_reason() = result->type(cond_ind);
            __stop_info[cond_ind].stop_condition_index() = cond_ind;
        }
        return result;
    }

    double OrbitSolver::stopping_age(double age,
                                     const BinarySystem &system,
                                     const std::list<double> &required_ages)
    {
#ifndef NDEBUG
        std::cerr << "Determining next stop age: " << std::endl;
#endif
        double result = system.next_stop_age();
#ifndef NDEBUG
        std::cerr << "Next system stop age: " << result << std::endl;
#endif
        if(required_ages.size() == 0) return result;

        static std::list<double>::const_iterator
            next_required_age = required_ages.begin();
        if(age <= required_ages.front())
            next_required_age = required_ages.begin();
        if(
            next_required_age != required_ages.end()
            &&
            age == *next_required_age
        )
            ++next_required_age;
        if(
            next_required_age != required_ages.end()
            &&
            result > *next_required_age
        )
            result = *next_required_age;
#ifndef NDEBUG
        std::cerr << "Required ages change that to: " << result << std::endl;
#endif
        return result;
    }

    void OrbitSolver::reached_stopping_condition(
        double stop_age,
        StoppingConditionType stop_reason
    )
    {
#ifndef NDEBUG
        std::cerr << "Stopped due to condition at t = "
                  << stop_age
                  << std::endl;
#endif
        for(
            std::vector<StopInformation>::const_iterator stop_i =
                __stop_info.begin();
            stop_i != __stop_info.end();
            ++stop_i
        ) {
            if(
                stop_i->is_crossing()
                &&
                (
                    stop_i->stop_age() < stop_age
                    ||
                    (
                        stop_i->stop_reason() == stop_reason
                        &&
                        at_exact_condition(stop_age, *stop_i)
                    )
                )
            ) {
#ifndef NDEBUG
                std::cerr << "Triggered condition: "
                          << __stopping_conditions->describe(
                              stop_i->stop_condition_index()
                          )
                          << std::endl;
#endif
                __stopping_conditions->reached(
                    stop_i->deriv_sign_at_crossing(),
                    stop_i->stop_condition_index()
                );
            }
        }
    }

    ///TODO: can be dramatically simplified
    void OrbitSolver::adjust_expansion_order(
        BinarySystem &system,
        const std::valarray<double> &orbit,
        Core::EvolModeType evolution_mode,
        bool must_increase
    )
    {
#ifndef NDEBUG
        std::cerr << "Adjusting expansion order at t ="
                  << system.age();
        if(must_increase)
            std::cerr << " upward!";
        std::cerr << std::endl;
#endif
        assert(evolution_mode == Core::BINARY);

        unsigned
            current_expansion_order = system.expansion_order(),
            required_expansion_order =
                TidalPotentialTerms::required_expansion_order(orbit[1]);

        if(must_increase)
            required_expansion_order = std::max(required_expansion_order,
                                                current_expansion_order + 1);

        if(required_expansion_order != current_expansion_order) {
            if(required_expansion_order > current_expansion_order)
                __last_order_upgrade_age = system.age();
            system.change_expansion_order(required_expansion_order);
        }

#ifdef NDEBUG
        if(__print_progress)
#endif
            std::cerr << "At e(t = "
                      << system.age()
                      << ") = "
                      << orbit[1]
                      << " adjusted expansion order from "
                      << current_expansion_order
                      << " to "
                      << required_expansion_order
                      << std::endl;
    }

    void OrbitSolver::reset(BinarySystem &system)
    {
        __tabulated_ages.clear();
        __tabulated_evolution_modes.clear();
        system.reset_evolution();
    }

    OrbitSolver::OrbitSolver(double max_age,
                             double required_precision,
                             bool print_progress) :
        __end_age(max_age),
        __precision(required_precision),
        __last_order_upgrade_age(-Core::Inf),
        __stopping_conditions(NULL),
        __print_progress(print_progress)
    {
#ifndef NDEBUG
        assert(max_age>0);
#endif
    }

    void OrbitSolver::operator()(BinarySystem &system,
                                 double max_step,
                                 const std::list<double> &required_ages,
                                 double max_runtime,
                                 unsigned max_time_steps,
                                 double min_extremum_search_step)
    {
#ifndef NDEBUG
        std::cerr << "Calculating evolution from t = " << system.age()
                  << " to t = " << __end_age << std::endl;
#endif
        time(&__evolution_start_time);
        __runtime_limit = max_runtime;
        __num_step_limit = max_time_steps;

        __min_extremum_search_step = min_extremum_search_step;

        double stop_evol_age = __end_age;

        reset(system);
        StoppingConditionType stop_reason = NO_STOP;
        double last_age = system.age();
        std::valarray<double> orbit;

        Core::EvolModeType evolution_mode = system.fill_orbit(orbit);

        if(evolution_mode == Core::BINARY) {
            adjust_expansion_order(system,
                                   orbit,
                                   evolution_mode);
            __last_order_upgrade_age = -Core::Inf;
        }

        while(last_age < stop_evol_age) {
            double next_stop_age = std::min(stopping_age(last_age,
                                                         system,
                                                         required_ages),
                                            stop_evol_age);
            __stopping_conditions = get_stopping_condition(system);
#ifndef NDEBUG
            std::cerr << "Next stop age: " << next_stop_age << std::endl;
            StopInformation stop_information =
#endif
                evolve_until(system,
                             next_stop_age,
                             orbit,
                             stop_reason,
                             max_step,
                             evolution_mode);

            Core::EvolModeType old_evolution_mode = evolution_mode;
#ifndef NDEBUG
            std::cerr << "Stop information: "
                      << stop_information
                      << std::endl;
            unsigned old_locked_zones = system.number_locked_zones();
#endif
            last_age = next_stop_age;
            if(last_age < stop_evol_age) {
                if(stop_reason == NO_STOP) {
                    system.reached_critical_age(last_age);
                } else if(
                    stop_reason == LARGE_EXPANSION_ERROR
                    ||
                    stop_reason == SMALL_EXPANSION_ERROR
                ) {
                    adjust_expansion_order(
                        system,
                        orbit,
                        evolution_mode,
                        stop_reason == LARGE_EXPANSION_ERROR
                    );
                } else
                    reached_stopping_condition(last_age, stop_reason);
            }
            evolution_mode = system.evolution_mode();
#ifndef NDEBUG
            std::valarray<double> old_orbit(orbit);
#endif
            system.fill_orbit(orbit);

            if(evolution_mode == Core::BINARY) {
                if(old_evolution_mode != Core::BINARY) {
                    adjust_expansion_order(system,
                                           orbit,
                                           evolution_mode);
                    __last_order_upgrade_age = -Core::Inf;
                }
                if(
                        stop_reason != SMALL_EXPANSION_ERROR
                        &&
                        stop_reason != BREAK_LOCK
                        &&
                        stop_reason != WIND_SATURATION
                ) {
                    system.initialize_locks(__precision);
                    system.fill_orbit(orbit);
                }
            }

#ifndef NDEBUG
            std::cerr
                << "At t=" << last_age
                << ", changing evolution mode from " << old_evolution_mode
                << " with " << old_locked_zones
                << " zones locked to " << evolution_mode
                << " with " << system.number_locked_zones()
                << " zones locked."
                << std::endl
                << "Transforming orbit from: " << old_orbit
                << " to " << orbit << std::endl;
#endif
            delete __stopping_conditions;
            __stopping_conditions = NULL;
        }
    }

} //End Evolve namespace.
