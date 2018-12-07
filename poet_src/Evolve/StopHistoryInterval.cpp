/**\file
 *
 * \brief Declares some of the methods of the StopHistoryInterval class.
 *
 * \ingroup OrbitSolver_group
 */

#define BUILDING_LIBRARY
#include "StopHistoryInterval.h"

namespace Evolve {

    void StopHistoryInterval::advance_iterator_set(
            std::list<double>::const_iterator &age_i,
            std::list< std::valarray<double> >::const_iterator &cond_i,
            std::list< std::valarray<double> >::const_iterator &deriv_i)
    {
        ++age_i; ++cond_i; ++deriv_i;
        if(age_i == __history_age_end) {
            assert(cond_i == __stop_cond_history_end
                   &&
                   deriv_i == __stop_deriv_history_end);
            age_i = __discarded_age_begin;
            cond_i = __stop_cond_discarded_begin;
            deriv_i = __stop_deriv_discarded_begin;
        }
    }

    void StopHistoryInterval::retreat_iterator_set(
            std::list<double>::const_iterator &age_i,
            std::list< std::valarray<double> >::const_iterator &cond_i,
            std::list< std::valarray<double> >::const_iterator &deriv_i)
    {
        if(age_i == __discarded_age_begin) {
            assert(cond_i == __stop_cond_discarded_begin
                   &&
                   deriv_i == __stop_deriv_discarded_begin);
            age_i = __history_age_end;
            cond_i = __stop_cond_history_end;
            deriv_i = __stop_deriv_history_end;
        }
        --age_i; --cond_i; --deriv_i;
    }

    StopHistoryInterval::StopHistoryInterval(
        size_t num_points,
        std::list<double>::const_iterator first_age,
        std::list<double>::const_iterator history_age_end,
        std::list<double>::const_iterator discarded_age_begin,
        std::list< std::valarray<double> >::const_iterator first_stop_cond,
        std::list< std::valarray<double> >::const_iterator
        stop_cond_history_end,
        std::list< std::valarray<double> >::const_iterator
        stop_cond_discarded_begin,
        std::list< std::valarray<double> >::const_iterator first_stop_deriv,
        std::list< std::valarray<double> >::const_iterator
        stop_deriv_history_end,
        std::list< std::valarray<double> >::const_iterator
        stop_deriv_discarded_begin
    ) :
        __num_points(num_points),
        __point_i(0),
        __first_age(first_age),
        __last_age(first_age),
        __history_age_end(history_age_end),
        __discarded_age_begin(discarded_age_begin),
        __age_i(first_age),
        __first_stop_cond(first_stop_cond),
        __last_stop_cond(first_stop_cond),
        __stop_cond_history_end(stop_cond_history_end),
        __stop_cond_discarded_begin(stop_cond_discarded_begin),
        __stop_cond_i(first_stop_cond),
        __first_stop_deriv(first_stop_deriv),
        __last_stop_deriv(first_stop_deriv),
        __stop_deriv_history_end(stop_deriv_history_end),
        __stop_deriv_discarded_begin(stop_deriv_discarded_begin),
        __stop_deriv_i(first_stop_deriv)
    {
        if(num_points==0)
            throw Core::Error::BadFunctionArguments(
                "Attempt to contsruct a StopHistoryInterval of size 0."
            );
        for(size_t i=0; i<num_points-1; i++)
            advance_iterator_set(__last_age,
                                 __last_stop_cond,
                                 __last_stop_deriv);
        reset();
    }

    StopHistoryInterval::StopHistoryInterval(const StopHistoryInterval &orig) :
        __num_points(orig.__num_points),
        __point_i(orig.__point_i),
        __first_age(orig.__first_age),
        __last_age(orig.__last_age),
        __history_age_end(orig.__history_age_end),
        __discarded_age_begin(orig.__discarded_age_begin),
        __age_i(orig.__age_i),
        __first_stop_cond(orig.__first_stop_cond),
        __last_stop_cond(orig.__last_stop_cond),
        __stop_cond_history_end(orig.__stop_cond_history_end),
        __stop_cond_discarded_begin(orig.__stop_cond_discarded_begin),
        __stop_cond_i(orig.__stop_cond_i),
        __first_stop_deriv(orig.__first_stop_deriv),
        __last_stop_deriv(orig.__last_stop_deriv),
        __stop_deriv_history_end(orig.__stop_deriv_history_end),
        __stop_deriv_discarded_begin(orig.__stop_deriv_discarded_begin),
        __stop_deriv_i(orig.__stop_deriv_i)
    {}

    void StopHistoryInterval::reset()
    {
        __point_i = 0;
        __age_i = __first_age;
        __stop_cond_i = __first_stop_cond;
        __stop_deriv_i = __first_stop_deriv;
    }

    StopHistoryInterval &StopHistoryInterval::operator++()
    {
        ++__point_i;
        if(__point_i>__num_points)
            throw Core::Error::Runtime(
                "Attempting to increment two points past the end of a "
                "StopHistoryInterval!"
            );
        advance_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
        return *this;
    }

    StopHistoryInterval StopHistoryInterval::operator++(int)
    {
        StopHistoryInterval result(*this);
        operator++();
        return result;
    }

    StopHistoryInterval &StopHistoryInterval::operator--()
    {
        if(__point_i==0) throw Core::Error::Runtime(
            "Attempting to go before the beginning of a "
            "StopHistoryInterval."
        );
        --__point_i;
        retreat_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
        return *this;
    }

    StopHistoryInterval StopHistoryInterval::operator--(int)
    {
        StopHistoryInterval result(*this);
        operator--();
        return result;
    }

    StopHistoryInterval &StopHistoryInterval::operator<<(size_t n)
    {
        for(size_t i=0; i<n; i++) {
            retreat_iterator_set(__first_age, __first_stop_cond,
                    __first_stop_deriv);
            retreat_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
            retreat_iterator_set(__last_age, __last_stop_cond,__last_stop_deriv);
        }
        return *this;
    }

    StopHistoryInterval &StopHistoryInterval::operator>>(size_t n)
    {
        for(size_t i=0; i<n; i++) {
            advance_iterator_set(__last_age, __last_stop_cond, __last_stop_deriv);
            advance_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
            advance_iterator_set(__first_age, __first_stop_cond, __first_stop_deriv);
        }
        return *this;
    }

    StopHistoryInterval &StopHistoryInterval::operator=(
            const StopHistoryInterval &rhs)
    {
        __num_points=rhs.__num_points;
        __point_i=rhs.__point_i;
        __first_age=rhs.__first_age;
        __last_age=rhs.__last_age;
        __history_age_end=rhs.__history_age_end;
        __discarded_age_begin=rhs.__discarded_age_begin;
        __age_i=rhs.__age_i;
        __first_stop_cond=rhs.__first_stop_cond;
        __last_stop_cond=rhs.__last_stop_cond;
        __stop_cond_history_end=rhs.__stop_cond_history_end;
        __stop_cond_discarded_begin=rhs.__stop_cond_discarded_begin;
        __stop_cond_i=rhs.__stop_cond_i;
        __first_stop_deriv=rhs.__first_stop_deriv;
        __last_stop_deriv=rhs.__last_stop_deriv;
        __stop_deriv_history_end=rhs.__stop_deriv_history_end;
        __stop_deriv_discarded_begin=rhs.__stop_deriv_discarded_begin;
        __stop_deriv_i=rhs.__stop_deriv_i;
        return *this;
    }

    bool StopHistoryInterval::operator==(const StopHistoryInterval &rhs)
    {
        return __num_points==rhs.__num_points &&
            __point_i==rhs.__point_i &&
            __first_age==rhs.__first_age &&
            __last_age==rhs.__last_age &&
            __history_age_end==rhs.__history_age_end &&
            __discarded_age_begin==rhs.__discarded_age_begin &&
            __age_i==rhs.__age_i &&
            __first_stop_cond==rhs.__first_stop_cond &&
            __last_stop_cond==rhs.__last_stop_cond &&
            __stop_cond_history_end==rhs.__stop_cond_history_end &&
            __stop_cond_discarded_begin==rhs.__stop_cond_discarded_begin &&
            __stop_cond_i==rhs.__stop_cond_i &&
            __first_stop_deriv==rhs.__first_stop_deriv &&
            __last_stop_deriv==rhs.__last_stop_deriv &&
            __stop_deriv_history_end==rhs.__stop_deriv_history_end &&
            __stop_deriv_discarded_begin==rhs.__stop_deriv_discarded_begin &&
            __stop_deriv_i==rhs.__stop_deriv_i;
    }

    void StopHistoryInterval::grow_left(size_t n)
    {
        __num_points+=n;
        for(size_t i=0; i<n; i++) {
            retreat_iterator_set(__first_age, __first_stop_cond,
                    __first_stop_deriv);
            retreat_iterator_set(__age_i, __stop_cond_i, __stop_deriv_i);
        }
    }

    void StopHistoryInterval::grow_right(size_t n)
    {
        __num_points+=n;
        for(size_t i=0; i<n; i++)
            advance_iterator_set(__last_age, __last_stop_cond,__last_stop_deriv);
    }

    std::ostream &operator<<(std::ostream &os, StopHistoryInterval interval)
    {
        return os;
        std::streamsize orig_precision=os.precision();
        std::ios_base::fmtflags orig_flags=os.flags();
        os.setf(std::ios_base::scientific);
        os.precision(16);
        os << std::endl;
        size_t current_index=interval.current_point_index();
        os << std::setw(20) << "Age:";
        for(interval.reset(); !interval.end(); interval++) {
            if(current_index==interval.current_point_index())
                os << "|";
            os << std::setw(25) << interval.age();
            if(current_index==interval.current_point_index())
                os << "|";
        }
        os << std::endl;
        for(size_t cond_ind=0; cond_ind<interval.number_conditions();
                cond_ind++) {
            os << std::setw(13) << "Condition[" << std::setw(5)
                << cond_ind << "]:";
                for(interval.reset(); !interval.end(); interval++) {
                    if(current_index==interval.current_point_index())
                        os << "|";
                    os << std::setw(25)
                        << interval.stop_condition_value(cond_ind);
                    if(current_index==interval.current_point_index())
                        os << "|";
                }
            os << std::endl;
            os << std::setw(13) << "Derivative[" << std::setw(5)
                << cond_ind << "]:";
            for(interval.reset(); !interval.end(); interval++) {
                if(current_index==interval.current_point_index())
                    os << "|";
                os << std::setw(25) << interval.stop_condition_deriv(cond_ind);
                if(current_index==interval.current_point_index())
                    os << "|";
            }
            os << std::endl;
        }
        os.precision(orig_precision);
        os.flags(orig_flags);
        return os;
    }

} //End Evolve namespace.
