/**\file
 *
 * \brief Declares the StopHistoryInterval class.
 *
 * \ingroup OrbitSolver_group
 */

#ifndef __STOP_HISTORY_INTERVAL_H
#define __STOP_HISTORY_INTERVAL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Error.h"
#include <valarray>
#include <list>
#include <ostream>
#include <cassert>
#include <iomanip>

namespace Evolve {

    ///\brief A collection of accepted and discarded evolution steps which
    ///contain some reason to stop.
    ///
    ///Generally tracks iterators over the stored evolution history and
    ///discarded values during the evolution calculations. One particular
    ///wrinkle is that this class deals with the possibility that step sizes can
    ///be driven to very small values by some other precision requirements
    ///leading to stopping condition values that differ only by numerical errors
    ///in their calculations causing real extrema to be missed or creating fake
    ///ones. To handle this, values that differ by less than some precsibed
    ///amount (see numerical_roundoff_factor argument to the constructor) are
    ///treated as a single value equal to their average.
    ///
    ///\ingroup OrbitSolver_group
    class LIB_LOCAL StopHistoryInterval {
    private:
        size_t __num_points,///< Number of points in the interval
               __point_i;///< The index of the current point

        ///The first age in the interval.
        std::list<double>::const_iterator __first_age,

            ///The last age in the interval.
            __last_age,

            ///\brief The one past last element of the history of stoppnig
            ///condition ages.
            __history_age_end,

            ///The first age in the discarded stopping conditions.
            __discarded_age_begin,

            ///The age of the current point.
            __age_i;

        ///The first stopping condition value in the interval
        std::list< std::valarray<double> >::const_iterator __first_stop_cond,

            ///The last stopping condition value in the interval
            __last_stop_cond,

            ///The one past last element of the history of stoppnig conditions
            __stop_cond_history_end,

            ///The first of the discarded stopping condition values.
            __stop_cond_discarded_begin,

            ///The current stopping condition value.
            __stop_cond_i,

            ///The first stopping condition derivative in the interval
            __first_stop_deriv,

            ///The last stopping condition derivative in the interval
            __last_stop_deriv,

            ///The one past last element of the history of stoppnig conditions
            __stop_deriv_history_end,

            ///The first of the discarded stopping condition values.
            __stop_deriv_discarded_begin,

            ///The current stop condition derivative
            __stop_deriv_i;

        ///\brief Increments all the iterators passed as arguments.
        ///
        ///Takes care of the switch from history to discarded if necessary.
        void advance_iterator_set(std::list<double>::const_iterator &age_i,
                std::list< std::valarray<double> >::const_iterator &cond_i,
                std::list< std::valarray<double> >::const_iterator &deriv_i);

        ///\brief Decrements all the iterators passed as arguments.
        ///
        ///Takes care of the switch from discarded to history if necessary.
        void retreat_iterator_set(std::list<double>::const_iterator &age_i,
                std::list< std::valarray<double> >::const_iterator &cond_i,
                std::list< std::valarray<double> >::const_iterator &deriv_i);

    public:
        ///Construct an interval of steps with some reason to stop the evolution.
        StopHistoryInterval(
            ///The number of points in the interval.
            size_t num_points=0,

            ///An iterator pointing to the first age in the interval.
            std::list<double>::const_iterator
            first_age=std::list<double>::const_iterator(),

            ///An iterator pointing to one past the last age in the history.
            std::list<double>::const_iterator
            history_age_end=std::list<double>::const_iterator(),

            ///An iterator pointing to the first age in the discarded list.
            std::list<double>::const_iterator
            discarded_age_begin=std::list<double>::const_iterator(),

            ///An iterator pointing to the first stopping condition in the
            ///interval.
            std::list< std::valarray<double> >::const_iterator
            first_stop_cond
            =
            std::list< std::valarray<double> >::const_iterator(),

            ///An iterator pointing to one past the last stopping condition
            ///in the history.
            std::list< std::valarray<double> >::const_iterator
            stop_cond_history_end
            =
            std::list< std::valarray<double> >::const_iterator(),

            ///An iterator pointing to the first stopping condition in the
            ///discarded list.
            std::list< std::valarray<double> >::const_iterator
            stop_cond_discarded_begin
            =
            std::list< std::valarray<double> >::const_iterator(),

            ///An iterator pointing to the first stopping derivative in the
            ///interval.
            std::list< std::valarray<double> >::const_iterator
            first_stop_deriv
            =
            std::list< std::valarray<double> >::const_iterator(),

            ///An iterator pointing to one past the last stopping derivative
            ///in the history.
            std::list< std::valarray<double> >::const_iterator
            stop_deriv_history_end
            =
            std::list< std::valarray<double> >::const_iterator(),

            ///An iterator pointing to the first stopping derivative in the
            ///discarded list.
            std::list< std::valarray<double> >::const_iterator
            stop_deriv_discarded_begin
            =
            std::list< std::valarray<double> >::const_iterator()
        );

        ///Copy orig to *this.
        StopHistoryInterval(const StopHistoryInterval &orig);

        ///Makes the current point the first point in the interval
        void reset();

        ///Advances to the next point in the interval.
        StopHistoryInterval &operator++();

        ///Advances to the next point in the interval.
        StopHistoryInterval operator++(int);

        ///Advances to the next point in the interval.
        StopHistoryInterval &operator--();

        ///Advances to the next point in the interval.
        StopHistoryInterval operator--(int);

        ///\brief Moves the entire interval, along with the current point left n
        ///points.
        ///
        ///Gains n new points at the front and loses n at the back.
        ///
        ///If there are not enough points in the history undefined behavior
        ///results.
        StopHistoryInterval &operator<<(size_t n);

        ///\brief Moves the entire interval, along with the current point right n
        ///points.
        ///
        ///Gains n new points at the back and loses n at the front.
        ///
        ///If there are not enough points in the discarded list undefined
        ///behavior results.
        StopHistoryInterval &operator>>(size_t n);

        ///Copies rhs to this.
        StopHistoryInterval &operator=(const StopHistoryInterval &rhs);

        ///Checks if the RHS is the same interval and is at the same point in it.
        bool operator==(const StopHistoryInterval &rhs);

        ///Adds the n points before the first point to the interval.
        void grow_left(size_t n=1);

        ///Adds the n points before the first point to the interval.
        void grow_right(size_t n=1);

        ///Returns the number of points in the interval.
        size_t num_points() {return __num_points;}

        ///Returns the index of the current point within the interval.
        size_t current_point_index() {return __point_i;}

        ///Returns the number of conditions at the first point.
        size_t number_conditions() {return __first_stop_cond->size();}

        ///\brief Returns true iff this is the invalid point marking the end of
        ///the interval.
        bool end() {return __point_i==__num_points;}

        ///Returns the age of the first point in the interval.
        double first_age() const {return *__first_age;}

        ///Returns the age of the last point in the interval.
        double last_age() const {return *__last_age;}

        ///Returns the age of the current point.
        double age() const {return *__age_i;}

        ///\brief Returns the value of the stop condition with the given index
        ///for the first point in the interval.
        double first_stop_condition_value(size_t condition_index) const
        {return (*__first_stop_cond)[condition_index];}

        ///\brief Returns the value of the stop condition with the given index
        ///for the last point in the interval.
        double last_stop_condition_value(size_t condition_index) const
        {return (*__last_stop_cond)[condition_index];}

        ///\brief Returns the value of the stop condition with the given index
        ///for the current point.
        double stop_condition_value(size_t condition_index) const
        {return (*__stop_cond_i)[condition_index];}

        ///\brief Returns the derivative of the stop condition with the given
        ///index for the first point in the interval.
        double first_stop_condition_deriv(size_t condition_index) const
        {return (*__first_stop_deriv)[condition_index];}

        ///\brief Returns the derivative of the stop condition with the given
        ///index for the last point in the interval.
        double last_stop_condition_deriv(size_t condition_index) const
        {return (*__last_stop_deriv)[condition_index];}

        ///\brief Returns the derivative of the stop condition with the given
        ///index for the current point.
        double stop_condition_deriv(size_t condition_index) const
        {return (*__stop_deriv_i)[condition_index];}
    }; //End StopHistoryInternval class.

    ///Civilized output of a StopHistoryInterval.
    LIB_LOCAL std::ostream &operator<<(std::ostream &os,
                                       StopHistoryInterval interval);

}//End Evolve namespace.

#endif
