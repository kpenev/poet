/**\file
 *
 * \brief Declares the StopInformation class.
 *
 * \ingroup OrbitSolver_group
 */

#ifndef __STOP_INFORMATION_H
#define __STOP_INFORMATION_H

#include "../Core/SharedLibraryExportMacros.h"
#include "StoppingCondition.h"

namespace Evolve {

    ///\brief The information about why and where the evolution should stop.
    ///
    ///\ingroup OrbitSolver_group
    class LIB_LOCAL StopInformation {
    private:
        ///The target stopping age in Gyr.
        double __stop_age,

               ///\brief The precision up to which the reason to stop is satisfied.
               ///
               ///If we are stopping because some stopping condition is actually
               ///satisfied, this is the value of the condition. If we are
               ///stopping because of an extremum this is the difference between
               ///the estimated value of the condition at the extremum and the
               ///value at the closest step.
               __stop_condition_precision;

        ///The reason for stopping.
        StoppingConditionType __stop_reason;

        bool 
            ///Is the reason for stopping that a condition actually crossed zero?
            __is_crossing, 

            ///Did the last step actually cross zero?
             __crossed_zero;

        ///The index of the condition which caused us to stop.
        size_t __stop_condition_index;

        ///The sign of the derivative at zero-crossing.
        short __deriv_sign_at_crossing;
    public:
        ///\brief Create an object with the information about why we evolution
        ///should be stopped.
        StopInformation(
                ///The target stopping age in Gyr.
                double stop_age=Core::Inf,
                
                ///The precision up to which the reason to stop is satisfied.
                double stop_precision=Core::NaN,

                ///The reason for stopping.
                StoppingConditionType stop_reason=NO_STOP,

                ///Is the reason for stopping that a condition actually crossed
                ///zero?
                bool is_crossing = false,

                ///Did we stop after a zero-crossing.
                bool crossed_zero = false,
                
                ///The index of the condition which caused us to stop.
                size_t stop_condition_index = 0,
                
                ///The sign of the derivative of the condition at zero-crossing
                ///(undefined for extrema.
                short deriv_sign_at_crossing = 0
        ) :
            __stop_age(stop_age),
            __stop_condition_precision(stop_precision),
            __stop_reason(stop_reason),
            __is_crossing(is_crossing),
            __crossed_zero(crossed_zero),
            __stop_condition_index(stop_condition_index),
            __deriv_sign_at_crossing(is_crossing ? deriv_sign_at_crossing : 0) {}

        ///Copy orig to *this.
        StopInformation(const StopInformation &orig) :
            __stop_age(orig.__stop_age),
            __stop_condition_precision(orig.__stop_condition_precision),
            __stop_reason(orig.__stop_reason),
            __is_crossing(orig.__is_crossing),
            __crossed_zero(orig.__crossed_zero),
            __stop_condition_index(orig.__stop_condition_index),
            __deriv_sign_at_crossing(orig.__deriv_sign_at_crossing)
        {}

        ///The target stopping age in Gyr.
        double stop_age() const {return __stop_age;}

        ///The target stopping age in Gyr.
        double &stop_age() {return __stop_age;}

        ///\brief The precision up to which the reason to stop is satisfied.
        double stop_condition_precision() const
        {return __stop_condition_precision;}

        ///\brief The precision up to which the reason to stop is satisfied.
        double &stop_condition_precision()
        {return __stop_condition_precision;}

        ///The reason for stopping.
        StoppingConditionType stop_reason() const {return __stop_reason;}

        ///The reason for stopping.
        StoppingConditionType &stop_reason() {return __stop_reason;}

        ///Is the reason for stopping that a condition actually crossed zero?
        bool is_crossing() const {return __is_crossing;}
        ///Is the reason for stopping that a condition actually crossed zero?
        bool &is_crossing() {return __is_crossing;}

        ///The index of the condition which caused us to stop.
        size_t stop_condition_index() const {return __stop_condition_index;}

        ///The index of the condition which caused us to stop.
        size_t &stop_condition_index() {return __stop_condition_index;}

        ///\brief Did we stop after the stopping condition crossed zero (always
        ///false for extrema).
        bool crossed_zero() const {return __crossed_zero;}

        ///\brief Did we stop after the stopping condition crossed zero (always
        ///false for extrema).
        bool &crossed_zero() {return __crossed_zero;}

        ///Copy rhs to *this.
        StopInformation &operator=(const StopInformation &rhs)
        {
            __stop_age=rhs.__stop_age;
            __stop_condition_precision=rhs.__stop_condition_precision;
            __stop_reason=rhs.__stop_reason;
            __is_crossing=rhs.__is_crossing;
            __crossed_zero=rhs.__crossed_zero;
            __stop_condition_index=rhs.__stop_condition_index;
            __deriv_sign_at_crossing=rhs.__deriv_sign_at_crossing;
            return *this;
        }

        ///The sign of the derivative at zero-crossing.
        short deriv_sign_at_crossing() const
        {
            assert(__is_crossing);

            return __deriv_sign_at_crossing;
        }

        ///The sign of the derivative at zero-crossing.
        short &deriv_sign_at_crossing() {return __deriv_sign_at_crossing;}

    };//End StopInformation class.

    ///Civilized output of a StopInformation object.
    LIB_LOCAL std::ostream &operator<<(std::ostream &os,
                                       const StopInformation &stop);

}//End Evolve namespaec.
#endif
