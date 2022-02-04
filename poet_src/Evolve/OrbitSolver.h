/**\file
 *
 * \brief Defines the OrbitSolver class, the various stopping conditions and
 * a number of other classes used while calculating the orbital evolution.
 *
 * \ingroup OrbitSolver_group
 */

#ifndef __ORBIT_SOLVER_H
#define __ORBIT_SOLVER_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/AstronomicalConstants.h"
#include "../Core/Common.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/Common.h"
#include "BinarySystem.h"
#include "CombinedStoppingCondition.h"
#include "ExternalStoppingConditions.h"
#include "StopInformation.h"
#include "StopHistoryInterval.h"
#include <math.h>
#include <list>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <sstream>
#include <iostream>
#include <limits>
#include <ctime>

namespace Evolve {

    typedef int (*GSL_ODE_TYPE)(double, const double*, double*, void*);
    typedef int (*GSL_JAC_TYPE)(double, const double*, double*, double*, void*);
    typedef bool (*STOP_EVOL_TYPE)(double, const double*, void*);

    ///\brief A wrapper tha allows the stellar system differential equation to be
    ///passed to the GSL ODE solver.
    LIB_LOCAL int stellar_system_diff_eq(
            ///System age in Gyr.
            double age,

            ///A C-style array of the variables evolved under the current set of
            ///equations.
            const double *parameters,

            ///A C-style array of the rates of change of the orbital_parameters.
            double *derivatives,

            ///A C-style array of pointers giving the stellar system, the
            ///current evolution mode, the wind saturation state to assume,
            ///the spin-orbit lock of the star and the dissipation rates.
            void *system);

    ///\brief A wrapper tha allows the stellar system jacobian to be passed
    ///to the GSL ODE solver.
    LIB_LOCAL int stellar_system_jacobian(
            ///System age in Gyr.
            double age,

            ///A C-style array of the variables evolved under the current set of
            ///equations.
            const double *parameters,

            ///A C-style array of the Jacobian matrix entries.
            double *param_derivs,

            ///A C-style array of the partial age derivatives of the differential
            ///equations. Age is in Gyr.
            double *age_derivs,

            ///A C-style array of pointers giving the stellar system and the
            ///current evolution mode.
            void *system_mode);

    ///\brief Infomation about an extremum of a function.
    ///
    ///\ingroup OrbitSolver_group
    class LIB_LOCAL ExtremumInformation {
    private:
        ///The value of the argument where the extremum occurs.
        double __x,

               ///The value of the function at the extremum.
               __y;
    public:
        ///Create an object to hold information about a function extremum.
        ExtremumInformation(
                ///The value of the argument where the extremum occurs.
                double x=Core::Inf,

                ///The value of the function at the extremum.
                double y=Core::NaN) : __x(x), __y(y) {}

        ///The value of the argument where the extremum occurs.
        double x() const {return __x;}
        ///The value of the argument where the extremum occurs.
        double &x() {return __x;}

        ///The value of the function at the extremum.
        double y() const {return __y;}
        ///The value of the function at the extremum.
        double &y() {return __y;}
    };//End ExtremumInformation class.

    ///\brief Solves the system of ODEs describing the evolution of a
    ///single planet around a single star.
    ///
    ///\ingroup OrbitSolver_group
    class LIB_PUBLIC OrbitSolver {
    private:
        double
            ///The last  age for which evolution is required.
            __end_age,

            ///The precision required of the solution
            __precision,

            ///The last age at which the eccentricity order was increased
            __last_order_upgrade_age,

            ///Extremum searching cannot limit the step size below this
            __min_extremum_search_step;

        ///The ages at which solution is tabulated
        std::list<double> __tabulated_ages;

        ///The evolution mode corresponding to the matching tabulated age.
        std::list<Core::EvolModeType> __tabulated_evolution_modes;

        ///\brief The number of points at the start of the history to skip when
        ///lookng for a zero crossing for each condition
        std::valarray<size_t> __skip_history_zerocrossing;

        ///The age after which to look for extrema for each condition
        std::valarray<double> __skip_history_extremum;

        ///The ages at which the stop condition history is kept
        std::list<double> __stop_history_ages,

            ///\brief The ages of steps which were discarded becauset they are
            ///past a zero or an extremum of a stopping condition .
            __discarded_stop_ages;

        std::list< std::valarray<double> >
            __orbit_history,///< Past orbits
            __orbit_deriv_history, ///< Past orbital derivatives
            __stop_cond_history,///< Past values of the stop conditions

            ///Past values of the stop condition derivatives
            __stop_deriv_history,

            ///\brief Discarded values of the stop conditions.
            ///
            ///Useful for interpolating to zeroes and extrema.
            __stop_cond_discarded,

            ///\brief Discarded derivatives of the stop conditions.
            ///
            ///Useful for interpolating to zeroes and extrema.
            __stop_deriv_discarded;

        ///The current set of stopping conditions.
        StoppingCondition *__stopping_conditions;

        ///The full information about the next age where each stopping
        ///condition indicates evolution should be stopped and why.
        std::vector<StopInformation> __stop_info;

        ///See print_progress argument of constructor.
        bool __print_progress;

        ///When did the currently running evolution start.
        time_t __evolution_start_time;

        ///Max number of seconds current evolution is allowed to run.
        double __runtime_limit;

        ///Max number of steps allowed to be stored in history and/or discarded.
        unsigned __num_step_limit;

#ifndef NDEBUG
        ///\brief Generates a nicely formatted table of the contents of the
        ///discarded and history stopping condition information.
        ///
        ///Only exists if NDEBUG is not \#defined.
        void output_history_and_discarded(std::ostream &os);
#endif

        ///Removes all stored discarded stop condition information.
        void clear_discarded();

        ///\brief Adds an entry in the discarded ages, stop conditions and
        ///derivatives.
        ///
        ///Makes sure the ages remain ordered.
        void insert_discarded(double age,
                const std::valarray<double> &current_stop_cond,
                const std::valarray<double> &current_stop_deriv);

        ///Adds the last step to the evolution.
        void add_to_evolution(
                ///The age of the evolving stellar system.
                double age,

                ///The evolution mode represented in orbit.
                Core::EvolModeType evolution_mode,

                ///The binary system being evolved.
                BinarySystem &system);

        ///\brief Rewinds the evlution to the last step before the given age and
        ///returns the age of that step.
        ///
        ///Sets the orbit to what it was at that step and removes any items from
        ///the histories and tabulations that are later than max_age.
        double go_back(double max_age,
                       BinarySystem &system,
                       std::valarray<double> &orbit);


        ///Clears the current stopping condition history.
        void clear_history();

        ///\brief Selects a history interval for interpolating to a reason to
        ///stop the evolution.
        ///
        ///Finds the smallest possible interval that contains a zero crossing/or
        ///an extremum straddling the history and discarded stop conditions,
        ///containing at most the specified number of points (could be less if
        ///there are not enough points). The interval is also guaranteed to
        ///contain at least one point in the history and one point in the
        ///discarded list.
        StopHistoryInterval select_stop_condition_interval(bool crossing,
                size_t cond_ind, size_t max_points) const;

        ///\brief Estimates the value and age of an extremum of a stopping
        ///condition for which no derivative information is available.
        ///
        ///If no extremum is indicated by the points, or if it is in the wrong
        ///direction to cause a zero-crossing, returns the result of the default
        ///constructor of ExtremumInformation.
        ExtremumInformation extremum_from_history_no_deriv(
                size_t condition_index
        ) const;

        ///\brief Estimates the value and age of an extremum of a stopping
        ///condition for which derivative information is available.
        ///
        ///If no extremum is indicated by the points, or if it is in the wrong
        ///direction to cause a zero-crossing, returns the result of the default
        ///constructor of ExtremumInformation.
        ExtremumInformation extremum_from_history(
            size_t condition_index,
            double min_extremum_x
        ) const;

        ///\brief Estimates the age at which a stopping condition with no
        ///derivative information crossed zero.
        ///
        ///If no zero-crossing is indicated, Inf is returned.
        double crossing_from_history_no_deriv(
                ///The index of the condition to search.
                size_t condition_index) const;

        ///\brief Estimates the age at which a stopping condition with
        ///derivative information crossed zero.
        ///
        ///If no zero-crossing is indicated, Inf is returned.
        double crossing_from_history(
                ///The index of the condition to search.
                size_t condition_index) const;

        ///\brief Initializes the skip_history_zerocrossing and
        ///skip_history_extremum arrays appropriately after a mode change.
        void initialize_skip_history(
            const StoppingCondition &stop_cond,
            const std::valarray<double> &stop_cond_values,
            StoppingConditionType stop_reason
        );

        ///\brief Is the condition causing a stop match to within the required
        ///precision?
        bool at_exact_condition(double previous_age,
                                const StopInformation &stop_info);

        ///\brief Updates the skip_history_zerocrossing and
        ///skip_history_extremum arrays appropriately after an acceptable step.
    	void update_skip_history(
                const std::valarray<double> &current_stop_cond,
                const StopInformation &stop_info);

        ///\brief Return true iff the step with the given stop information is
        ///acceptable.
        bool acceptable_step(double current_age,
                             double previous_age,
                             const StopInformation &stop_info);

        ///\brief Updates stop_cond_history and stop_deriv_history after a
        ///GSL step, returning if/where the evolution needs to stop.
        ///
        ///Returns the index of the condition with the closest estimated age
        ///for zero crossing or an extremum exists which might have
        ///crossed zero.
        StopInformation update_stop_condition_history(
            ///System age in Gyr.
            double age,

            ///The values of the current evolution variables.
            const std::valarray<double> &orbit,

            ///The rates of change of the evolution variables per Gyr.
            const std::valarray<double> &derivatives,

            ///The current evolution mode.
            Core::EvolModeType evolution_mode,

            ///Current eccentricity expansion order
            unsigned current_expansion_order,

            ///For the first call of this function for an evolution
            ///stretch, this should indicate the reason why the previous
            ///stretch was stopped. For subsequent calls during the same
            ///evolution stretch it should be NO_STOP.
            StoppingConditionType stop_reason=NO_STOP
        );

        ///Handle the situation when the last step has to be rejected.
        void reject_step(
            ///The age up to which the rejected step reached, on exit gets
            ///updated with the last acceptable age found.
            double &age,

            ///The reason we stopped. If a zero-crossing, the precision is
            ///updated.
            StopInformation &stop,

            ///The binary system being evolved.
            BinarySystem &system,

            ///The current orbit of the system, gets updated with the last
            ///acceptable orbit found.
            std::valarray<double> &orbit,

            ///Gets updated with the maximum age up to which evolution is
            ///allowed to proceed.
            double &max_next_t,

            ///Gets updated with the step size to use for the next step.
            double &step_size
#ifndef NDEBUG
            ///The reason for rejecting the step.
            , std::string reason
#endif
        );

        ///\brief Evolves a system until either some age cut-off is reached or
        ///some stopping condition crosses zero.
        ///
        ///
        ///Appends each accepted step to tabulated_ages,
        ///tabulated_evolution_mode, tabulated_orbit and tabulated_deriv.
        ///
        ///The return value is true if the last step finished after the
        ///stopping condition crossed zero and false if it ended before that.
        StopInformation evolve_until(
            ///The planet-star system to evolve.
            BinarySystem &system,

            ///The age at which to stop this part of the evolution. On
            ///exit, it is overwritten with the age of the last accepted
            ///step.
            double &max_age,

            ///The initial conditions. The contents depends on the value
            ///of evolution_mode. See #BinarySystem.differential
            ///equations for details.
            ///
            ///On exit, it is overwritten with the orbit of the last
            ///accepted step.
            std::valarray<double> &orbit,

            ///On input should be the reason why the last evolution
            ///stopped. It should be NO_STOP if this is the first piece
            ///of evolution being calculated. On exit it is overwritten
            ///with the value appropriate for the next run.
            StoppingConditionType &stop_reason,

            ///The maximum step that GSL is allowed to take.
            double max_step,

            ///The evolution mode for this part of the evolution.
            Core::EvolModeType evolution_mode
        );

        ///\brief Returns the stopping conditions which end the given
        ///evolution mode and update __stop_info.
        CombinedStoppingCondition *get_stopping_condition(
            ///The system being evolved.
            BinarySystem &system
        );

        ///\brief The age at which the evolution should stop next if no other
        ///stopping condition occurs.
        double stopping_age(
                ///The age from which the next part of the evolution starts.
                double age,

                ///The stellar system being evolved.
                const BinarySystem &system,

                ///A sorted list of ages which must be stopped at.
                const std::list<double> &required_ages);

        ///\brief Handle a stop in the evolution due to at least one
        ///condition reaching a critical value.
        ///
        ///Assumes that __stop_info was updated per the last accepted step.
        void reached_stopping_condition(
            ///The age of the last accepted evolution step.
            double stop_age,

            ///The type of condition which caused the stoppage.
            StoppingConditionType stop_reason
        );

        ///\brief Increase/decrease the eccentricity expansion order until error
        ///is acceptable and return the new order.
        void adjust_expansion_order(
            ///The system being evolved.
            BinarySystem &system,

            ///The current orbital elements.
            const std::valarray<double> &orbit,

            ///The evolution mode for this part of the evolution.
            Core::EvolModeType evolution_mode,

            ///If true the eccentircity order is increased by at least 1.
            bool must_increase = false
        );


        ///Clears any previously calculated evolution.
        void reset(BinarySystem &system);
    public:
        ///\brief Prepare to solve for the orbital evolution.
        OrbitSolver(
            ///The end age for the evolution.
            double max_age,

            ///The precision which to require of the solution.
            double required_precision,

            ///Should output be created to track the progress in calculating the
            ///evolution?
            bool print_progress=false
        );

        ///\brief Actually solves the given differential equation with the given
        ///boundary conditions.
        void operator()(
            ///The stellar system to calculate the evolution for
            BinarySystem &system,

            ///The maximum size of the time steps allowed (useful if finer
            ///sampling of the output than default is necessary).
            double max_step=Core::Inf,

            ///A sorted list of ages to include in the tabulated evolution.
            const std::list<double> &required_ages=std::list<double>(),

            ///The maximum amount of wall-clock time in seconds to allow for
            ///calculating the evolution. If it runs out, whatever portion is
            ///calculated is preserved.
            double max_runtime=0,

            ///The maximum number of time steps to store in history (+discarded)
            unsigned max_time_steps=0,

            ///The minimum time step to impose when searching for extrema of
            ///stopping conditions. Note that the GSL ODE solver and searching
            ///for zero crossing of stopping conditions can still impose smaller
            ///steps.
            double min_extremum_search_step=1e-5
        );

        ///The ages at which evolution has been tabulated so far.
        const std::list<double> &evolution_ages() const
        {return __tabulated_ages;}

        ///The tabulated evolution modes so far.
        const std::list<Core::EvolModeType> &mode_evolution() const
        {return __tabulated_evolution_modes;}

        ///Clean up.
        ~OrbitSolver()
        {if(__stopping_conditions) delete __stopping_conditions;}

    }; //End OrbitSolver class.

}//End Evolve namespace.

#endif
