#ifndef __INTERP_SOLUTION_ITERATOR_H
#define __INTERP_SOLUTION_ITERATOR_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SerializableSpline1dInterpolant.h"

#include "Common.h"
#include <valarray>

namespace Core {
    ///\brief An iterator over a set of solutions to an interpolating
    ///function.
    ///
    ///Iterates over all abscissas where the intepolating function takes a
    ///pre-determined value, merging solutions that are too close together.
    ///
    ///Works by adding solutions between later and later pairs of consecutive
    ///interpolation nodes whenever it runs out of found solutions and the
    ///next solution is requested.
    class LIB_LOCAL InterpSolutionIterator : 
        public std::iterator<std::input_iterator_tag, double, ptrdiff_t, 
        const double*, const double&>
    {
    private:
        ///The ALGLIB spline.
        const alglib_impl::spline1dinterpolant *spline;

        ///The node up to which solutions have been reported
        int node_index;

        ///\brief Controls when solutions are considered distinct.
        ///
        ///Solutions that differ by less than min_diff fraction of the
        ///distance between two nodes are treated as a single solution.
        double min_diff,

               ///Iterate over abscissas when the interpolation=y.
               y;

        ///An iterator over the list of solutions found so far. 
        std::list<double>::const_iterator solution_iter;

        ///The list of solutions found so far.
        std::list<double> solutions;

        ///\brief Whether we have gone past the last solution or before
        ///the first.
        bool is_out_of_range;

        ///\brief Find another set of solutions.
        ///
        ///Adds the solutions between the the present node and the next
        ///(or the closest node after it if there are none) to the end of
        ///solutions, incrementing the node to one past the node from
        ///which solutions were added.
        void get_solutions();
    public:
        ///Default constructor of a non meaningful object
        InterpSolutionIterator() {};

        ///Copy constructor
        InterpSolutionIterator(const InterpSolutionIterator &rhs);

        ///Start iterating over the solutions of the given spline.
        InterpSolutionIterator(
            const alglib::spline1dinterpolant &spline_var,
            double offset, double min_sol_distance=1e-8
        );

        ///Go to the next solution.
        InterpSolutionIterator &operator++();

        ///Go to the next solution.
        InterpSolutionIterator operator++(int);

        ///Go to the previous solution.
        InterpSolutionIterator &operator--();

        ///Go to the previous solution.
        InterpSolutionIterator operator--(int);

        ///Returns the current solution
        const double &operator*() const;

        ///\brief Checks if this iterator is at the same solution of the 
        ///same spline as rhs.
        bool operator==(const InterpSolutionIterator &rhs) const;

        ///\brief The opposite of operator==.
        bool operator!=(const InterpSolutionIterator &rhs) const;

        ///\brief Whether we have gone past the last solution or before
        ///the first.
        bool out_of_range() const; 
    }; //End InterpSolutionIterator class.

} //End Core namespace.

#endif


