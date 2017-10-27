/**\file
 *
 * \brief Declaration of a class that handles multithreaded stellar evolution
 * interpolation.
 *
 * \ingroup StellarEvolution_group
 */

#include "InterpolationQuantities.h"
#include "../Core/InterpolatingFunctionALGLIB.h"

#ifndef WINDOWS
    #include <pthread.h>
#endif

#include <vector>
#include <list>

namespace StellarEvolution {

    ///\brief A class that handles a queue of interpolation tasks. Also
    ///functions as an iterator over the results.
    class InterpolationQueue {
    private:
        ///The array of x values to use in the interpolation.
        std::list<const double *> __x;

        ///The array of y values to use in the interpolation.
        std::list<const double *> __y;

        ///The number of points to use in the interpolation.
        std::list<size_t> __npoints;

        ///The number of nodes to use in the interpolation.
        std::list<int> __nodes;

        ///The smoothing to use in the interpolation.
        std::list<double> __smoothing;
        
        ///The quantity ID of the quantity being interpolated.
        std::list<int> __quantity_id;

        ///The index of the grid point being interpolated.
        std::list<int> __grid_index;

        ///The interpolation results.
        std::vector<Core::InterpolatingFunctionALGLIB*> __result;

#ifndef WINDOWS
        ///\brief A pthread mutex used to ensure that only one thread is
        ///extracting the next quantity for interpolation.
        pthread_mutex_t __sync_mutex;
#endif

#ifndef NDEBUG
        std::list<int>::const_iterator __quantity_id_iter;
#endif

        ///\brief Interpolate the first quantity and discard it from __x,
        ///__y, __npoints, __nodes and __smoothing.
        void interpolate_thread();

    public:
        ///Create an empty queue.
        InterpolationQueue()
        {
#ifndef WINDOWS
            pthread_mutex_init(&__sync_mutex, NULL);
#endif
        }

        ///Add an interpolation taks to the queue.
        void push_back(const double *x,
                       const double *y,
                       size_t npoints,
                       int nodes,
                       double smoothing,
                       int quantity_id,
                       int grid_index);

        ///Carry out the interpoltaions using the specified number of
        ///simultaneous threads.
        void calculate(unsigned num_threads);

        ///The quantity ID of the currently selected result.
        int quantity_id() const {return __quantity_id.front();}

        ///The grid index of the currently selected result.
        int grid_index() const {return __grid_index.front();}

        ///The currently selected interupolation result.
        Core::InterpolatingFunctionALGLIB* result() const
        {return __result[__result.size() - __quantity_id.size()];}

        ///\brief Move to the next result in the queue (earlier results are
        ///no longer accessible.
        ///
        ///Does not actually delete the result, since it is intended to be
        ///grabbed by the user.
        void pop_front();

        ///Are there any un-popped results?
        operator bool() const {return __quantity_id.size();}

        friend void *do_interpolation(void *);

    };//End InterpolationQueue class.

} //End StellarEvolution namespace.
