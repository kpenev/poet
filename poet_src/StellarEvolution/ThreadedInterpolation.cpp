/**\file
 *
 * \brief The implementation of some of the methods of the InterpolationQueue
 * class.
 *
 * \ingroup StellarEvolution_group
 */

#include "ThreadedInterpolation.h"

namespace StellarEvolution {
    void *do_interpolation (void *queue)
    {
        reinterpret_cast<InterpolationQueue*>(queue)->interpolate_thread();
        return NULL;
    }

    void InterpolationQueue::push_back(const double *x,
                                       const double *y,
                                       size_t npoints,
                                       int nodes,
                                       double smoothing,
                                       int quantity_id,
                                       int grid_index)
    {
        __x.push_back(x);
        __y.push_back(y);
        __npoints.push_back(npoints);
        __nodes.push_back(nodes);
        __smoothing.push_back(smoothing);
        __quantity_id.push_back(quantity_id);
        __grid_index.push_back(grid_index);
    }

    void InterpolationQueue::interpolate_thread()
    {
        while(true) {
            pthread_mutex_lock(&__sync_mutex);

            if(__x.size() == 0) {
                assert(__y.size() == 0);
                assert(__npoints.size() == 0);
                assert(__nodes.size() == 0);
                assert(__smoothing.size() == 0);
                pthread_mutex_unlock(&__sync_mutex);
                pthread_exit(NULL);
                assert(false);
            }

            const double *x = __x.back(),
                         *y = __y.back();
            size_t npoints = __npoints.back();
            int nodes = __nodes.back();
            double smoothing = __smoothing.back();

            __x.pop_back();
            __y.pop_back();
            __npoints.pop_back();
            __nodes.pop_back();
            __smoothing.pop_back();
            int destination = __x.size();

#ifndef NDEBUG
            std::clog
                << "Interpolating "
                << QUANTITY_NAME[*__quantity_id_iter]
                << " using " << nodes << " nodes and "
                << smoothing << " smoothing."
                << std::endl;
            ++__quantity_id_iter;
#endif

            pthread_mutex_unlock(&__sync_mutex);
            Core::InterpolatingFunctionALGLIB *quantity = 
                new Core::InterpolatingFunctionALGLIB(x,
                                                      y,
                                                      npoints,
                                                      NULL,
                                                      smoothing,
                                                      nodes);

            pthread_mutex_lock(&__sync_mutex);
#ifndef NDEBUG
            std::clog << "Created result @: " << quantity << std::endl;
#endif
            __result[destination] = quantity;
            pthread_mutex_unlock(&__sync_mutex);
        }
    }

    void InterpolationQueue::calculate(unsigned num_threads)
    {
        pthread_attr_t thread_attributes;
        pthread_attr_init(&thread_attributes);
        pthread_attr_setdetachstate(&thread_attributes,
                                    PTHREAD_CREATE_JOINABLE);
        std::vector<pthread_t> threads(num_threads);
        __result.resize(__x.size());
#ifndef NDEBUG
        __quantity_id_iter = __quantity_id.begin();
#endif
        for(unsigned i = 0; i < num_threads; ++i)
            pthread_create(&threads[i],
                           &thread_attributes,
                           do_interpolation,
                           reinterpret_cast<void*>(this));
        for(unsigned i = 0; i < num_threads; ++i)
            pthread_join(threads[i], NULL);
    }

    void InterpolationQueue::pop_front()
    {
#ifndef NDEBUG
        assert(*this);
        std::clog << "quantity ID size: " << __quantity_id.size() << std::endl;
        std::clog << "grid index size: " << __grid_index.size() << std::endl;
        std::clog << "result size: " << __result.size() << std::endl;
        std::clog << "Popping result @: " << __result[__quantity_id.size()-1] << std::endl;
#endif
        __quantity_id.pop_front();
        __grid_index.pop_front();
    }
}//End StellarEvolution namespace.
