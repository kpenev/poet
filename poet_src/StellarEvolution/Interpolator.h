/**\file
 *
 * \brief Defines the StellarEvolution class needed for interpolating among
 * stellar evolution tracks.
 * 
 * \ingroup StellarEvolution_group
 */

#ifndef __INTERPOLATOR_H
#define __INTERPOLATOR_H

#include "EvolvingStellarQuantity.h"
#include "InterpolationQuantities.h"
#include "ThreadedInterpolation.h"
#include "../Core/StellarZone.h"
#include "../Core/Error.h"
#include <valarray>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/valarray.hpp>
#include <boost/serialization/vector.hpp>

namespace StellarEvolution {

    ///The primordial Helium fraction of the universe.
    const double Yprimordial = 0.249;

    ///The Helium fraction with which the Sun formed.
    const double Yprotosun = 0.2612;

    ///The metal fraction with which the Sun formed.
    const double Zprotosun = 0.0150;

    ///The hydrogen fraction with which the Sun formed.
    const double Xprotosun = 1.0 - Yprotosun - Zprotosun;

    ///\brief A class that interpolates among stellar evolution tracks.
    ///
    ///Uses a  set of pre-computed evolution tracks to generate 
    ///inrpolating functions that represents reasonably well the evolution of
    ///various properties of an arbitrary mass star as a function of age.
    ///
    ///At present the only implementing class is MESA::Interpolator based on
    ///a set of MESA tracks.
    ///
    ///\ingroup StellarSystem_group
    class Interpolator {
    private:
        friend class boost::serialization::access;

        ///Serialize the current interpolation.
        template<class Archive> void serialize(
            ///The archive to serialize to.
            Archive & ar, 

            ///Version number. Ignored!
            const unsigned int
        ); 

        std::valarray<double> 
            ///The stellar masses for which evolution tracks are available in
            /// \f$M_\odot\f$
            __track_masses,

            ///\brief The stellar metallicities for which evolution tracks
            ///are available in [Fe/H]
            __track_metallicities;

        ///\brief The interpolated stellar evolution quantities for each
        ///track.
        ///
        ///See ::QuantityID for the order.
        std::vector< std::vector<const OneArgumentDiffFunction*> >
            __interpolated_quantities;

        ///Was the interpolation of the corresponding quantity vs. log(age)?
        std::vector<bool> __vs_log_age;

        ///Was the interpolation of the log(corresponding quantity)?
        std::vector<bool> __log_quantity;

        ///The age at which the core starts forming in Gyr.
        double __core_formation;

        ///Return the index of the first non-zero value in the argument.
        int find_first_core_index(
            const std::valarray<double> &core_mass
        ) const;

        ///Perform all queued interpolations.
        void perform_queued_interpolations(
            ///The queue of pending interpolations to calculate.
            InterpolationQueue &interpolation_queue,

            ///The number of threads to use for simultaneous processing
            unsigned num_threads = 1
        );
    public:
        ///\brief Construct an object that can be set to interpolate between
        ///tabulated evolution tracks.
        Interpolator() :
            __track_masses(),
            __track_metallicities(),
            __interpolated_quantities(NUM_QUANTITIES),
            __core_formation(Core::NaN)
        {}

        ///\brief Creates a fully functional stellar evolution interpolator.
        Interpolator(
            ///The stellar masses (in \f$M_\odot\f$) for which evolution
            ///tracks are tabulated.
            const std::valarray<double> &tabulated_masses,

            ///The stellar metallicities (in [Fe/H]) for which evolution
            ///tracks are tabulated.
            const std::valarray<double> &tabulated_metallicities,

            ///A set of ages for each track in Gyr on the grid defined by 
            ///\p track_masses and \p track_metallicities. The mass index
            ///varies faster.
            const std::list< std::valarray<double> > &tabulated_ages,

            ///A set of stellar quantities for each age of each track. See
            ///StellarEvolution::QuantityID for the list, order and units
            ///of the quantities.
            const std::vector< std::list< std::valarray<double> > > 
            &tabulated_quantities,

            ///\brief How much to smooth each quantity when fitting. Use NaN
            ///for no smoothing. Corresponds entry by entry with
            /// \p tabulated_quantities.
            const std::vector<double> &smoothing,

            ///How many nodes to use when smoothing each quantity (ignored if
            ///the corresponding smoothing entry is NaN - no smoothing).
            ///Corresponds entry by entry with \p tabulated_quantities and
            ///\p smoothing.
            const std::vector<int> &nodes,

            ///Should interpolation be done vs. log(age) instead of age for
            ///each quantity?
            const std::vector<bool> &vs_log_age,

            ///Should interpolation be done of log(quantity) instead of
            ///quantity for each quantity?
            const std::vector<bool> &log_quantity,

            ///The number of simultaneosly running threads to use for the
            ///interpolation.
            unsigned num_threads
        )
        {
            create_from(tabulated_masses,
                        tabulated_metallicities,
                        tabulated_ages,
                        tabulated_quantities,
                        smoothing,
                        nodes,
                        vs_log_age,
                        log_quantity,
                        num_threads);
        }

        ///Fully setup an object created by the default constructor.
        void create_from(
            ///The stellar masses (in \f$M_\odot\f$) for which evolution
            ///tracks are tabulated.
            const std::valarray<double> &tabulated_masses,

            ///The stellar metallicities (in [Fe/H]) for which evolution
            ///tracks are tabulated.
            const std::valarray<double> &tabulated_metallicities,

            ///A set of ages for each track in Gyr on the grid defined by 
            ///\p track_masses and \p track_metallicities. The mass index
            ///varies faster.
            const std::list< std::valarray<double> > &tabulated_ages,

            ///A set of stellar quantities for each age of each track. See
            ///StellarEvolution::QuantityID for the list, order and units
            ///of the quantities.
            const std::vector< std::list< std::valarray<double> > > 
            &tabulated_quantities,

            ///\brief How much to smooth each quantity when fitting. Use NaN
            ///for no smoothing. Corresponds entry by entry with
            /// \p tabulated_quantities.
            const std::vector<double> &smoothing,

            ///How many nodes to use when smoothing each quantity (ignored if
            ///the corresponding smoothing entry is NaN - no smoothing).
            ///Corresponds entry by entry with \p tabulated_quantities and
            ///\p smoothing.
            const std::vector<int> &nodes,

            ///Should interpolation be done vs. log(age) instead of age for
            ///each quantity?
            const std::vector<bool> &vs_log_age,

            ///Should interpolation be done of log(quantity) instead of
            ///quantity for each quantity?
            const std::vector<bool> &log_quantity,

            ///The number of simultaneosly running threads to use for the
            ///interpolation.
            unsigned num_threads
        );

        ///\brief Return a single quantity interpolation to a given mass and
        ///[Fe/H].
        ///
        ///The result must be destroyed when it becomes obsolete.
        EvolvingStellarQuantity *operator()(
            ///The quantity for which to set-up the interpolation.
            QuantityID quantity,

            ///The stellar mass to which to interpolate in \f$M_\odot\f$.
            double mass,

            ///The stellar [Fe/H] to which to interpolate.
            double feh
        ) const;

        ///The age at which the core begins to form in Gyr.
        virtual double core_formation_age() const {return __core_formation;}

        ///\brief Serializes the interpolation state to file.
        ///
        ///Only call this on objects initialized with the default
        ///constructor.
        virtual void save_state(
            ///The name of the file to save the state to.
            const std::string &filename="../interp_state_data"
        ) const;

        ///\brief Loads data from serialization.
        ///
        ///Only call this on objects NOT initialized using the default
        ///constructor (otherwise it has no data to save). Serializes state
        ///to file.
        ///
        ///Recursively saves data of YRECEvolution and every class it depends
        ///on:
        /// - StellarEvolution
        /// - InterpolatingFunctionALGLIB
        /// - OneArgumentDiffFunction,
        /// - OneArgumentFunction
        /// - spline1dinterpolant
        /// - _spline1dinterpolant_owner.
        ///
        ///In _spline1dinterpolant_owner, serialize() serializes everything
        ///EXCEPT p_struct->x.data and p_struct->y.data, because those are
        ///just copies of the original data on which the spline was based and
        ///are not necessary for evaluating the spline.
        virtual void load_state(
            ///The name of a file previously created using save_state()
            const std::string &filename="../interp_state_data"
        );

        ///\brief Free all evolution tracks, rendering all created quantities
        ///unuseable!
        ///
        ///Essentially an only explicitly called destructor. The reason for
        ///requiring explicit destruction is that generated quantities cannot
        ///be used after this method is invoked, and not invoking it simply
        ///causes a benign memory leak.
        void delete_tracks();

        virtual ~Interpolator() {}
    }; //End of Interpolator class declaration.

    ///Serialize the current interpolation.
    template<class Archive> void Interpolator::serialize(Archive & ar,
                                                         const unsigned int)
    {
        if(!std::isfinite(__core_formation)) __core_formation = -1;
        ar & __track_masses;
        ar & __track_metallicities;

        ar & __interpolated_quantities;

        ar & __vs_log_age;
        ar & __log_quantity;

        ar & __core_formation;
        if(__core_formation < 0) __core_formation = Core::Inf;
    }

    ///\brief Return the metallicity interpolation parameter corresponding to
    ///the given [Fe/H] value.
    double metallicity_from_feh(double feh);

    ///\brief Return the [Fe/H] value corresponding to the given metallicity
    //interpolation parameter.
    double feh_from_metallicity(double metallicity);

}//End of StellarEvolution namespace

#endif
