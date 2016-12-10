/**\file
 *
 * \brief Declares a class implementing the intepolation of a single stellar
 * quantity from stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __EVOLVING_STELLAR_QUANTITY_H
#define __EVOLVING_STELLAR_QUANTITY_H

#include "AllowedGridGrowth.h"
#include "RemoveLogDeriv.h"
#include "SumDerivatives.h"
#include "../Core/Functions.h"
#include "../Core/Error.h"
#include <valarray>
#include <list>
#include <string>
#include <iostream>
#include <fstream>

namespace StellarEvolution {

    using Core::OneArgumentDiffFunction;
    using Core::FunctionDerivatives;
    using Core::InterpSolutionIterator;

    ///\brief A class for stellar properties that depend on age.
    ///
    ///\ingroup StellarSystem_group
    class EvolvingStellarQuantity : public OneArgumentDiffFunction {
    private:
        double 
            ///The mass to which to interpolate in \f$M_\odot\f$.
            __mass,

            ///The mass to which to interpolate in Solar metallicities.
            __metallicity;


        double 
            ///The minimum age for which this quantity is defined in Gyr.
            __min_age,

            ///The maximum age for which this quantity is defined in Gyr.
            __max_age; 

        ///Whether the tracks have log(age) instead of age as their argument.
        bool __use_log_age;

        ///Should the quantity be assumed zero below the minimum track age.
        bool __initially_zero;

        std::valarray<double> 
            ///\brief The masses of the evolution tracks.
            __track_masses,

            ///The metallicities of the evolution tracks.
            __track_metallicities;

        std::valarray<double> 
            ///\brief The minimum interpolation age for the current star to 
            ///which each track can contribute.
            __min_interp_ages,

            ///\brief The maximum interpolation age for the current star to 
            ///which each track can contribute.
            __max_interp_ages;

        ///\brief The ages at which the interpolation grid needs to be
        ///re-determined.
        std::vector<double> __interp_grid_change_ages;

        ///\brief The entry in ::__interp_grid_change_ages up to which the
        ///current interpolation grid is valid.
        mutable std::vector<double>::const_iterator __next_grid_change_age;

        mutable size_t 
            ///\brief The index of the smallest track mass not smaller than
            /// ::__mass.
            __mass_index_above,

            ///The index of the largest track mass not exceeding ::__mass.
            __mass_index_below,

            ///\brief The index of the smallest track metallicity not smaller
            ///than ::__metallicity.
            __metallicity_index_above,

            ///\brief The index of the largest track metallicity not
            ///exceeding ::__metallicity.
            __metallicity_index_below;

        ///\brief The model tracks for the evolution of the quantity on the
        ///grid defined by ::__track_masses and ::__track_metallicities.
        ///
        ///The mass index varies faster.
        std::vector<const OneArgumentDiffFunction *> __evolution_tracks;

        ///\brief How far to extrapolate past the last tabulated age for
        ///a model.
        ///
        ///Models are included in the mass interpolation only if
        ///the required age is no larger than the maximum tabulated age 
        ///times this factor.
        double __extrapolate;

        mutable alglib::real_1d_array
            ///The current track masses participating in the interpolation.
            __interp_masses,

            ///\brief The current track metallicities participating in the
            ///interpolation.
            __interp_metallicities;

        mutable size_t 
            ///\brief The index within ::__track_masses of the lowest mass
            ///currently participating in the interpolation.
            __min_interp_mass_index,

            ///\brief The index within ::__track_masses of the highest mass
            ///currently participating in the interpolation.
            __max_interp_mass_index,

            ///\brief The index within ::__track_metallicities of the lowest
            ///metallicity currently participating in the interpolation.
            __min_interp_metallicity_index,

            ///\brief The index within ::__track_metallicities of the highest
            ///metallicity currently participating in the interpolation.
            __max_interp_metallicity_index;

        ///\brief Return the index within ::__evolution_tracks for the given 
        ///mass and metallicity indices.
        inline size_t track_index(
            ///The index within ::__track_masses of the desired mass.
            size_t mass_index,

            ///The index within ::__track_metallicities of the desired
            ///metallicity.
            size_t metallicity_index) const
        {return metallicity_index * __track_masses.size() + mass_index;}

        ///\brief Answer if a given track can participate in interpolating to 
        ///the given age.
        inline bool track_in_range(
            ///The index of the track to check within ::__evolution_tracks.
            size_t track_i,

            ///The age to which interpolation is desired.
            double age) const
        {return (__min_interp_ages[track_i] <= age
                 &&
                 __max_interp_ages[track_i] >= age);}

        ///\brief Answer if a given track can participate in interpolating to 
        ///the given age.
        inline bool track_in_range(
            ///The index of the mass of the track to check within
            ///::__track_masses.
            size_t mass_i,

            ///The index of the metallicity of the track to check within
            ///::__track_metallicities.
            size_t metallicity_i,

            ///The age to which interpolation is desired.
            double age) const
        {return track_in_range(track_index(mass_i, metallicity_i), age);}

        ///\brief Verify that the stellar mass and metallicity are within 
        ///range of the evolution tracks.
        void check_grid_range() const;

        ///\brief The two indices within the given sorted array defining the
        ///closed internal containing value.
        ///
        ///If value is not exactly equal to an array entry, the two indices 
        ///are consecutive, if the value is exactly equal to an entry, the
        ///two indices are the same.
        void find_cell(
            ///The boundaries of the grid cells in a single dimension.
            const std::valarray<double> &boundaries,

            ///The value whose cell we are looking for.
            double value,

            ///The index of the largest boundary <= \p value.
            size_t &below_index,

            ///The index of the smallest boundary >= \p value.
            size_t &above_index
        ) const;

        ///Fill the ::__min_interp_ages and ::__max_interp_ages members.
        void set_interp_age_ranges();

        ///\brief Interpolate the quantity for the given track to the given 
        ///age, returning NaN if out of age range.
        ///
        ///If derivatives is not NULL initializes that to a pointer to a
        ///derivatives at the current age structure.
        double evaluate_track(
            ///The age at which to evaluate the track. If ::__use_log_age is
            ///true, the track is evaluated at log(\p track_age), otherwise
            ///it is directly passed to the track. In particular, the caller
            ///should have already transformed this argument to the correct
            ///interpolation parameter.
            double age,

            ///The track to evaluate. Should be an entry from
            ///::__evolution_tracks.
            const OneArgumentDiffFunction &track,

            ///If not NULL, \p *derivatives is set to a newly allocated
            ///derivative instance of one of the children FunctionDerivatives
            ///classes.
            const FunctionDerivatives **derivatives) const;

        ///\brief Figure out in which directions we can expand a
        ///mass-metallicity interpolation grid, assuming a single direction
        ///expansion.
        void check_grid_expansion_directions(
            ///The current state of the grid expansion possibilities. On
            ///output, directions in which growth is no longer allowed are
            ///disabled.
            AllowedGridGrowth &grow,

            ///The age for which we will interpolate.
            double age
        ) const;

        ///Expand the current interpolation grid in one of the allowed
        ///directions.
        void expand_grid(
            ///The current state of the grid expansion possibilities.
            const AllowedGridGrowth &grow,

            ///The age to which we will interpolate.
            double age
        ) const;

        ///Find the best sub-grid of tracks to interpolate on.
        void update_interpolation_grid() const;

        ///\brief Interpolate the quantity to the desired age.
        ///
        ///If derivatives is not NULL initializes that to a pointer to a
        ///derivatives at the current age structure.
        double interpolate(
            double age, 
            const FunctionDerivatives **derivatives=NULL
        ) const;

        ///\brief Return the interpoltaion parameter for the given age for
        ///the current star.
        template<typename VALUE_TYPE>
            inline const VALUE_TYPE &age_to_interp_param(
                const VALUE_TYPE &age
            ) const
            {return age_to_interp_param(age, __mass, __metallicity);}

        ///\brief Return the age for the given interpoltaion parameter for
        ///the current star.
        template<typename VALUE_TYPE>
            inline const VALUE_TYPE &interp_param_to_age(
                const VALUE_TYPE &interp_param
            ) const
            {
                return interp_param_to_age(interp_param,
                                           __mass,
                                           __metallicity);
            }

    protected:
        ///\brief Return the interpoltaion parameter given age, mass and
        ///metallicity.
        ///
        ///Must be an increasing monotonic function
        virtual double age_to_interp_param(
            ///The age for which the interpolation parameter is needed in
            ///Gyrs.
            double age,

            ///The stellar mass for which the interpolation parameter is
            ///needed in \f$M_\odot\f$.
            double mass,

            ///The stellar metallicity for which the interpolation
            ///parameter is needed in \f$[Fe/H]\f$.
            double metallicity) const;

        ///\brief Return the age in Gyrs given an interpolation parameter, 
        ///mass, and metallicity.
        ///
        ///Must be an increasing monotonic function
        virtual double interp_param_to_age(
            ///The interpolation parameter for which the age is needed.
            double interp_param,

            ///The stellar mass for which the interpolation parameter is
            ///needed in \f$M_\odot\f$.
            double mass,

            ///The stellar metallicity for which the interpolation
            ///parameter is needed in \f$[Fe/H]\f$.
            double metallicity) const;

    public:
        ///\brief Default constructor (only useful for derived classes which
        ///do not use the interpolation).
        EvolvingStellarQuantity() {}

        ///Create an evolving quantity that interpolates to the given mass.
        EvolvingStellarQuantity(
            ///The stellar mass to interpolate to in \f$M_\odot\f$
            double mass, 

            ///The stellar metallicity (\f$[Fe/H]f$) to interpolate to.
            double metallicity, 

            ///The masses for which evolution tracks are given in
            /// \f$M_\odot\f$
            const std::valarray<double> &track_masses,

            ///The metallicities (\f$[Fe/H]\f$) for which evolution tracks
            ///are given.
            const std::valarray<double> &track_metallicities,

            ///The evolution tracks of the relevant quantity on the grid
            ///defined by \p track_masses and \p track_metallicities. The
            //mass index varies faster.
            const std::vector<const OneArgumentDiffFunction *> 
            &evolution_tracks,

            ///Whether the track uses log(age) as the independent argument
            ///instead of age.
            bool log_age=true,

            ///Whether this is a quantity that is identically zero below some
            ///age and turns on afterwards
            bool starts_zero=false);

        ///\brief Prepare the quantity for interpolation around the given
        ///age.
        ///
        ///After calling this method, requesting values or derivatives
        ///outside the range of the continuous region containing this age
        ///(see ::discontinuities) fails an assert.
        virtual void select_interpolation_region(double age) const;

        ///Return the value the quantity takes at the given age in Gyr.
        virtual double operator()(double age) const
        {return interpolate(age);}

        ///Return the age derivative of the quantity at the given age in Gyr.
        virtual const FunctionDerivatives *deriv(double age) const;

        ///The largest age for which the quantity can be interpolated in Gyr.
        virtual double range_high() const {return __max_age;}

        ///\brief The smallest age for which the quantity can be interpolated
        ///in Gyr.
        virtual double range_low() const {return __min_age;}

        ///The ages at which the quantity may be discontinuous.
        virtual const std::vector<double> &discontinuities() const
        {return __interp_grid_change_ages;}

        ///\brief The upper bound of the current interpolation region (over
        ///which the quantity is guaranteed continuous).
        ///
        ///The current interpolation region is set either
        ///by ::set_starting_interpolation_age or
        ///by ::enable_next_interpolation_region.
        virtual double next_discontinuity() const
        {return *__next_grid_change_age;}

        ///\brief Set up the interpolation over the next interpolation region
        ///(between consecutive discontinuities.)
        virtual void enable_next_interpolation_region() const;

        ///An iterator over the ages (in Gyr) where the quantity takes the
        ///given y value.
        InterpSolutionIterator crossings(double =0) const
        {
            throw Core::Error::Runtime(
                "Called EvolvingStellarQuantity::crossings, which are ill"
                "defined."
            );
        }
    }; //End EvolvingStellarQuantity class declaration.

} //End StellarEvolution namespace

#endif
