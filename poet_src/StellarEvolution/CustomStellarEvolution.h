/**\file
 *
 * \brief Defines the class which implements a custom single track stellar
 * evolution.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __CUSTOM_STELLAR_EVOLUTION_H
#define __CUSTOM_STELLAR_EVOLUTION_H

#include "StellarEvolution.h"
#include "IOUtil.h"
#include <fstream>
#include <string>
#include <limits>

namespace StellarEvolution {

    ///A namespace to isolate all custom stellar evolution related entities.
    ///
    ///\ingroup StellarSystem_group
    namespace CustomStellarEvolution {

        ///Tags for the columns in input stellar tracks.
        enum Columns {
            ///\brief Moment of inertia of the convective envelope of the 
            ///star in \f$M_\odot R_\odot^2\f$.
            ICONV,

            ///\brief Moment of inertia of the radiative core of the star in 
            /// \f$M_\odot R_\odot^2\f$.
            IRAD,

            RSTAR,///< Radius of the star in \f$R_\odot\f$.

            ///\brief Radius of the stellar core in \f$R_\odot\f$ (low mass
            ///stars only).
            RRAD,

            ///\brief Mass of the stellar core in \f$M_\odot\f$ (low mass
            ///stars only).
            MRAD,

            ///Luminosity of the star in \f$L_\odot\f$.
            LSTAR,

            ///Age of the star in Gyr.
            AGE,

            ///A column which is not needed to interpolate the evolution.
            SKIP,

            ///The number of different input quantities supported.
            NUM_TRACK_QUANTITIES = SKIP
        };

        ///\brief A stellar evolution interpolator using only a single track,
        ///assumed to apply to all stars.
        class Interpolator : public StellarEvolution::Interpolator {
        private:
            ///\brief Reads the stellar evolution track from a file returning
            ///the array of the quantities indexed by the column tag.
            std::valarray< std::valarray<double> > read_track(
                ///The name of the file containing the track.
                const std::string &filename,

                ///A list of the columns in the track
                const std::vector<Columns> &format
            ) const;
        public:
            ///\brief Creates an interpolator from a given track.
            Evolution(
                ///The filename of the track.
                const std::string &filename,

                ///A list of the columns in the track
                const std::vector<Columns> &format,

                ///The smoothing to apply to each column (NaN for no 
                ///smoothing).
                ///
                ///Entries for AGE and SKIP are ignored.
                const std::vector<double> &smoothing,

                ///The number of nodes to use when smoothing each column.
                ///
                ///This is ignored for columns for which the corresponding
                ///entry in smoothing is NaN
                ///
                ///Negative entries result in using the smaller of the
                ///absolute value of the entry and three times the number of
                ///track points.
                const std::vector<int> &nodes,

                ///Should interpolation be done vs. log(age) instead of age
                ///for each quantity?
                const std::vector<bool> &vs_log_age,

                ///Should interpolation be done of log(quantity) instead of
                ///quantity for each quantity?
                const std::vector<bool> &log_quantity
            );
        };

    } //End CustomStellarEvolution namespace.

} //End StellarEvolution namespace.

#endif
