/**\file
 *
 * \brief Defines the classes for generating stellar evolution interpolators
 * from the MESA tracks. 
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __MESAIO_H
#define __MESAIO_H

#include "Interpolator.h"
#include "../../IO/IOUtil.h"
#include "../../Core/Common.h"
#include "../../Core/AstronomicalConstants.h"
#include "../../Core/Error.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cassert>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>

namespace StellarEvolution {

    ///\brief A namespace to isolate all MESA related entities, in order to avoid
    ///conflicts with other StellarEvolution implentations (e.g. YREC).
    ///
    ///\ingroup StellarSystem_group
    ///
    ///\todo Tune default nodes and smoothing for good interpolation.
    namespace MESA {

        ///Names for the interesting columns in a MESA track.
        enum Column {
            ///The total mass of the star in \f$M_\odot\f$.
            MTRACK,

            ///Stellar age in years.
            AGE,

            ///Log10 of the radius of the star in \f$R_\odot\f$.
            LOG_RSTAR,

            ///The radius of the star in \f$R_\odot\f$.
            RSTAR,

            ///Log10 of the luminosity of the star in \f$\L_\odot\f$.
            LOG_LSTAR,

            ///The luminosity of the star in \f$\L_\odot\f$.
            LSTAR,

            ///Mass of the radiative core in \f$M_\odot\f$.
            MRAD,

            ///Radius of the radiative core in \f$R_\odot\f$.
            RRAD,

            ///\brief Moment of inertia of the convective envelope in
            /// \f$\mathrm{kg}\cdot\mathrm{m}^2\cdot\mathrm{rad}/\mathrm{s}\f$.
            ICONV,

            ///\brief Moment of inertia of the radiative core in
            /// \f$\mathrm{kg}\cdot\mathrm{m}^2\cdot\mathrm{rad}/\mathrm{s}\f$.
            IRAD,


            ///The total number of interesting columns 
            NUM_COLUMNS
        };

        ///\brief A class which parses the header of a MESA evolution track.
        ///
        ///\ingroup StellarSystem_group
        class Header {
            private:
                ///\brief The names of the columns in the track, indexed by
                ///MESA::Column.
                std::vector<std::string> __column_names;

                ///\brief The column numbers of the interesting quantities,
                ///indexed by MESA::Column.
                std::vector<int> __column_numbers;

                ///\brief Checks that the next line in the input stream
                ///consists of sequential numbers starting with 1.
                ///
                ///Skips leading empty lines, incrementing line_number
                ///appropriately.
                void read_column_numbers(std::istream &track,
                        const std::string &filename, unsigned &line_number);

                ///Sets the column names.
                void set_column_names();
            public:
                ///Parse the header information from the given track stream.
                Header(std::ifstream &track, const std::string &filename);

                ///The column number corresponding to the given quantity.
                int get_column(MESA::Column quantity) const
                {return __column_numbers[quantity];}

                ///Returns all column numbers at once.
                const std::vector<int> &get_all_columns() const
                {return __column_numbers;}
        };

        

        ///\brief An iterator over the list of extracted tracks.
        ///
        ///\ingroup StellarSystem_group
        class EvolutionIterator {
        public:
            ///\brief Create an iterator, which must have all its *_iter 
            ///members set before it can be used.
            EvolutionIterator() : quantity_iter(NUM_QUANTITIES) {}

            ///Iterator over the masses of the tracks.
            std::list<double>::iterator mass_iter;

            ///Iterator over the masses of the tracks.
            std::list<double>::iterator metallicity_iter;

            ///Iterator over the array of ages of the tracks.
            std::list< std::valarray<double> >::iterator age_iter;

            ///Iterators over the arrays of the track quantities.
            std::vector< std::list< std::valarray<double> >::iterator >
                quantity_iter;

            ///Copy orig to *this.
            EvolutionIterator(const EvolutionIterator &orig) :
                mass_iter(orig.mass_iter),
                metallicity_iter(orig.metallicity_iter),
                age_iter(orig.age_iter),
                quantity_iter(orig.quantity_iter)
            {}

            ///Copy rhs to *this.
            EvolutionIterator &operator=(const EvolutionIterator &rhs);

            ///Advance all iterators to the next track.
            EvolutionIterator &operator++();

            ///Advance all iterators to the next track.
            EvolutionIterator operator++(int)
            {EvolutionIterator result(*this); ++(*this); return result;}

            ///\brief Is RHS at the same position as this?
            ///
            ///Assumes that the iteration is over the same list of tracks.
            bool operator==(const EvolutionIterator &rhs)
            {return mass_iter == rhs.mass_iter;}

            ///\brief Is RHS at a different position than this?
            ///
            ///Assumes that the iteration is over the same list of tracks.
            bool operator!=(const EvolutionIterator &rhs)
            {return !((*this) == rhs);}
        };

        ///\brief A stellar evolution interpolator based on the MESA tracks.
        ///
        ///\ingroup StellarSystem_group
        class Interpolator : public StellarEvolution::Interpolator {
        private:

            ///\brief The value at the indexed from Column is the
            ///StellarEvolution::QuantityID to which this column is
            ///converted.
            ///
            ///Columns not corresponding to a quantity are set to
            ///StellarEvolution::NUM_QUANTITIES.
            static const std::vector<QuantityID> __column_to_quantity;

            ///\brief The default amount of smoothing to use for each 
            ///quantity. See StellarEvolution::Interpolator::create_from.
            static const std::vector<double> __default_smoothing;

            ///\brief The default number of node to use for each 
            ///quantity. See StellarEvolution::Interpolator::create_from.
            static const std::vector<int> __default_nodes;

            ///The masses of the available tracks in the order read.
            std::list<double> __mass_list;

            ///The metallicities of the available tracks in the order read.
            std::list<double> __metallicity_list;

            ///The ages at which each track is tabulated.
            std::list< std::valarray<double> > __track_ages;

            ///\brief A structure holding all interesting quantities from the
            ///MESA tracks except age.
            ///
            ///The innermost valarray is a particular parameter for a
            ///particular track as a function of age.
            ///
            ///The middle list holds the same quantity for all tracks and is
            ///in the same order as #__mass_list.
            ///
            ///The outside array is indexed by quantity identified by
            ///StellarEvolution::QuantityID.
            std::vector< std::list< std::valarray<double> > >
                __track_quantities;

            ///\brief Parse the mass (in \f$M_\odot\f$) and [Fe/H] from a
            ///track filename.
            ///
            ///If the filaneme follows the expected pattern, add the parsed 
            ///values to ::__mass_list and ::__metallicity_list respectively 
            ///and return true. If the filename is not formatted as expected
            ///return false and leave ::__mass_list and ::__metallicity_list
            ///unchanged.
            bool parse_model_file_name(const std::string &filename);

            ///Reads a single evolution track file
            void read_model_file(const std::string &filename);

#ifndef NDEBUG
            ///Output the current masses metallicities and age ranges.
            void log_current_age_ranges() const;
#endif

            ///\brief Verify that the track masses and metallicities form a
            ///grid and return the grid.
            void get_mass_metallicity_grid(
                ///Output argument: the list of stellar masses in the grid
                ///(sorted and unique values only).
                std::valarray<double> &masses,

                ///Output argument: the list of stellar metallicites in the
                ///grid (sorted and unique values only).
                std::valarray<double> &metallicities
            );

            ///\brief Returns an EvolutionIterator pointing to the beginning
            ///of all quantities.
            EvolutionIterator begin();

            ///\brief Returns an EvolutionIterator pointing to the end of all
            ///quantities.
            EvolutionIterator end();

            ///Moves source to right before destination
            void move(EvolutionIterator &dest, EvolutionIterator &source);

            ///Sorts the data by mass and metallicity.
            void sort_tracks();
        public:
            ///Default constructor, use load_state to get a working
            ///interpolator.
            Interpolator() {};

            ///\brief Creates a stellar evolution interpolator based on
            ///evolution tracks computed with MESA.
            Interpolator(
                ///The directory containing the MESA evolution tracks
                const std::string &model_directory,

                ///How much to smooth each stellar evolution quantity when
                ///fitting. Use NaN for no smoothing. See
                ///StellarEvolution::QuantityID for the order of the
                ///quantities.
                const std::vector<double> &smoothing = __default_smoothing,

                ///How many nodes to use when smoothing each stellar
                ///evolution quantity. Ignored if smooth_conv_inertia is
                ///NaN - no smoothing.
                ///
                ///Negative values result in using min(-nodes[i],
                ///number of tabulated ages for each track).
                const std::vector<int> &nodes = __default_nodes
            );
        };

    } //End MESA namespace.

} //End StellarEvolution namespace.

///Civilized output of mesa column names.
std::ostream &operator<<(std::ostream &os,
                         StellarEvolution::MESA::Column col);

#endif
