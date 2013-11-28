/**\file
 *
 * \brief Defines the classes for generating stellar evolution interpolators
 * from the MESA tracks. 
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __MESAIO_H
#define __MESAIO_H

#include "StellarEvolution.h"
#include "IOUtil.h"
#include "Common.h"
#include "AstronomicalConstants.h"
#include "Error.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <dirent.h>
#include <string>
#include <assert.h>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>

///\brief A namespace to isolate all MESA related entities, in order to avoid
///conflicts with other StellarEvolution implentations (e.g. YREC).
///
///\ingroup StellarSystem_group
///
///\todo Tune default nodes and smoothing for good interpolation.
namespace MESA {

	///Names for the interesting columns in a MESA track.
	enum Column {
		AGE,      ///< Stellar age in years.
		LOG_RSTAR,///< Log10 of the radius of the star in \f$R_\odot\f$.
		RSTAR=LOG_RSTAR,///< We convert Log10(R*) to R*
		LOG_LSTAR,///< Log10 of the luminosity of the star in \f$\L_\odot\f$.
		LSTAR=LOG_LSTAR,///< We convert Log10(L*) to L*
		MRAD,     ///< Mass of the radiative core in \f$M_\odot\f$.
		RRAD,     ///< Radius of the radiative core in \f$R_\odot\f$.

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
			///The mass of the star in the track being read.
			double __track_mass;

			///\brief The names of the columns in the track, indexed by
			///MESA::Column.
			std::vector<std::string> __column_names;

			///\brief The column numbers of the interesting quantities,
			///indexed by MESA::Column.
			std::vector<int> __column_numbers;

			///\brief Checks that the next line in the input stream consists
			///of sequential numbers starting with 1.
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

			///The stellar mass (in \f$M_\odot\f$) of the track.
			double get_mass() const {return __track_mass;}

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
		///\brief Create an iterator, which must have all its *_iter members
		///set before it can be used.
		EvolutionIterator() : quantity_iter(NUM_COLUMNS) {}

		///Iterator over the masses of the tracks.
		std::list<double>::iterator mass_iter;

		///Iterators over the arrays of the track quantities.
		std::valarray< std::list< std::valarray<double> >::iterator >
			quantity_iter;

		///Copy orig to *this.
		EvolutionIterator(const EvolutionIterator &orig) :
			mass_iter(orig.mass_iter), quantity_iter(orig.quantity_iter) {}

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
		{return mass_iter==rhs.mass_iter;}

		///\brief Is RHS at a different position than this?
		///
		///Assumes that the iteration is over the same list of tracks.
		bool operator!=(const EvolutionIterator &rhs)
		{return !((*this)==rhs);}
	};

	///\brief A stellar evolution interpolator based on the MESA tracks.
	///
	///\ingroup StellarSystem_group
	class Evolution : public StellarEvolution {
	private:
		///The masses of the available tracks.
		std::list<double> __mass_list;

		///\brief A structure holding all interesting quantities from all
		///MESA tracks.
		///
		///The innermost valarray is a particular parameter for a
		///particular track as a function of age. 
		///
		///The middle list holds the same quantity for all tracks and is
		///in the same order as #__mass_list.
		///
		///The outside array is indexed by quantity identified by
		///MESA::Column.
		std::valarray< std::list< std::valarray<double> > >
			track_quantities;

		///Reads a single evolution track file
		void read_model_file(const std::string &filename);

		///\brief Returns an EvolutionIterator pointing to the beginning
		///of all quantities.
		EvolutionIterator begin();

		///\brief Returns an EvolutionIterator pointing to the end of all
		///quantities.
		EvolutionIterator end();

		///Moves source to right before destination
		void move(EvolutionIterator &dest, EvolutionIterator &source);

		///Sorts the data by mass.
		void sort_masses();
	public:
		///Default constructor, use load_state to get a working
		///interpolator.
		Evolution() {};

		///\brief Creates a stellar evolution interpolator based on
		///evolution tracks computed with MESA.
		Evolution(
				///The directory containing the MESA evolution tracks
				const std::string &model_directory,

				///How much to smooth the moment of inertia of the
				///convective zone when fitting. Use NaN for no smoothing.
				double smooth_conv_inertia=0,

				///How much to smooth the moment of inertia of the
				///radiative zone of the star when fitting. Use NaN for no
				///smoothing.
				double smooth_rad_inertia=2,

				///How much to smooth the mass in the radiative zone when
				///fitting. Use NaN for no smoothing.
				double smooth_rad_mass=2,
				
				///How many nodes to use when smoothing the moment of inertia
				///of the convective zone (ignored if #smooth_conv_inertia is
				///NaN - no smoothing).
				///
				///Negative values result in using min(-conv_inertia_nodes,
				///number of tabulated ages for each track).
				int conv_inertia_nodes=-1000,

				///How many nodes to use when smoothing the moment of inertia
				///of the radiative zone (ignored if #smooth_rad_inertia is
				///NaN - no smoothing).
				///
				///Negative values result in using min(-rad_inertia_nodes,
				///number of tabulated ages for each track).
				int rad_inertia_nodes=-1000,

				///How many nodes to use when smoothing the mass of the
				///radiative zone (ignored if #smooth_rad_inertia is NaN - no
				///smoothing).
				///
				///Negative values result in using min(-rad_mass_nodes,
				///number of tabulated ages for each track).
				int rad_mass_nodes=-1000);
	};
};

///Civilized output of mesa column names.
std::ostream &operator<<(std::ostream &os, MESA::Column col);

#endif
