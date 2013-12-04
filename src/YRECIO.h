/**\file
 *
 * \brief Defines the classes for generating stellar evolution interpolators
 * from the YREC tracks. 
 * 
 * \ingroup StellarSystem_group
 *
 * \todo Make it follow the same namespace scheme as MESA
 */

#ifndef __YRECIO_H
#define __YRECIO_H

#include "StellarEvolution.h"
#include "Common.h"
#include "AstronomicalConstants.h"
#include "Error.h"
#include <fstream>
#include <dirent.h>
#include <string>
#include <assert.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>

///\brief A class which parses the header of a YREC evolution track.
///
///\ingroup StellarSystem_group
class YRECHeader {
private:
	///The masss of the tracks provided on construction in \f$M_\odot\f$.
	double track_mass;

	int age_col, ///< The index of the age column in the track.
		radius_col, ///< The index of the radius column in the track.

		///\brief The index of the mass in the convective envelope column
		///in the track.
		envelope_mass_col,

		///The index of the core-envelope boundary column in the track.
		rad_conv_boundary_col, 

		///The index of the radiative core inertia column in the track.
	    rad_inertia_col,

		///The index of the convective envelope inertia column in the track.
		conv_inertia_col,

		///The index of the lg(luminosity) column in the track.
		log_luminosity_col;
public:
	///Parse the header information from the given track stream.
	YRECHeader(std::ifstream &track, const std::string &filename);

	///The stellar mass (in \f$M_\odot\f$) of the track.
	double get_mass() const {return track_mass;}

	///The column index within the track that contains the track ages.
	int get_age_col() const {return age_col;}

	///The column index within the track that contains the stellar radii.
	int get_radius_col() const {return radius_col;}

	///The column index within the track that contains the lg(luminosity).
	int get_log_luminosity_col() const {return log_luminosity_col;}

	///The column index for the envelope mass.
	int get_envelope_mass_col() const {return envelope_mass_col;}

	///The column index for the convective-radiative boundary radius.
	int get_core_boundary_col() const {return rad_conv_boundary_col;}

	///The column index for the moment of inertia of the radiative core.
	int get_rad_inertia_col() const {return rad_inertia_col;}

	///The column index for the moment of inertia of the convective envelope.
	int get_conv_inertia_col() const {return conv_inertia_col;}
};

///\brief An iterator over the list of extracted tracks.
///
///\ingroup StellarSystem_group
class EvolutionIterator {
public:
	///\brief Create an iterator, which must have all its *_iter members set
	///before it can be used.
	EvolutionIterator() {}

	///Iterator over the masses of the tracks.
	std::list<double>::iterator mass_iter;

	///Iterator over the arrays of ages of the tracks.
	std::list< std::valarray<double> >::iterator age_iter,

		///Iterator over the arrays of stellar radii of the tracks.
		radius_iter,

		///Iterator over the arrays of stellar lg(luminosity) of the tracks.
		luminosity_iter,

		///Iterator over the arrays of core masses of the tracks.
		rad_mass_iter,

		///\brief Iterator over the arrays of core-envelope boundaries of the
		///tracks.
		core_boundary_iter,

		///\brief Iterator over the arrays of convective envelope moments of
		///inertia of the tracks.
		conv_inertia_iter,
		
		///\brief Iterator over the arrays of radiative core moments of
		///inertia of the tracks.
		rad_inertia_iter;

	///Copy orig to *this.
	EvolutionIterator(const EvolutionIterator &orig);

	///Copy rhs to *this.
	EvolutionIterator &operator=(const EvolutionIterator &rhs);

	///Advance all iterators to the next track.
	EvolutionIterator &operator++();

	///Advance all iterators to the next track.
	EvolutionIterator operator++(int);

	///\brief Is RHS at the same position as this?
	///
	///Assumes that the iteration is over the same list of tracks.
	bool operator==(const EvolutionIterator &rhs);

	///\brief Is RHS at a different position than this?
	///
	///Assumes that the iteration is over the same list of tracks.
	bool operator!=(const EvolutionIterator &rhs) {return !((*this)==rhs);}
};

///\brief A stellar evolution interpolator based on the YREC tracks.
///
///It relies on 10 tracks, all at solar metallicitity and none of them going
///past 10Gyr stellar age. The masses of the tracks are (in \f$M_\odot\f$):
///0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.15, and 1.2. 
///
///\ingroup StellarSystem_group
class YRECEvolution : public StellarEvolution {
private:
	///The masses of the available tracks.
	std::list<double> mass_list;

	std::list< std::valarray<double> > 
		ages, ///< The ages in each track.
		radii, ///< The stellar radii in each track.	
		luminosities, ///< The luminosities in each track
		rad_masses,///< The masses of the core in each track

		///The core-envelope boundary radius in each track
		core_boundaries,

		///The moments of inertia of the convective envelope
		conv_inertias,
		
		///The moments of inertia of the radiative zone of the star
		rad_inertias;

	///Reads a single evolution track file
	void read_model_file(const std::string &filename);

	///\brief Returns an EvolutionIterator pointing to the beginning of all
	///quantities.
	EvolutionIterator begin();

	///\brief Returns an EvolutionIterator pointing to the end of all
	///quantities.
	EvolutionIterator end();

	///Moves source to right before destination
	void move(EvolutionIterator &dest, EvolutionIterator &source);

	///Sorts the data by mass.
	void sort_masses();
public:
	///Default constructor, use load_state to get a working interpolator.
	YRECEvolution(){};

	///\brief Creates a stellar evolution interpolator based on evolution
	///tracks computed with YREC.
	YRECEvolution(
		///The directory containing the YREC evolution tracks
		const std::string &model_directory,

		///How much to smooth the moment of inertia of the convective zone
		///when fitting.
		double smooth_conv_inertia=0,

		///How much to smooth the moment of inertia of the radiative zone of
		///the star when fitting.
		double smooth_rad_inertia=2,

		///How much to smooth the mass in the radiative zone when fitting.
		double smooth_rad_mass=2,
		
		///How many nodes to use when smoothing the moment of inertia of the
		///convective zone (ignored if #smooth_conv_inertia is NaN - no
		///smoothing).
		///
		///Negative values result in using min(-conv_inertia_nodes, number of
		///tabulated ages for each track).
		int conv_inertia_nodes=-1000,

		///How many nodes to use when smoothing the moment of inertia of the
		///radiative zone (ignored if #smooth_rad_inertia is NaN - no
		///smoothing).
		///
		///Negative values result in using min(-rad_inertia_nodes, number of
		///tabulated ages for each track).
		int rad_inertia_nodes=-1000,

		///How many nodes to use when smoothing the mass of the radiative
		///zone (ignored if #smooth_rad_inertia is NaN - no smoothing).
		///
		///Negative values result in using min(-rad_mass_nodes, number of
		///tabulated ages for each track).
		int rad_mass_nodes=-1000);
};

#endif
