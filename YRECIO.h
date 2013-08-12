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

class YRECHeader {
private:
	///The masss for which the tracks provided on construction refers to.
	double track_mass;

	///The various columns necessary to calculate the orbital evolution.
	int age_col, radius_col, envelope_mass_col, rad_conv_boundary_col, 
	    rad_inertia_col, conv_inertia_col, log_luminosity_col;
public:
	///Parse the header information from the given track stream.
	YRECHeader(std::ifstream &track, const std::string &filename);

	///Return the stellar mass (in solar masses) for which this tracks 
	///applies
	double get_mass() const {return track_mass;}

	///Return the column index within the track stream that contains the
	///track ages.
	int get_age_col() const {return age_col;}

	///Return the column index within the track stream that contains the
	///track stellar radii.
	int get_radius_col() const {return radius_col;}

	///Return the column index within the track stream that contains the
	///track stellar log luminosities.
	int get_log_luminosity_col() const {return log_luminosity_col;}

	///Returns the column index for the envelope mass.
	int get_envelope_mass_col() const {return envelope_mass_col;}

	///Returns the column index for the convective-radiative boundary 
	///radius.
	int get_core_boundary_col() const {return rad_conv_boundary_col;}

	///Returns the column index for the moment of inertia of the 
	///radiative core.
	int get_rad_inertia_col() const {return rad_inertia_col;}

	///Returns the column index for the moment of inertia of the 
	///convective envelope.
	int get_conv_inertia_col() const {return conv_inertia_col;}
};

class EvolutionIterator {
public:
	EvolutionIterator() {}
	std::list<double>::iterator mass_iter;
	std::list< std::valarray<double> >::iterator age_iter, radius_iter,
		luminosity_iter, rad_mass_iter, core_boundary_iter,
		conv_inertia_iter, rad_inertia_iter;
	EvolutionIterator(const EvolutionIterator &orig);
	EvolutionIterator &operator=(const EvolutionIterator &rhs);
	EvolutionIterator &operator++();
	EvolutionIterator operator++(int);
	bool operator==(const EvolutionIterator &rhs);
	bool operator!=(const EvolutionIterator &rhs) {return !((*this)==rhs);}
};

class YRECEvolution : public StellarEvolution {
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & boost::serialization::base_object<StellarEvolution>(*this);
		ar & mass_list;
		ar & ages & radii & luminosities & rad_masses &
			core_boundaries & conv_inertias & rad_inertias;
	}
private:
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

	///Returns an EvolutionIterator pointing to the beginning of all
	///quantities.
	EvolutionIterator begin();

	///Returns an EvolutionIterator pointing to the end of all
	///quantities.
	EvolutionIterator end();

	///Moves source to right before destination
	void move(EvolutionIterator &dest, EvolutionIterator &source);

	///Sorts the data by mass.
	void sort_masses();
public:
	YRECEvolution(){};
	///Creates a stellar evolution variable based on evolution tracks coputed
	///with YREC.
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
		double smooth_rad_mass=2);

	/*Loads data from serialization. Only call this on objects initialized
	 * with the default constructor. */
	void load_state(std::string filename="../interp_state_data");

	/*Only call this on objects NOT initialized using the default constructor
	 * (otherwise it has no data to save). Serializes state to file.
	 * Recursively saves data of YRECEvolution and every class it depends on:
	 * StellarEvolution, InterpolatingFunctionALGLIB, OneArgumentDiffFunction,
	 * OneArgumentFunction, spline1dinterpolant, and
	 * _spline1dinterpolant_owner. In _spline1dinterpolant_owner, serialize()
	 * serializes everything EXCEPT p_struct->x.data and p_struct->y.data,
	 * because I have no idea what they are. .*/
	void save_state(std::string filename="../interp_state_data") const;
};

#endif
