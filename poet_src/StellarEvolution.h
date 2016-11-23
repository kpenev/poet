/**\file
 *
 * \brief Defines the StellarEvolution class needed for interpolating among
 * stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __STELLAR_EVOLUTION_H
#define __STELLAR_EVOLUTION_H

#include "StellarZone.h"
#include "EvolvingStellarQuantity.h"
#include "Error.h"
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

///\brief A class that interpolates among stellar evolution tracks.
///
///Uses a  set of pre-computed evolution tracks to generate 
///inrpolating functions that represents reasonably well the evolution of
///various properties of an arbitrary mass star as a function of age.
///
///At present the only implementing class is YRECEvolution based on a set of
///YREC tracks.
///
///\ingroup StellarSystem_group
class StellarEvolution {
private:
	///Serialize the found interpolation.
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int) {
		if(!std::isfinite(core_formation)) core_formation=-1;
		ar & track_masses;
		ar & interpolated_radius & interpolated_luminosity
		   & interpolated_conv_inertia & interpolated_rad_inertia
		   & interpolated_rad_mass & interpolated_core_env_boundary;
		ar & mass_break & low_age_scaling & high_age_scaling
		   & extrapolate_low & extrapolate_high & core_formation;
		if(core_formation<0) core_formation=Inf;
		std::cerr << "core_formation=" << core_formation << std::endl;
	}

	///The stellar masses for which evolution tracks are available in
	/// \f$M_\odot\f$
	const std::valarray<double> *track_masses;

	///\brief Interpolates the evulotion of the radius of the star for each
	///track in \f$R_\odot\f$.
	std::list<const OneArgumentDiffFunction*> interpolated_radius,

		///\brief Interpolates the evulotion of the convective zone moment of
		///inertia for each track in \f$M_\odot \cdot R_\odot^2\f$.
		interpolated_conv_inertia,
		
		///\brief Interpolates the evulotion of the radiative zone moment of
		///inertia for each track\f$M_\odot \cdot R_\odot^2\f$.
		interpolated_rad_inertia, 

		///\brief Interpolates the evulotion of the radiative zone mass for
		///each track in \f$M_\odot\f$.
		interpolated_rad_mass,
		
		///\brief Interpolates the evulotion of the convective-radiative
		///boundary for each track in \f$R_\odot\f$.
		interpolated_core_env_boundary,

		///\brief Interpolates the evulotion of the luminosity for each
		///track in \f$L_\odot\f$.
		interpolated_luminosity;

	double
		///The break between high and low mass models in \f$M_\odot\f$.
		mass_break,

		///The age of a low mass model is scaled by mass to this power
		low_age_scaling,

		///The age of a low mass model is scaled by mass to this power
		high_age_scaling,

		///\brief How far to extrapolate low mass models.
		///
		///Exrapolating a low mass model by this fractional amount of the
		///final age is allowed.
		extrapolate_low,

		///\brief How far to extrapolate low mass models.
		///
		///Exrapolating a high mass model by this fractional amount of the
		///final age is allowed.
		extrapolate_high,
		
		///The age at which the core starts forming in Gyr.
		core_formation;

public:
	///\brief Construct an object that can be set to interpolate between
	///tabulated evolution tracks.
	StellarEvolution() : track_masses(NULL) {}

	///\brief Creates a fully functional stellar evolution interpolator.
	StellarEvolution(
		///The stellar masses (in \f$M_\odot\f$) for which evolution tracks
		///are tabulated.
		const std::valarray<double> &tabulated_masses,
		
		///A set of ages for each track in Gyr.
		const std::list< std::valarray<double> > &tabulated_ages,

		///A set of stellar radii for each age of each track in
		/// \f$R_\odot\f$.
		const std::list< std::valarray<double> > &tabulated_radii,

		///A set of moments of inertia of the stellar convective zone for 
		///each age of each track in \f$ M_\odot \cdot R_\odot^2\f$.
		const std::list< std::valarray<double> > &tabulated_conv_inertia,

		///A set of moments of inertia of the radiative zone for each age
		///of each track in \f$ M_\odot \cdot R_\odot^2\f$.
		const std::list< std::valarray<double> > &tabulated_rad_inertia,

		///A set of masses of the stellar radiative zone for each age of
		///each track in \f$M_\odot\f$. 
		const std::list< std::valarray<double> > &tabulated_rad_mass,

		///A set of radii (in \f$R_\odot\f$) for the convective-envelope 
		///boundary for each age of each track.
		const std::list< std::valarray<double> > 
			&tabulated_core_env_boundary,

		///How much to smooth the stellar radius when fitting. Use NaN for no
		///smoothing.
		double smooth_radius,

		///How much to smooth the moment of inertia of the convective zone
		///when fitting.
		double smooth_conv_inertia,

		///How much to smooth the moment of inertia of the stellar core
		///when fitting.
		double smooth_rad_inertia,

		///How much to smooth the mass in the radiative zone when fitting.
		double smooth_rad_mass,

		///How much to smooth the radius of the radiative zone when fitting.
		double smooth_core_env_boundary, 

		///How many nodes to use when smoothing the stellar radius (ignored
		///if smooth_radius is NaN - no smoothing).
		int radius_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///convective zone (ignored if smooth_conv_inertia is NaN - no
		///smoothing).
		int conv_inertia_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///radiative zone (ignored if smooth_rad_inertia is NaN - no
		///smoothing).
		int rad_inertia_nodes,

		///How many nodes to use when smoothing the mass of the radiative
		///zone (ignored if smooth_rad_inertia is NaN - no smoothing).
		int rad_mass_nodes,

		///How many nodes to use when smoothing the radius of the radiative
		///zone (ignored if smooth_core_env_boundary is NaN - no smoothing).
		int core_env_boundary_nodes, 

		///A set of lg(luminosities) (in \f$L_\odot\f$) for each age of
		///each track. Can be omitted if lominosity interpolation is not
		///necessary.
		const std::list< std::valarray<double> > &tabulated_luminosities=
			std::list< std::valarray<double> >(),

		///How much to smooth the luminosities when fitting.
		double smooth_luminosities=NaN,

		///How many nodes to use when smoothing the luminosities (ignored if
		///smooth_luminosities is NaN - no smoothing).
		int luminosities_nodes=0,

		///The mass above which the stars are considered 
		///high mass in \f$M_\odot\f$.
		double max_low_mass=1.075,

		///When interpolating the age of each low mass model is scaled by
		///the mass to negative this power in order to make mass
		///dependence smoother.
		double low_mass_age_scaling=2.5,

		///When interpolating the age of each high mass model is scaled 
		///by the mass to negative this power in order to make mass
		///dependence smoother.
		double high_mass_age_scaling=3.0,

		///Low mass models are included in the mass interpolation only if
		///the required age is no larger than the maximum tabulated age 
		///times this factor.
		double low_mass_extrapolate=1.01,

		///High mass models are included in the mass interpolation only
		///if the required age is no larger than the maximum tabulated 
		///age times this factor.
		double high_mass_extrapolate=1.01) {
			interpolate_from(tabulated_masses, 
					tabulated_ages, 
					tabulated_radii, 
					tabulated_conv_inertia,
					tabulated_rad_inertia, 
					tabulated_rad_mass,
					tabulated_core_env_boundary, 

					smooth_radius,
					smooth_conv_inertia, 
					smooth_rad_inertia,
					smooth_rad_mass,
					smooth_core_env_boundary,

					radius_nodes,
					conv_inertia_nodes,
					rad_inertia_nodes,
					rad_mass_nodes,
					core_env_boundary_nodes,

					tabulated_luminosities, 
					smooth_luminosities,
					luminosities_nodes,

					max_low_mass, low_mass_age_scaling,
					high_mass_age_scaling, low_mass_extrapolate,
					high_mass_extrapolate);}

	///Fully setup an object created by the default constructor.
	void interpolate_from(
		///The stellar masses (in \f$M_\odot\f$) for which evolution tracks
		///are tabulated.
		const std::valarray<double> &tabulated_masses,
		
		///A set of ages for each track in Gyr.
		const std::list< std::valarray<double> > &tabulated_ages,

		///A set of stellar radii for each age of each track in
		/// \f$R_\odot\f$.
		const std::list< std::valarray<double> > &tabulated_radii,

		///A set of moments of inertia of the stellar convective zone for 
		///each age of each track in \f$M_\odot\cdot R_\odot^2\f$.
		const std::list< std::valarray<double> > &tabulated_conv_inertia,

		///A set of moments of inertia of the radiative zone for each age
		///of each track in \f$M_\odot\cdot R_\odot^2\f$.
		const std::list< std::valarray<double> > &tabulated_rad_inertia,

		///A set of masses of the stellar radiative zone for each age of
		///each track in \f$M_\odot\f$
		const std::list< std::valarray<double> > &tabulated_rad_mass,

		///A set of radii (in \f$R_\odot\f$) for the convective-envelope 
		///boundary for each age of each track.
		const std::list< std::valarray<double> > 
			&tabulated_core_env_boundary,

		///How much to smooth the stellar radius when fitting. Use NaN for no
		///smoothing.
		double smooth_radius,

		///How much to smooth the moment of inertia of the convective zone
		///when fitting.
		double smooth_conv_inertia,

		///How much to smooth the moment of inertia of the radiative zone
		///when fitting.
		double smooth_rad_inertia,

		///How much to smooth the mass in the radiative zone when fitting.
		double smooth_rad_mass,

		///How much to smooth the radius of the radiative zone when fitting.
		double smooth_core_env_boundary, 

		///How many nodes to use when smoothing the stellar radius (ignored
		///if smooth_radius is NaN - no smoothing).
		int radius_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///convective zone (ignored if smooth_conv_inertia is NaN - no
		///smoothing).
		int conv_inertia_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///radiative zone (ignored if smooth_rad_inertia is NaN - no
		///smoothing).
		int rad_inertia_nodes,

		///How many nodes to use when smoothing the mass of the radiative
		///zone (ignored if smooth_rad_inertia is NaN - no smoothing).
		int rad_mass_nodes,

		///How many nodes to use when smoothing the radius of the radiative
		///zone (ignored if smooth_core_env_boundary is NaN - no smoothing).
		int core_env_boundary_nodes, 

		///A set of lg(luminosities) (in \f$L_\odot\f$) for each age of
		///each track. Can be omitted if lominosity interpolation is not
		///necessary.
		const std::list< std::valarray<double> > &tabulated_luminosities=
			std::list< std::valarray<double> >(),

		///How much to smooth the luminosities when fitting.
		double smooth_luminosities=NaN,

		///How many nodes to use when smoothing the luminosities (ignored if
		///smooth_luminosities is NaN - no smoothing).
		int luminosities_nodes=0,
			
		///The mass above which the stars are considered 
		///high mass in \f$M_\odot\f$
		double max_low_mass=1.075,

		///When interpolating the age of each low mass model is scaled by
		///the mass to negative this power in order to make mass
		///dependence smoother.
		double low_mass_age_scaling=2.5,

		///When interpolating the age of each high mass model is scaled 
		///by the mass to negative this power in order to make mass
		///dependence smoother.
		double high_mass_age_scaling=3.0,

		///Low mass models are included in the mass interpolation only if
		///the required age is no larger than the maximum tabulated age 
		///times this factor.
		double low_mass_extrapolate=1.01,

		///High mass models are included in the mass interpolation only
		///if the required age is no larger than the maximum tabulated 
		///age times this factor.
		double high_mass_extrapolate=1.01);

	///\brief The moment of inertia of a stellar zone in
	/// \f$M_\odot\cdot R_\odot^2\f$.
	///
	///Returns a single argument function which gives the moment of
	///inertia of the specified zone of a star of the specified mass as a
	///function of age.
	///
	///The result must be destroyed when it becomes obsolete.
	virtual const EvolvingStellarQuantity *interpolate_moment_of_inertia(
			double stellar_mass,
			StellarZone zone,
			double present_age=-1) const;

	///\brief The stellar radius in \f$R_\odot\f$.
	///
	///Returns a single argument function which gives the radius of a 
	///star of the specified mass as a function of age.
	///
	///The result must be destroyed when it becomes obsolete.
	virtual const EvolvingStellarQuantity *interpolate_radius(
			double stellar_mass,
			double present_age=-1) const;

	///\brief The stellar luminosity in \f$L_\odot\f$.
	//
	///Returns a single argument function which gives the luminosity of a 
	///star of the specified mass as a function of age.
	///
	///The result must
	///be destroyed when it becomes obsolete.
	virtual const EvolvingStellarQuantity *interpolate_luminosity(
			double stellar_mass,
			double present_age=-1) const;

	///\brief The mass of a stellar zone in \f$M_\odot\f$
	///
	///Returns a single argument function which gives the mass of
	///of the specified zone of a star of the specified mass as a
	///function of age.
	///
	///The result must be destroyed when it becomes obsolete.
	virtual const EvolvingStellarQuantity *interpolate_zone_mass(
			double stellar_mass,
			StellarZone zone) const;

	///\brief The core-envelope boundary in \f$R_\odot\f$.
	///
	///Returns a single argument function which gives the radius of the 
	///convective-radiative boundary for a star of the specified mass as
	///a function of age.
	///
	///The result must be destroyed when it becomes obsolete. 
	virtual const EvolvingStellarQuantity *interpolate_core_boundary(
			double stellar_mass,
			double present_age=-1) const;

	///The age at which the core begins to form in Gyr.
	virtual double core_formation_age() const {return core_formation;}

	///\brief Returns the mass that separates low from high mass stars in
	/// \f$M_\odot\f$.
	double get_mass_break() const {return mass_break;}

	///\brief Serializes the interpolation state to file.
	///
	///Only call this on objects initialized with the default constructor.
	virtual void save_state(
			const std::string &filename="../interp_state_data") const;

	///\brief Loads data from serialization.
	///
	///Only call this on objects NOT initialized using the default
	///constructor (otherwise it has no data to save). Serializes state to
	///file.
	///
	///Recursively saves data of YRECEvolution and every class it depends on:
	/// - StellarEvolution
	/// - InterpolatingFunctionALGLIB
	/// - OneArgumentDiffFunction,
	/// - OneArgumentFunction
	/// - spline1dinterpolant
	/// - _spline1dinterpolant_owner.
	///
	///In _spline1dinterpolant_owner, serialize() serializes everything
	///EXCEPT p_struct->x.data and p_struct->y.data, because those are just
	///copies of the original data on which the spline was based and are not
	///necessary for evaluating the spline.
	virtual void load_state(
			const std::string &filename="../interp_state_data");

	virtual ~StellarEvolution() {}
};

#endif
