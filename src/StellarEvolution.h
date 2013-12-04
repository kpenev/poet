/**\file
 *
 * \brief Defines classes needed for interpolating among stellar evolution
 * tracks.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __STELLAR_EVOLUTION_H
#define __STELLAR_EVOLUTION_H

#include "StellarZone.h"
#include "Functions.h"
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

///\brief A class that calculates derivatives with respect to an argument for
///functions of the log(argument).
///
///\ingroup StellarSystem_group
class LogArgDerivatives : public FunctionDerivatives {
private:
	///The value of the argument at which derivatives are calculated.
	double x;

	///The value (zeroth derivative).
	mutable double value;

	///\brief Was the interpolation done against the logarithm of the
	///argument and hence the derivative needs to be corrcted.
	bool correct_log_arg;

	///\brief The currently computed derivatives.
	///
	///With respect to log(arg) if correct_log_arg is true, or with respect
	///to arg if not.
	mutable std::vector<double> underlying_deriv_values,

			///Previously calculated values of the derivatives.
			///
			///These are reused if requested multiple times.
			deriv_values;

	///\brief Actually corrects for differentiating w.r.t. log(arg) instead
	///of arg.
	///
	///Uses the pre-computed array of the derivatives up to the given order
	///in underlying_deriv_values with respect to ln(x) to return the
	///order-th derivative with respect to x.
	double transform_log_arg_deriv(unsigned order) const;
protected:
	///\brief Should be overwritten to calculate the derivatives with respect
	///to either arg or log(arg) as specified on construction.
	virtual double calc_deriv(unsigned deriv_order) const =0;
public:
	///\brief Create a derivative for functions of possibly log(arg).
	///
	///The created object corrects for the fact that the underlying
	///derivative (defined by the calc_deriv method) is with respect to the
	///logarithm of the argument if arg_val is not NaN. No correction if it
	///is Nan.
	LogArgDerivatives(double arg_val=NaN) :
		x(arg_val), value(NaN), correct_log_arg(!std::isnan(arg_val)) {}

	///Returns the deriv_order-th derivative of the quantity
	double order(unsigned deriv_order=1) const;
};

///\brief Derivative class for stellar quantities which are interpolated
///both in age and mass.
///
///\ingroup StellarSystem_group
class InterpolatedDerivatives : public LogArgDerivatives {
private:
	///The mass of the star in \f$M_\odot\f$.
	double stellar_mass;

	///The age derivatives for each stellar model.
	std::valarray<const FunctionDerivatives *> *interp_deriv;

	///The masses of the stelar models in \f$M_\odot\f$
	std::valarray<double> *interp_masses;

	///Whether to delete the derivatives and masses it was created with
	bool delete_interp_data;
protected:
	///Returns the deriv_order-th derivative of the quantity
	double calc_deriv(unsigned deriv_order) const;
public:
	///\brief Create an object that interpolates derivatives from evolution
	///tracks.
	///
	///If age is specified the input derivatives are assumed to be with
	///respect to ln(age), while derivatives always with respect to age are
	///output.
	InterpolatedDerivatives(double mass,
			std::valarray<const FunctionDerivatives*> *derivatives,
			std::valarray<double> *masses, double age=NaN,
			bool delete_data=false) :
		LogArgDerivatives(age), stellar_mass(mass),
		interp_deriv(derivatives), interp_masses(masses),
		delete_interp_data(delete_data) {}

	///Deletes the interpolation data if so specified on creation.
	~InterpolatedDerivatives()
	{
		if(delete_interp_data) {
			for(size_t i=0; i<interp_deriv->size(); i++) 
				delete (*interp_deriv)[i];
			delete interp_deriv; 
			delete interp_masses;
		}
	}
};

///\brief Makes a derivative with respect to linear argument from a
///derivative with respect to log(argument).
class RemoveLogDeriv : public LogArgDerivatives {
private:
	///The original logarithmic derivative 
	const FunctionDerivatives *__log_deriv;

	///Whether to delete the underlying log-derivative on destruction.
	bool __delete_deriv;
protected:
	///Returns the deriv_order-th derivative of the quantity
	double calc_deriv(unsigned deriv_order) const
	{return __log_deriv->order(deriv_order);}
public:
	///Create a linear derivative from a log one.
	RemoveLogDeriv(double age, const FunctionDerivatives *log_deriv,
			bool delete_deriv) :
		LogArgDerivatives(age), __log_deriv(log_deriv), 
		__delete_deriv(delete_deriv) {}

	///Deletes the input logarithmic derivative if so specified on creation.
	~RemoveLogDeriv()
	{if(__delete_deriv) delete __log_deriv;}
};

///\brief Derivative class for stellar quantities which age scaled
///quantities for a tabulated mass.
///
///\ingroup StellarSystem_group
class ScaledDerivatives : public LogArgDerivatives {
private:
	///The derivative of the underlying quantity
	const FunctionDerivatives *underlying_deriv;

	///The scaling applied to the age of the underlying quantity
	double scaling;

	///Whether to delete the input derivative when the object is destroyed.
	bool delete_underlying,
		 
		 ///\brief Was the interpolation done against log(argument) and hence
		 ///the derivative needs to be corrcted.
		 correct_log_arg;
protected:
	///\brief Returns the deriv_order-th derivative of the quantity.
	///
	///It should return the derivative with respect to either log(arg) or arg
	///as specified on construction.
	double calc_deriv(unsigned deriv_order) const;
public:
	///\brief Construct a derivative from the derivative of a function whose
	///argument is scaled by the given factor.
	///
	///If age is given and not NaN,
	///the underlying derivatives are assumed to be with respect to ln(age)
	///and a correction is made to return derivatives with respect to age.
	ScaledDerivatives(const FunctionDerivatives *deriv, double factor,
			double age=NaN, bool delete_deriv=false) :
		LogArgDerivatives(age), underlying_deriv(deriv), scaling(factor),
		delete_underlying(delete_deriv), correct_log_arg(!std::isnan(age)) {}

	///Deletes the derivative if so specified on creation.
	~ScaledDerivatives() {
		if(delete_underlying) delete underlying_deriv;
	}

};

///\brief Derivative class for a quantity that is the sum of two other
///quantities.
///
///\ingroup StellarSystem_group
class SumDerivatives : public FunctionDerivatives {
private:
	///The derivatives of the first quantity in the sum.
	const FunctionDerivatives *q1_deriv,

		  ///The derivatives of the second quantity in the sum.  
		  *q2_deriv;

	///Whether to delete the input derivative when the object is destroyed.
	bool destroy_derivs;
public:
	///Create a derivative object for a sum of two quantities: q1+q2.
	SumDerivatives(
			///Pointer to the derivative of the first quantity (q1).
			const FunctionDerivatives *derivative1,

			///Pointer to the derivative of the second quantity (q2).
			const FunctionDerivatives *derivative2,
			
			///Delete the input derivatives on destruction?
			bool delete_inputs=true)
		: q1_deriv(derivative1), q2_deriv(derivative2),
		destroy_derivs(delete_inputs) {}
	
	///The deriv_order-th derivative.
	double order(unsigned deriv_order=1) const
	{return q1_deriv->order(deriv_order)+q2_deriv->order(deriv_order);}

	///Clean up.
	~SumDerivatives()
	{if(destroy_derivs) {delete q1_deriv; delete q2_deriv;}}
};

///\brief The derivatives of an identically zero quantity.
///
///\ingroup StellarSystem_group
class ZeroDerivatives : public FunctionDerivatives {
public:
	///Create a derivative of an identically zero quantity.
	ZeroDerivatives() {}

	///The deriv_order-th derivative.
	double order(unsigned =1) const {return 0;}
};

///\brief A class for stellar properties that depend on age.
///
///\ingroup StellarSystem_group
class EvolvingStellarQuantity : public OneArgumentDiffFunction {
private:
	///The minimum age for which this quantity is defined in Gyr.
	double min_age, 

		   ///The maximum age for which this quantity is defined in Gyr.
		   max_age; 

	///Whether the tracks have log(age) instead of age as their argument.
	bool use_log_age,
		 
		 ///Should the quantity be assumed zero below the minimum track age.
		 initially_zero;

	///The masses of the evolution tracks below the high low mass split 
	std::list<double> low_masses;

	///The tracks for the model masses below the high--low mass cut
	std::list<const OneArgumentDiffFunction *> low_mass_tracks;

	///The evolution track for the high mass model with the closest mass.
	const OneArgumentDiffFunction *closest_high_mass_track;

	///The mass to which to interpolate in \f$M_\odot\f$.
	double stellar_mass,

		   ///The mass of the high mass track closest to the stellar mass in
		   /// \f$M_\odot\f$
	       closest_high_mass,

		   ///\brief Age of low mass tracks is scaled by this power in order
		   ///to make mass dependence smoother.
		   age_scaling_low_mass,


		   ///\brief Age of high mass tracks is scaled by this power in order
		   ///to make mass dependence smoother.
		   age_scaling_high_mass,

		   ///\brief How far to extrapolate low mass models.
		   ///
		   ///Low mass models are included in the mass interpolation only if
		   ///the required age is no larger than the maximum tabulated age 
		   ///times this factor.
		   extrapolate_low_mass,

		   ///\brief How far to extrapolate low mass models.
		   ///
		   ///High mass models are included in the mass interpolation only
		   ///if the required age is no larger than the maximum tabulated 
		   ///age times this factor.
		   extrapolate_high_mass;

	
	///Whether the star is low mass (as opposed to high mass)
	bool is_low_mass;

	///Prepares the class to interpolate to a mass above the high--low cutoff
	void init_high_mass(
			///The masses for which evolution tracks are given
			const std::valarray<double> &masses_of_tracks,

			///The evolution tracks of the relevant quantity
			const std::list<const OneArgumentDiffFunction *> 
			&evolution_tracks);

	///Prepares the class to interpolate to a mass below the high--low cutoff
	void init_low_mass(
			///The masses for which evolution tracks are given.
			const std::valarray<double> &masses_of_tracks,
			////The evolution tracks of the relevant quantity.
			const std::list<const OneArgumentDiffFunction *> 
			&evolution_tracks,
			///The mass above which the stars are considered high mass.
			double max_low_mass);

	///\brief Interpolate the quantity for the given track to the given age,
	///returning NaN if out of age range.
	///
	///If derivatives is not NULL initializes that to a pointer to a
	///derivatives at the current age structure.
	double evaluate_track(double age,
			const OneArgumentDiffFunction &track,
			const FunctionDerivatives **derivatives) const;

	///\brief Interpolate the quantity to the desired age assuming we are in
	///the high mass regime.
	///
	///If derivatives is not NULL initializes that to a pointer to a
	///derivatives at the current age structure.
	double high_mass_interp(double age, 
			const ScaledDerivatives **derivatives=NULL) const;

	///\brief Interpolate the quantity to the desired age assuming we are in
	///the low mass regime.
	///
	///If derivatives is not NULL initializes that to a pointer to a
	///derivatives at the current age structure.
	double low_mass_interp(double age, 
			const FunctionDerivatives **derivatives=NULL) const;

	///Calculates the derivatives of the quantity for a high mass star.
	const CubicSplineDerivatives *high_mass_deriv(double age) const;

	///Calculates the derivatives of the quantity for a low mass star.
	const InterpolatedDerivatives *low_mass_deriv(double age) const;
public:
	///\brief Construct an object that can be set to interpolate between
	///tabulated evolution tracks of a quantity.
	EvolvingStellarQuantity() {};

	///Create an evolving quantity that interpolates to the given mass.
	EvolvingStellarQuantity(
			///The stellar mass to interpolate to in \f$M_\odot\f$
			double mass, 

			///The masses for which evolution tracks are given in
			/// \f$M_\odot\f$
			const std::valarray<double> &masses_of_tracks,

			///The evolution tracks of the relevant quantity
			const std::list<const OneArgumentDiffFunction *> 
			&evolution_tracks,

			///Whether the track uses log(age) as the independent argument
			///instead of age.
			bool log_age=true,

			///The mass above which the stars are considered high mass in
			/// \f$M_\odot\f$.
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
			double high_mass_extrapolate=1.01,
			
			///Whether this is a quantity that is identically zero below some
			///age and turns on afterwards
			bool starts_zero=false);


	///Return the value the quantity takes at the given age in Gyr.
	virtual double operator()(double age) const;

	///Return the age derivative of the quantity at the given age in Gyr.
	virtual const FunctionDerivatives *deriv(double age) const;

	///The largest age for which the quantity can be interpolated in Gyr.
	virtual double range_high() const {return max_age;}

	///The smallest age for which the quantity can be interpolated in Gyr.
	virtual double range_low() const {return min_age;}

	///An iterator over the ages (in Gyr) where the quantity takes the given
	///y value.
	InterpSolutionIterator crossings(double =0) const
	{throw Error::Runtime("Called EvolvingStellarQuantity::crossings, "
			"which are ill defined.");}
};

class ZeroQuantity : public EvolvingStellarQuantity {
public:
	///Return the value the quantity takes at the given age.
	double operator()(double) const {return 0;}

	///Return the age derivative of the quantity at the given age.
	const FunctionDerivatives *deriv(double) const
	{return new ZeroDerivatives;}

	///The largest age for which the quantity can be interpolated
	double range_high() const {return Inf;}

	///The smallest age for which the quantity can be interpolated.
	double range_low() const {return -Inf;}

	///An iterator over the ages where the quantity takes the given y value.
	InterpSolutionIterator crossings(double =0) const
	{throw Error::Runtime("Called ZeroQuantity::crossings, "
			"which are ill defined.");}
};

///\brief A clas for stellar quantities that are the sum of two other
///quantities.
///
///\ingroup StellarSystem_group
class SumQuantity : public EvolvingStellarQuantity {
private:
	///This quantity will be q1+q2
	const EvolvingStellarQuantity *q1, *q2;

	///Whether to destroy the input quantities on destruction
	bool destroy_qs;
public:
	///Create a quantity that is (*quantity1)-(*quantity2)
	SumQuantity(const EvolvingStellarQuantity *quantity1,
			const EvolvingStellarQuantity *quantity2,
			bool delete_inputs=false)
		: q1(quantity1), q2(quantity2), destroy_qs(delete_inputs) {}

	///Return the value the quantity takes at the given age.
	double operator()(double age) const
	{return (*q1)(age)+(*q2)(age);}

	///Return the age derivative of the quantity at the given age.
	const FunctionDerivatives *deriv(double age) const
	{return new SumDerivatives(q1->deriv(age), q2->deriv(age), true);}

	///The largest age for which the quantity can be interpolated
	double range_high() const
	{return std::min(q1->range_high(), q2->range_high());}

	///The smallest age for which the quantity can be interpolated.
	double range_low() const
	{return std::max(q1->range_low(), q2->range_low());}

	///An iterator over the ages where the quantity takes the given y value.
	InterpSolutionIterator crossings(double =0) const
	{throw Error::Runtime("Called EvolvingStellarQuantity::crossings, "
			"which are ill defined.");}

	///Clean up.
	~SumQuantity()
	{if(destroy_qs) {delete q1; delete q2;}}
};

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
		ar & track_masses;
		ar & interpolated_radius & interpolated_luminosity &
			interpolated_conv_inertia & interpolated_rad_inertia &
			interpolated_rad_mass & interpolated_core_env_boundary;
		ar & mass_break & low_age_scaling & high_age_scaling &
			extrapolate_low & extrapolate_high & core_formation;
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
		///if #smooth_radius is NaN - no smoothing).
		int radius_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///convective zone (ignored if #smooth_conv_inertia is NaN - no
		///smoothing).
		int conv_inertia_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///radiative zone (ignored if #smooth_rad_inertia is NaN - no
		///smoothing).
		int rad_inertia_nodes,

		///How many nodes to use when smoothing the mass of the radiative
		///zone (ignored if #smooth_rad_inertia is NaN - no smoothing).
		int rad_mass_nodes,

		///How many nodes to use when smoothing the radius of the radiative
		///zone (ignored if #smooth_core_env_boundary is NaN - no smoothing).
		int core_env_boundary_nodes, 

		///A set of lg(luminosities) (in \f$L_\odot\f$) for each age of
		///each track. Can be omitted if lominosity interpolation is not
		///necessary.
		const std::list< std::valarray<double> > &tabulated_luminosities=
			std::list< std::valarray<double> >(),

		///How much to smooth the luminosities when fitting.
		double smooth_luminosities=NaN,

		///How many nodes to use when smoothing the luminosities (ignored if
		///#smooth_luminosities is NaN - no smoothing).
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
		///if #smooth_radius is NaN - no smoothing).
		int radius_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///convective zone (ignored if #smooth_conv_inertia is NaN - no
		///smoothing).
		int conv_inertia_nodes,

		///How many nodes to use when smoothing the moment of inertia of the
		///radiative zone (ignored if #smooth_rad_inertia is NaN - no
		///smoothing).
		int rad_inertia_nodes,

		///How many nodes to use when smoothing the mass of the radiative
		///zone (ignored if #smooth_rad_inertia is NaN - no smoothing).
		int rad_mass_nodes,

		///How many nodes to use when smoothing the radius of the radiative
		///zone (ignored if #smooth_core_env_boundary is NaN - no smoothing).
		int core_env_boundary_nodes, 

		///A set of lg(luminosities) (in \f$L_\odot\f$) for each age of
		///each track. Can be omitted if lominosity interpolation is not
		///necessary.
		const std::list< std::valarray<double> > &tabulated_luminosities=
			std::list< std::valarray<double> >(),

		///How much to smooth the luminosities when fitting.
		double smooth_luminosities=NaN,

		///How many nodes to use when smoothing the luminosities (ignored if
		///#smooth_luminosities is NaN - no smoothing).
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
