#ifndef __TIDAL_DISSIPATION_H
#define __TIDAL_DISSIPATION_H

#include "DissipatingBody.h"
#include "AstronomicalConstants.h"
#include "Common.h"
#include <cmath>

/**\file
 *
 * \brief Declares the TidalDissipation class and supporting constants.
 *
 * \ingroup StellarSystem_group
 */

///\brief Isolates constants related to the tidal dissipation.
///
///\ingroup StellarSystem_group
namespace Dissipation {
	///The quantities which evolve due to tidal dissipation
	enum Quantity {
		///The rate at which energy is deposited into the body.
		///Units: \f$\frac{M_\odot R_\odot^2 rad^2}{day^2\,Gyr}\f$
		POWER, 

		///The torque exerted on the body in the x direction.
		///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$. 
		TORQUEX, 

		///The torque exerted on the body in the z direction.
		///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
		TORQUEZ,

		///Minus the rate of change of the semimajor axis in AU/Gyr.
		SEMIMAJOR_DECAY,

		///The rate of decrease of the the angle between the spin angular
		///momentum of the body and the orbital angular momentum. In units of
		/// \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
		INCLINATION_DECAY,

		///The rate at which the eccentricity decays per Gyr.
//		ECCENTRICITY_DECAY,

		///The total number of dissipation quantitise supported.
		NUM_QUANTITIES};

	///\brief All evolving quantities also have derivatives.
	enum Derivative {

		///The quantity itself, undifferentiated.
		NO_DERIV, 

		///The derivative w.r.t. age, but the only time dependence is through
		///the modified phase lag function. The full age derivative can be
		///reconstructed lated by adding the radius derivative times the rate
		///of change of the radius and similarly for the moment of inertia.
		AGE,

		///The derivative w.r.t. the radius of the body in \f$R_\odot\f$.
		RADIUS,

		///The derivative w.r.t. the moment of inertia o the body in
		/// \f$M_\odot R_\odot^2\f$.
		MOMENT_OF_INERTIA,
				
		///The derivative w.r.t. the spin angular momentum in
		/// \f$M_\odot R_\odot^2 rad/day\f$.
		SPIN_ANGMOM,
		
		///The derivative w.r.t. the semimajor axis in AU.
		SEMIMAJOR,

		///The derivative w.r.t. the eccentricity.
//		ECCENTRICITY,

		///The derivative w.r.t. the angle between the spin angular momentum
		///of the body and the orbital angular momentum.
		INCLINATION,

		///The total number of derivatives supported
		NUM_DERIVATIVES};
};

///\brief The rates of change of various quantities due to tidal dissipation.
///
///\ingroup StellarSystem_group
class TidalDissipation {
private:
	///\brief The constant coefficiients in \f$\mathcal{U}_{m,m'}\f$ of Lai
	///(2012).
	///
	///The first index is m+2 (since m starts from -2) and the second index
	///is m'/2+1 since the only allowed values are -2, 0 and 1.
	static const double __Umm_coef[][3],

				 ///\brief \f$\kappa_{m,m'}^+/\kappa_{m,m'}\f$ as a function
				 ///of \f$m=-2 \ldots 2\f$.
				 __torque_x_plus_coef[],

				 ///\brief \f$\kappa_{m,m'}^-/\kappa_{m,m'}\f$ as a function
				 ///of \f$m=-2 \ldots 2\f$.
				 __torque_x_minus_coef[];

	///The orbital frequency in rad/day.
	double __orbital_frequency,
		   
		   ///The product of the masses of the two bodies in \f$M_\odot^2\f$.
		   __mass_product,
		   
		   ///\brief The product divided by the sum of the masses of the two
		   ///bodies in \f$M_\odot\f$.
		   __reduced_mass,
		   
		   ///\brief The orbital energy in units of 
		   /// \f$\frac{M_\odot R_\odot^2 rad^2}{day^2}\f$.
		   __orbital_energy,
		   
		   ///\brief The orbital angular momentum in units of
		   /// \f$\frac{M_\odot R_\odot^2 rad}{day}\f$.
		   __orbital_angular_momentum,
		   
		   ///The semiamjor axis in units of \f$R_\odot\f$
		   __semimajor;

	///The \f$\mathcal{U}_{m,m'}\f$ quantities defined in Lai (2012).
	std::valarray< std::valarray<double> > __Umm;

	///\brief Rates of change of the various quantities and derivatives due
	///to tidal dissipation
	std::valarray<double> __dissipation_rate,

		///The spin angular momenta of the bodies
		__spin_angular_momentum,
		
		///The inclinations of the bodies.
		__inclination;

	///\brief Returns a reference to the entry __dissipation_rates 
	///corresponding to the given variable.
	inline double& rate_entry(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV)
	{return __dissipation_rate[quantity+Dissipation::NUM_QUANTITIES*(
			derivative+Dissipation::NUM_DERIVATIVES*body_index)];}

	///Computes the \f$\mathcal{U}_{m,m'}\f$ values. 
	void fill_Umm(
			///The angle between the spin angular momentum of the body and
			///the orbital angular momentum in radians.
			double inclination,
			
			///Whether to use the derivatives w.r.t. \f$\Theta\f$ instead of
			///the values
			bool deriv);

	///\brief Calculates the dimensionless x and y torques and power due to
	///tidal dissipation.
	///
	///The \f$\mathcal{U}_{m,m'}\f$ constants must already be filled.
	void calculate_torque_power(
			///The body doing the dissipating.
			const DissipatingBody &body,

			///If this values is not zero, the spin frequency of is /assumed
			//to approach the orbital frequency from below/above if
			///forcing_sign is +1/-1 instead, regardless of the currently set
			///spin frequency of the body. Result is undefined if the spin
			///frequency is not very close to the orbital frequency and yet
			///this value is non-zero.
			short forcing_sign,

			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,
			
			///Calculate the quantity itself or one of its derivatives. For
			///all derivatives except Dissipation::SEMIMAJOR,
			///Dissipation::SPIN_ANGMOM and Dissipation::MOMENT_OF_INERTIA it
			///calculates the derivative w.r.t. the specified quantity. For 
			///Dissipation::SEMIMAJOR the derivative w.r.t. the orbital 
			///frequency. For Dissipation::SPIN_ANGMOM and 
			/// Dissipation::MOMENT_OF_INERTIA it calculates the derivative 
			///w.r.t. the spin frequency.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV);

	///Calculates the value and all derivatives of the semimajor axis decay.
	void calculate_semimajor_decay(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index);

	///Calculates the value and all derivatives of the inclination decay.
	void calculate_inclination_decay(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index);

	///Calculates the value and all derivatives of the eccentricity decay.
	void calculate_eccentricity_decay(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index);

	///Calculates the given quantity from the tidal power and torque.
	void calculate_secondary_quantity(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV);

public:
	///\brief Calculates the rates of change of various quantities due to
	///tidal dissipation.
	///
	///For now only works for zero eccentricity, throws and ecception
	///otherwise.
	TidalDissipation(
			///The first dissipating body.
			const DissipatingBody &body1,

			///The second dissipating body.
			const DissipatingBody &body2,

			///The semimajor axis of the orbit in AU.
			double semimajor,

			///The eccentricity of the orbit.
			double eccentricity,
				
			///If this values is not zero, the spin frequency of body1 is
			///assumed to approach the orbital frequency from below/above if
			///forcing_sign is +1/-1 instead, regardless of the currently set
			///spin frequency of body1. Result is undefined if the spin
			///frequency is not very close to the orbital frequency and yet
			///this value is non-zero.
			short forcing_sign1=0,

			///If this values is not zero, the spin frequency of body2 is
			///assumed to approach the orbital frequency from below/above if
			///forcing_sign is +/- instead, regardless of the currently set
			///spin frequency of body2. Result is undefined if the spin
			///frequency is not very close to the orbital frequency and yet
			///this value is non-zero.
			short forcing_sign2=0);

	///\brief Retuns the rate of change of some quantity due to tidal
	///dissipation.
	double operator()(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV);
};

#endif
