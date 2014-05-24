#ifndef __TIDAL_DISSIPATION_H
#define __TIDAL_DISSIPATION_H

#include "DissipatingBody.h"
#include "AstronomicalConstants.h"
#include "OrbitalExpressions.h"
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

		///The rate of change of the orbital frequency in rad/day/Gyr.
		ORBIT_SPINUP,

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

///More civilized output for Dissipation::Quantity variables.
std::ostream &operator<<(std::ostream &os,
		const Dissipation::Quantity &quantity);

///More civilized output for Dissipation::Derivative variables.
std::ostream &operator<<(std::ostream &os,
		const Dissipation::Derivative &deriv);

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
	///to tidal dissipation for each body. 
	///
	///There are three separate entries for quantity:
	/// - the evolution rate only due to the locked terms, assuming the
	//    body's spin is approaching the lock from below.
	//  - the evolution rate due to only non-locked terms
	//  - the evolution rate only due to the locked terms, assuming the
	//    body's spin is approaching the lock from above.
	std::valarray<double> __dissipation_rate,

		///The spin angular momenta of the bodies
		__spin_angular_momentum,
		
		///The inclinations of the bodies.
		__inclination;

#ifdef DEBUG
	///Is this a valid tidal dissipation?
	bool __valid;
#endif

	///\brief Returns a reference to the entry __dissipation_rates 
	///corresponding to the given variable.
	inline double& rate_entry(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative,
			
			///Which lock state terms to access (+/-1 for body spin 
			///approaching the lock from above/below, 0 for not locked terms.
			short lock_direction)
	{
		return __dissipation_rate[
			lock_direction+1 + 3*(
				quantity + Dissipation::NUM_QUANTITIES*(
					derivative + Dissipation::NUM_DERIVATIVES*body_index
				)
			)
		];
	}

	///Returns a copy of the same value that the non-const rate entry does.
	inline double rate_entry(
			short body_index,
			Dissipation::Quantity quantity,
			Dissipation::Derivative derivative,
			short lock_direction) const
	{
		return __dissipation_rate[
			lock_direction+1 + 3*(
				quantity + Dissipation::NUM_QUANTITIES*(
					derivative + Dissipation::NUM_DERIVATIVES*body_index
				)
			)
		];
	}

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

			///Whether to assume that body 1 is in a spin orbit lock and
			///identify the term that is locked.
			SpinOrbitLockInfo &lock,

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

public:
	///Create an invalid tidal dissipation, which fails an assert if used.
	TidalDissipation()
#ifdef DEBUG
		: __valid(false)
#endif
		{}

	///See #init().
	TidalDissipation(
			const DissipatingBody &body1,
			const DissipatingBody &body2,
			double semimajor,
			double eccentricity,
			SpinOrbitLockInfo &lock1,
			SpinOrbitLockInfo &lock2)
	{init(body1, body2, semimajor, eccentricity, lock1, lock2);}

	///\brief Calculates the rates of change of various quantities due to
	///tidal dissipation.
	///
	///For now only works for zero eccentricity, throws and ecception
	///otherwise.
	///
	///If some term in the tidal dissipation equations is locked (see the
	///lock1 and lock2 arguments) the dissipation due to that term is kept
	///separate.
	void init(
			///The first dissipating body.
			const DissipatingBody &body1,

			///The second dissipating body.
			const DissipatingBody &body2,

			///The semimajor axis of the orbit in AU.
			double semimajor,

			///The eccentricity of the orbit.
			double eccentricity,
				
			///Whether to assume that body 1 is in a spin orbit lock and
			///identify the term that is locked. The direction is ignored. If
			///this variable converts to true, the angular momentum and spin 
			///frequency of body1 are never used.
			SpinOrbitLockInfo &lock1,

			///Whether to assume that body 2 is in a spin orbit lock and
			///identify the term that is locked. The direction is ignored. If
			///this variable converts to true, the angular momentum and spin 
			///frequency of body2 are never used.
			SpinOrbitLockInfo &lock2);

	///\brief Retuns the rate of change of some quantity due to tidal
	///dissipation.
	double operator()(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV,
			
			///Whether to return one of the locked rates or the non-locked
			///rate.
			short lock_dir=false) const;

	///\brief The semimajor axis used when creating this tidal dissipation
	///object in AU.
	double semimajor() const {return __semimajor;}

	///\brief The orbital frequency corresponding to the semimajor axis in
	///rad/day.
	double orbital_frequency() const {return __orbital_frequency;}

	///\brief The spin frequency of one of the bodies in rad/day.
	double spin_angular_momentum(short body_index) const
	{return __spin_angular_momentum[body_index];}
};

#endif
