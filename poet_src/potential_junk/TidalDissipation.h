#ifndef __TIDAL_DISSIPATION_H
#define __TIDAL_DISSIPATION_H

#include "DissipationQuantities.h"
#include "DissipatingBody.h"
#include "AstronomicalConstants.h"
#include "OrbitalExpressions.h"
#include "Common.h"
#include <cmath>
#include <vector>

/**\file
 *
 * \brief Declares the TidalDissipation class and supporting constants.
 *
 * \ingroup StellarSystem_group
 */

///\brief The rates of change of various quantities due to tidal dissipation.
///
///\ingroup StellarSystem_group
class TidalDissipation {
private:
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

	///\brief Rates of change of the various quantities and derivatives due
	///to tidal dissipation for each body. 
	///
	///There are three separate entries for quantity:
	/// - the evolution rate only due to the locked terms, assuming the
	///   body's spin is approaching the lock from below.
	/// - the evolution rate due to only non-locked terms
	/// - the evolution rate only due to the locked terms, assuming the
	///   body's spin is approaching the lock from above.
	std::valarray<double> __dissipation_rate,

		///The spin angular momenta of the bodies
		__spin_angular_momentum,
		
		///The inclinations of the bodies.
		__inclination;

	///\brief The terms which participate in the evolution and whether to
	///force them to be below or above synchroneity one set of terms for each
	///body.
	std::vector< std::vector< SpinOrbitLockInfo > > __spin_orbit_harmonics;

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

	///Fills the __spin_orbit_haromincs structure with all terms locked.
	void fill_spin_orbit_harmonics();

	///\brief Returns the forcing frequency that should be used for the given
	///term.
	///
	///Returns zero if the nominal forcing frequency is on the other side of
	///what is defined by the corresponding __spin_orbit_haromincs entry.
	double forcing_frequency(
			///For which body are we constructing the term. 
			short body_index,
			
			///The coefficient in front of the orbital frequency identifying
			///the tidal term.
			int orbital_multiplier,

			///The coefficient in front of the spin frequency of the body,
			///identifying the tidal term.
			int spin_multiplier,
			
			///The spin frequency multiplied by its multiplier.
			double multiplied_spin);

	///\brief Calculates the dimensionless x and y torques and power due to
	///tidal dissipation.
	///
	///The \f$\mathcal{U}_{m,m'}\f$ constants must already be filled.
	void calculate_torque_power(
			///The body doing the dissipating.
			const DissipatingBody &body,

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
		{fill_spin_orbit_harmonics();}

	///See #init().
	TidalDissipation(
			const DissipatingBody &body1,
			const DissipatingBody &body2,
			double semimajor,
			double eccentricity,
			const SpinOrbitLockInfo &lock1=SpinOrbitLockInfo(),
			const SpinOrbitLockInfo &lock2=SpinOrbitLockInfo())
	{
		fill_spin_orbit_harmonics(); 
		init(body1, body2, semimajor, eccentricity, lock1, lock2);
	}

	///\brief Calculates the rates of change of various quantities due to
	///tidal dissipation.
	///
	///For now only works for zero eccentricity, throws and ecception
	///otherwise.
	///
	///The lock states to assume for the various terms must already be set
	///correctly by calling one of the init_haromnics methods.
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
			///identify the term that is locked. Leave at default to use the
			//currently set locks.
			const SpinOrbitLockInfo &lock1=SpinOrbitLockInfo(),

			///Whether to assume that body 2 is in a spin orbit lock and
			///identify the term that is locked. Leave at default to use the
			//currently set locks.
			const SpinOrbitLockInfo &lock2=SpinOrbitLockInfo()
			);

	///\brief Rates of change of quantities due to tidal dissipation split
	///into locked and non-locked terms.
	double operator()(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV,
			
			///Whether to return one of the locked rates or the non-locked
			///rate.
			short lock_dir=0) const;

	///\brief Rates of change of quantities due to tidal dissipation for a 
	///particular combination of locked and non-locked terms.
	double operator()(
			///Which body's dissipation is needed (shold be 0 or 1).
			short body_index,

			///The quantity needed.
			Dissipation::Quantity quantity,
			
			///The fraction of above the lock terms to use (the below the
			///lock terms are assigned complimentary fraction).
			double above_fraction,

			///Return the quantity itself or one of its derivatives.
			Dissipation::Derivative derivative=Dissipation::NO_DERIV) const;

	///\brief The semimajor axis used when creating this tidal dissipation
	///object in AU.
	double semimajor() const {return __semimajor;}

	///\brief The orbital frequency corresponding to the semimajor axis in
	///rad/day.
	double orbit_frequency() const {return __orbital_frequency;}

	///The orbital energy in \f$M_\odot R_\odot^2 rad^2/day^2\f$.
	double orbit_energy() const {return __orbital_energy;}

	///\brief The spin angular momentum of one of the bodies in 
	/// \f$M_{odot}R_{odot}^2 rad/day\f$.
	double spin_angular_momentum(short body_index) const
	{return __spin_angular_momentum[body_index];}

	///\brief Initialize __spin_orbit_harmonics assuming that the given
	///harmonic approaches zero from the given direction.
	void init_harmonics(
			///Which body's harmonics are we setting.
			short body_index,

			///The multiple of the orbital frequency which is in sync.
			int orbital_frequency_multiplier,
			
			///The multiplier of the spin frequency which is in sync.
			int spin_frequency_multiplier,
			
			///The direction from which the harmonic approaches zero. Set to
			///zero for a locked term.
			short direction)
	{init_harmonics(body_index,
			SpinOrbitLockInfo(orbital_frequency_multiplier,
				spin_frequency_multiplier, direction));}

	///\brief Initialize __spin_orbit_harmonics. 
	void init_harmonics(
			///Which body's harmonics are we setting.
			short body_index,

			///The harmonic to assume zero. If a term corresponding to the
			///given multipliers is included in the dissipation, the locked
			///state of that term is set to whatever lock's is.
			const SpinOrbitLockInfo &lock);

	///\brief Initialize __spin_orbit_harmonics assuming the given spin to
	///orbital frequency ratio. 
	///
	///The ratio must not precisely match any of the harmonics. Use
	/// #init_harmonics(int, int) if some terms should have precisely zero
	///forcing frequency.
	void init_harmonics(
			///Which body's harmonics are we setting.
			short body_index,
			
			///The spin frequency is assumed to be this number times the
			///orbital freqency and the values of the lock are set
			///accordingly.
			double spin_to_orbital_ratio);

	///\brief The number of different spin-orbit harmonics contirubuting to
	///the dissipation.
	unsigned num_harmonics() const {return 4;}

	///\brief Returns a const reference to given harmonic.
	const SpinOrbitLockInfo &harmonic(short body_index, unsigned i)
		 const {return __spin_orbit_harmonics[body_index][i];}
};

#endif
