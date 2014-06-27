#ifndef __DISSIPATING_BODY_H
#define __DISSIPATING_BODY_H

#include "Common.h"

#include <cassert>

/**\file 
 *
 * \brief Declares the DissipatingBody class and supporting constants.
 *
 * \ingroup StellarSystem_group
 */


///\brief Isolate the constant identifying what to differentiate phase lags
///w.r.t.
///
///\ingroup StellarSystem_group
namespace PhaseLag {
	///What derivatives should a dissipating body be able to compute.
	enum Derivative {
		///No derivative, calculate the phase lag itself.
		NO_DERIV=0,

		///Derivative w.r.t. age (due to the body evolving) per Gyr.
		AGE,

		///Derivative w.r.t. the spin frequency in units of rad/day.
		SPIN_FREQUENCY,

		///Derivative w.r.t. the forcing frequency in units of rad/day.
		FORCING_FREQUENCY
	};
};

///\brief Defines a lock between the spin of a dissipating body and the
///orbit.
///
///With inclined and eccentric orbits, locks can occur at many different
//frequencies, not only when the orbital and spin periods are the same. In
///general almost any rational ratio can result in a lock if the dissipation
///has the appropriate frequency dependence.
///
///\ingroup StellarSystem_group
class SpinOrbitLockInfo {
private:
	///The mutiplier in front of the orbital frequency in the lock.
	int __orbital_freq_mult,

		///The multiplier in front of the spin frequency in the lock.
		__spin_freq_mult;

	///\brief Should a lock be assumed, and if so from which direction is it
	///approached?
	///
	///The values have the following meanings:
	/// - <0 	: the spin frequency of the body is slightly lagrer than
	///			  necessary for a lock (or the orbital freqency is slightly
	///		      smaller).
	///
	/// - 0 	: The spin frequency and the orbital frequency have precisely
	///			  the values to result in zero forcing frequency for this
	///			  term.
	///
	/// - >0 	: the spin frequency of the body is slightly smaller than
	///			  necessary for a lock (or the orbital freqency is slightly
	///		      larger).
	short __lock_direction;

public:
	///\brief Define which tidal dissipation term is in a lock.
	SpinOrbitLockInfo(
			///The multiple of the orbital frequency at the lock.
			int orbital_freq_mult=0,

			///The multiple of the spin frequency at the lock.
			int spin_freq_mult=0,

			///The direction from which the spin frequency is approaching the
			///lock. See #__lock_direction for the meaning of the values.
			short lock_direction=-1)
	{set_lock(orbital_freq_mult, spin_freq_mult, lock_direction);}

	///\brief Copy the original to this.
	SpinOrbitLockInfo(const SpinOrbitLockInfo &orig)
	{set_lock(orig.orbital_frequency_multiplier(),
			orig.spin_frequency_multiplier(),
			orig.lock_direction());}

	///\brief Make this the same lock as RHS.
	SpinOrbitLockInfo &operator=(const SpinOrbitLockInfo &rhs)
	{set_lock(rhs.orbital_frequency_multiplier(),
			rhs.spin_frequency_multiplier(),
			rhs.lock_direction()); return *this;}

	///\brief Define which tidal dissipation term is in a lock.
	void set_lock(
			///The multiple of the orbital frequency at the lock.
			int orbital_freq_mult,

			///The multiple of the spin frequency at the lock.
			int spin_freq_mult,

			///The sign of the forcing frequency for this term.
			///See #__lock_direction for the meaning of the values.
			short lock_direction);

	///\brief Spin frequency at exactly the lock that corresponds to the
	///given orbital frequency.
	double spin(double orbital_frequency) const
	{return (orbital_frequency*__orbital_freq_mult)/__spin_freq_mult;}

	///Is the given tidal dissipation term one of the locked terms?
	bool operator()(
			///The multiple of the orbital frequency to consider.
			int orbital_freq_mult,

			///The multiple of the spin frequency to consider.
			int spin_freq_mult) const
	{return (__lock_direction ? false :
			 term(orbital_freq_mult, spin_freq_mult));}

	///Returns true if the lock is referring to the given term, regardless of
	///whether it is locked or not.
	bool term(			
			///The multiple of the orbital frequency to consider.
			int orbital_freq_mult,

			///The multiple of the spin frequency to consider.
			int spin_freq_mult) const
	{return orbital_freq_mult*__spin_freq_mult==
		spin_freq_mult*__orbital_freq_mult;}

	///Should this lock be assumed.
	operator bool() const {return __lock_direction==0;}

	///The sign of the forcing frequnecy associated with this component.
	short lock_direction() const {return __lock_direction;}

	///Set the lock direction to the given value.
	void lock_direction(short value) {__lock_direction=value;}

	///The multiplier in front of the orbital frequency in the lock.
	int orbital_frequency_multiplier() const {return __orbital_freq_mult;}

	///The multiplier in front of the spin frequency in the lock.
	int spin_frequency_multiplier() const {return __spin_freq_mult;}

	///Are the two locks for the same frequency ratio and in the same
	///enabled/disabled state.
	bool operator==(const SpinOrbitLockInfo &rhs) const;
};

///Civilized output for locks.
std::ostream &operator<<(std::ostream &os, const SpinOrbitLockInfo &lock);

///\brief A base class for any body contributing to tidal dissipation.
///
///\ingroup StellarSystem_group
class DissipatingBody {
private:
	///\brief The spin angular momentum of the body in 
	/// \f$M_\odot R_\odot^2 rad/day\f$.
	double __angular_momentum,

		   ///Angle between the spin and orbital angular momentum in radians.
		   __inclination;
public:
	///\brief Initialize the angular momentum and the inclination of the
	///DissipatingBody.
	DissipatingBody(double angular_momentum=NaN, double inclination=NaN) :
		__angular_momentum(angular_momentum), __inclination(inclination) {}

	///\brief A function defining the dissipation efficiency of the body.
	///
	///In our formalism this function should return
	/// \f$\Delta_{m,m'}'\equiv\kappa_{m,m'}\sin(\Delta_{m,m'})\f$.
	///For small dissipation efficiency, this is the phase lag times the love
	///number. It must satisfy \f$\Delta_{m,m'}=\Delta_{-m,-m'}\f$.
	virtual double modified_phase_lag(
			///The \f$m\f$ index.
			int m, 

			///The forcing frequency (\f$m'\Omega-m\Omega_s\f$ in Lai, 2012)
			///in rad/day.
			double forcing_frequency,
			
			///If this function is discontinuous at zero and an exactly zero
			///forcing frequency is encountered, this argument determines if
			///the zero should be interpreted as an infinitesimal positive or
			///negative amount. If this argument evaluates to true, the
			///pre-set spin should be ignored and the spin at the lock used
			///instead.
			const SpinOrbitLockInfo &lock,

			///Whether to return the phase lag or one of its derivatives.
			PhaseLag::Derivative derivative=PhaseLag::NO_DERIV) const=0;

	///The radius of this body in \f$R_\odot\f$.
	virtual double radius() const =0;

	///The mass of this body in \f$M_\odot\f$.
	virtual double mass() const =0;

	///The moment of inertia of the body in \f$M_\odot R_\odot^2\f$.
	virtual double moment_of_inertia() const =0;

	///The angular momentum of the body in \f$M_\odot R_\odot^2 rad/day\f$.
	virtual double angular_momentum() const
	{
#ifdef DEBUG
		assert(!std::isnan(__angular_momentum));
#endif
		return __angular_momentum;
	};

	///\brief Set the angular momentum of the body to the given value
	///in \f$M_\odot R_\odot^2 rad/day\f$.
	virtual void angular_momentum(double value)
	{__angular_momentum=value;}

	///\brief The spin angular velocity of the body in rad/day.
	virtual double spin() const
	{
#ifdef DEBUG
		assert(!std::isnan(__angular_momentum));
#endif
		return __angular_momentum/moment_of_inertia();
	}

	///The angle between the spin and orbital angular momentum in radians.
	virtual double inclination() const
	{
#ifdef DEBUG
		assert(!std::isnan(__inclination));
#endif
		return __inclination;
	}

	///\brief Set the angle between the spin and orbital angular momentum to
	///the given value in radians.
	virtual void inclination(double value)
	{__inclination=value;}
};

#endif
