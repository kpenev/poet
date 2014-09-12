#ifndef __DISSIPATING_BODY_H
#define __DISSIPATING_BODY_H

#include "SpinOrbitLockInfo.h"
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
