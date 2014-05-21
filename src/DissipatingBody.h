#ifndef __DISSIPATING_BODY_H
#define __DISSIPATING_BODY_H

#include "Common.h"

#include <assert.h>

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
		NO_DERIV,

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
	///The mass of the body in \f$M_\odot\f$.
	double __mass,

		   ///The radius of the body in \f$R_\odot\f$.
		   __radius,

		   ///The moment of inertia of the body in \f$M_\odot R_\odot^2\f$.
		   __moment_of_inertia,

		   ///The spin frequency of the body in rad/day.
		   __spin,

		   ///Angle between the spin and orbital angular momentum in radians.
		   __inclination;
public:
	///Create a dissipating body with the given parameters.
	DissipatingBody(
			///The current mass of the body in \f$M_\odot\f$.
			double mass=NaN,

			///The current radius of the body in \f$R_\odot\f$.
			double radius=NaN,
			
			///The current spin frequency of the body in rad/day.
			double spin=NaN,

			///The current angle between the angular momentum of the body and
			///the orbit
			double inclination=NaN) :
		__mass(mass), __radius(radius), __spin(spin),
		__inclination(inclination)
	{
#ifdef DEBUG
		assert(!(mass<0));
		assert(!(radius<0));
#endif
	}

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
			///forcing frequency is encountered, this flag determines if the
			///zero should be interpreted as an infinitesimal positive or
			///negative amount above. It is safe to ignore this flag if the
			///function is continuous at zero forcing frequency.
			short forcing_sign,

			///Whether to return the phase lag or one of its derivatives.
			PhaseLag::Derivative derivativ=PhaseLag::NO_DERIV) const=0;

	///The current radius of this body in \f$R_\odot\f$.
	virtual double current_radius() const 
	{
#ifdef DEBUG
		assert(!std::isnan(__radius));
#endif
		return __radius;
	}

	///\brief Set the current radius of this body to the given value in
	/// \f$R_\odot\f$.
	virtual void current_radius(double radius) 
	{
#ifdef DEBUG
		assert(radius>0);
#endif
		__radius=radius;
	}

	///The current mass of this body in \f$M_\odot\f$.
	virtual double current_mass() const
	{
#ifdef DEBUG
		assert(!std::isnan(__mass));
#endif
		return __mass;
	}

	///Set the current mass of this body to the given value in \f$M_\odot\f$.
	virtual void current_mass(double mass)
	{
#ifdef DEBUG
		assert(mass>0);
#endif
		__mass=mass;
	}

	///The current moment of inertia of the body in \f$M_\odot R_\odot^2\f$.
	virtual double current_moment_of_inertia() const
	{
#ifdef DEBUG
		assert(!std::isnan(__moment_of_inertia));
#endif
		return __moment_of_inertia;
	}

	///\brief Set the current moment of inertia to the givan value in
	/// \f$M_\odot R_\odot^2\f$.
	virtual void current_moment_of_inertia(double moment_of_inertia)
	{
#ifdef DEBUG
		assert(moment_of_inertia>0);
#endif
		__moment_of_inertia=moment_of_inertia;
	}

	///The current spin frequency of the body in rad/day.
	virtual double current_spin() const
	{
#ifdef DEBUG
		assert(!std::isnan(__spin));
#endif
		return __spin;
	};

	///\brief Set the current spin frequency of the body to the given value
	///in rad/day.
	virtual void current_spin(double spin) {__spin=spin;}

	///\brief The current angle between the spin and orbital angular momentum
	///in radians.
	virtual double current_inclination() const
	{
#ifdef DEBUG
		assert(!std::isnan(__inclination));
#endif
		return __inclination;
	}

	///\brief Set the current angle between the spin and orbital angular 
	///momentum to the given value in radians.
	virtual void current_inclination(double inclination)
	{__inclination=inclination;}
};

#endif
