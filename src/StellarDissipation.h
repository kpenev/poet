#ifndef __STELLAR_Q_H
#define __STELLAR_Q_H

#include "Star.h"

/**\file
 *
 * \brief Defines the tidal dissipation properties of stars.
 *
 * \ingroup StellarSystem_group
 */

class Star : public StarBase {
private:
	///The tidal quality factor of the star.
	double __mod_phase_lag;

public:
	///Create a star with the given properties.
	Star(
			///Mass of the star
			double mass,

			///Tidal quality factor
			double tidal_quality, 

			///The normalization constant (K) in the magnetic wind equation
			//in \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
			/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
			double wind_strength,

			///The saturation frequency (rad/day) in the magnetic wind
			///equation.
			double wind_saturation,

			///The timescale on which the rotation of the core and the
			///envelope are coupled in Gyr.
			double coupling_timescale,
			
			///The frequency range (in rad/day) over which 1/Q decays to zero.
			///Deprecated: This parameter was initially used to avoid a
			///discontinuity in the equations when the orbit and the stallar
			///spin went through synchroneity. It can now be safely set to
			///zero.
			double,

			///The frequency of the stellar convective zone while the
			///disk is present in rad/day.
			double disk_lock_ang_vel,
			
			///The age (in Gyrs) when the circumstellar disk dissipates
			///and no longer keeps the surface rotating at a constant angular
			///velocity.
			double disk_lock_time,

			///A StellarEvolution interpolator.
			const StellarEvolution &evolution,

			///The present age of the star in Gyr (optional).
			double age=NaN,
			
			///The present surface rotation rate of the star in rad/day
			///(optional).
			double conv_spin=NaN,

			///The present core rotation rate of the star in rad/day
			///optional).
			double rad_spin=NaN) :
				StarBase(mass, wind_strength, wind_saturation,
						coupling_timescale, disk_lock_ang_vel,
						disk_lock_time, evolution, age,
						conv_spin, rad_spin),
				__mod_phase_lag(30.0/(16.0*M_PI*tidal_quality)) {}

	///\brief A function defining the dissipation efficiency of the body.
	///
	///In our formalism this function should return
	/// \f$\Delta_{m,m'}'\equiv\kappa_{m,m'}\sin(\Delta_{m,m'})\f$.
	///For small dissipation efficiency, this is the phase lag times the love
	///number. It must satisfy \f$\Delta_{m,m'}=\Delta_{-m,-m'}\f$.
	double modified_phase_lag(
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
			PhaseLag::Derivative derivative=PhaseLag::NO_DERIV) const;
};
#endif
