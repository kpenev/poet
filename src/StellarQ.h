#ifndef __STELLAR_Q_H
#define __STELLAR_Q_H

#include "Star.h"

/**\file
 *
 * \brief Defines the tidal dissipation parameter \f$Q^*(\omega)\f$.
 *
 * \ingroup StellarSystem_group
 */

class Star : public StarBase {
private:
	///The tidal quality factor of the star.
	double tidal_Q,

		   ///\deprecated
		   ///\brief The frequency range (in rad/day) over which 1/Q decays
		   ///to zero.
		   ///
		   ///This parameter was initially used to avoid a
		   ///discontinuity in the equations when the orbit and the stallar
		   ///spin went through synchroneity. It can now be safely set to
		   ///zero.
		   Q_transition_width;
public:
	///Create a star with the given properties.
	Star(
			///Mass of the star
			double current_mass,

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
			double dissipation_transition_width,

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
			double current_age=NaN,
			
			///The present surface rotation rate of the star in rad/day
			///(optional).
			double current_conv_spin=NaN,

			///The present core rotation rate of the star in rad/day
			///optional).
			double current_rad_spin=NaN) :
				StarBase(current_mass, wind_strength, wind_saturation,
						coupling_timescale, disk_lock_ang_vel,
						disk_lock_time, evolution, current_age,
						current_conv_spin, current_rad_spin),
				tidal_Q(tidal_quality), 
				Q_transition_width(dissipation_transition_width) {}
						
	///\brief The tidal quality factor of the star for tides having the
	///specified frequency (in rad/day).
	double get_tidal_Q(double tidal_frequency) const;

	///\brief The frequency derivative of tidal quality factor of the star
	///for tides having the specified frequency (in rad/day).
	double get_tidal_Q_deriv(double tidal_frequency) const;

	///\brief The transition width of the tidal quality factor.
	double get_trans_width() const {return Q_transition_width;}
};
#endif
