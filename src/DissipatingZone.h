#ifndef __DISSIPATING_ZONE_H

/**\file
 *
 * \brief Declares a class representing one zone of a body dissipative to
 * tidal distortions.
 */

#include "EccentricityExpansionCoefficients.h"
#include "DissipationQuantities.h"
#include "Common.h"
#include <valarray>

///\brief A layer of a system body for which the tidal bulge is not exactly
///in phase with the tidal potential.
class DissipatingZone {
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

	///The eccentricity expansion of \f$p_{m,s}\f$.
	static EccentricityExpansionCoefficients __pms;

	///The inclination with which __Ummp was last filled.
	double __Ummp_inclination,
		   
		   ///The spin frequency of the zone.
		   __spin_frequency,
		   
		   ///The absolute value of the angular momentum of the orbit
		   __orbital_angmom;

	///The \f$\mathcal{U}_{m,m'}\f$ quantities defined in Lai (2012).
	std::valarray< std::valarray<double> > __Ummp,

		///The derivatives of the \f$\mathcal{U}_{m,m'}\f$ quantities w.r.t.
		///the inclination.
		__Ummp_deriv;

	///\brief The dimensionless tidal power and its derivatives.
	///
	///Consists of pairs of numbers one for each derivative. The first number
	///of each pair is always filled and if the zone is in a lock it is the
	///tidal power calculated assuming the zone spin frequency approaches the
	///lock from below. The second number is filled only if the zone is in a
	///spin-orbit lock and is the tidal power assuming the zone spin
	///frequency approaches the lock from above.
	std::valarray<double> __power,

		///\brief The dimensionless tidal torque in the x direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_x,

		///\brief The dimensionless tidal torque in the y direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_y,

		///\brief The dimensionless tidal torque in the y direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_y,

		///\brief The dimensionless tidal torque in the z direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_z;

	///Computes the \f$\mathcal{U}_{m,m'}\f$ values and their derivatives. 
	void fill_Umm();

	///\brief Calculates \f$\sum_s W_{2,s}D_{m,s}(\Theta)p_{s,m'}\f$ (see
	///documentation) and its derivatives w.r.t. e and \f$\Theta\f$.
	///
	///fill_Umm should already have been called with the appropriate
	///inclination.
	void potential_term(
			///The eccentricity.
			double e,

			///The m index.
			int m,

			///The m' index.
			int mp,
			
			///The highest order of the eccentricity to include.
			unsigned e_order,

			///Set to the undifferentiated value.
			double &no_deriv,

			///Set to the inclination derivative.
			double &inclination_deriv,

			///Set to the eccentricity_derivative.
			double &eccentricity_deriv);

protected:
	///To what order should eccentricity expansion be performed for the given
	///value of the eccentricity.
	virtual unsigned eccentricity_order(double e) const =0;

public:
	DissipatingZone();

	///\brief Defines the current orbit, triggering re-calculation of all
	///quantities.
	void set_orbit(
			///The angular velocity of the orbit in rad/day.
			double orbital_frequency,

			///The eccentricity of the orbit
			double eccentricity,
			
			///The absolute value of the angular momentum of the orbit.
			double orbital_angmom);


	///The angle between the angular momenta of the zone and the orbit.
	virtual double inclination() const =0;

	///The argument of periapsis of this zone minus the reference zone's
	virtual double periapsis() const =0;

	///\brief The rate at which the periapsis of the orbit in this zone's
	///equatorial plane is changing.
	///
	///set_orbit() mest already have been called, and inclination() and
	///spin_frequency() must be current.
	double periapsis_evolution(
			///The torque on the orbit due to all other zones in this zone's
			///coordinate system.
			const Eigen::Vector3d &orbit_torque,

			///All torques acting on this zone (i.e. tidale, angular
			///momentum loss due to wind for the surface zone and coupling to
			///neightboring zones due to differential rotation).
			const Eigen::Vector3d &zone_torque,

			///If not Dissipation::NO_DERIV, the derivative of the rate with
			///respect to the given quantity is returned. For zone-specific
			///quantities, derivative with respect to this zone's quantity is
			///computed. If derivatives with respect to other zone's
			///quantities are required, those only come in through the orbit
			///torque and external torque, so pass the corresponding
			///derivative instead of the actual torques, and ignore this and
			///subsequent arguments.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			
			///This argument is required if dervi is neither NO_DERIV nor
			///PERIAPSIS, and should contain the derivative of the orbital
			///torque relative to the quantity identified by deriv.
			const Eigen::Vector3d &orbit_torque_deriv=
				Eigen::Vector3d(),

			///This argument is required if deriv is neither NO_DERIV nor
			///PERIAPSIS, and shoul contain the derivative of the zone
			///torque relative to the quantity identified by deriv.
			const Eigen::Vector3d &zone_torque_deriv=
				Eigen::Vector3d());

	///\brief The rate at which the inclination between this zone and the 
	///orbit is changing.
	///
	///set_orbit() mest already have been called, and inclination() and
	///spin_frequency() must be current.
	double inclination_evolution(
			///The torque on the orbit due to all other zones in this zone's
			///coordinate system.
			const Eigen::Vector3d &orbit_torque,

			///All torques acting on this zone (i.e. tidal, angular
			///momentum loss due to wind for the surface zone and coupling to
			///neightboring zones due to differential rotation).
			const Eigen::Vector3d &zone_torque,
			
			///If not Dissipation::NO_DERIV, the derivative of the rate with
			///respect to the given quantity is returned. For zone-specific
			///quantities, derivative with respect to this zone's quantity is
			///computed. If derivatives with respect to other zone's
			///quantities are required, those only come in through the orbit
			///torque and external torque, so pass the corresponding
			///derivative instead of the actual torques, and ignore this and
			///subsequent arguments.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			
			///This argument is required if dervi is neither NO_DERIV nor
			///PERIAPSIS, and shoul contain the derivative of the orbital
			///torque relative to the quantity identified by deriv.
			const Eigen::Vector3d &orbit_torque_deriv=
				Eigen::Vector3d(),

			///This argument is required if dervi is neither NO_DERIV nor
			///PERIAPSIS, and shoul contain the derivative of the zone
			///torque relative to the quantity identified by deriv.
			const Eigen::Vector3d &zone_torque_deriv=
				Eigen::Vector3d());


	///\brief Should return true iff the given term is presently locked.
	virtual bool locked(int orbital_frequency_multiplier,
			int spin_frequency_multiplier) const =0;

	///\brief Should return the tidal phase lag time the love number for the
	///given tidal term (or one of its derivatives).
	///
	///In case the specified term is in a lock, it should return the phase
	///lag for the case of the spin frequency approaching the lock from
	///below. The lag for spin frequency approaching from above should be
	///written to above_lock_value. If the term is not locked 
	///leave above_lock_value untouched.
	virtual double modified_phase_lag(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,
			
			///The current orbital frequency in rad/day
			double orbital_frequency,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::Derivative deriv,

			///If the lag of a locked term is calculated this should be set
			///to the lag assuming the spin frequency is just above the lock.
			///Otherwise, leave untouched.
			double &above_lock_value) const =0;

	///\brief Should return the corresponding component of the love
	///coefficient (Lai 2012 Equation 24).
	virtual double love_coefficient(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::Derivative deriv) const =0;

	///The moment of inertia of the zone or its age derivative.
	virtual double moment_of_inertia(
			///What to return:
			/// - 0 The moment of inertia in \f$M_\odot R_\odot^2\f$
			/// - 1 The rate of change of the moment of inertia in 
			///     \f$M_\odot R_\odot^2/Gyr\f$
			/// - 2 The second derivative in \f$M_\odot R_\odot^2/Gyr^2\f$
			int deriv_order=0);

	///\brief The spin frequency of the given zone.
	double spin_frequency() const {return __spin_frequency;}

	///The angular momentum of the given zone in \f$M_\odot R_\odot^2\f$.
	double angular_momentum() const 
	{return __spin_frequency*moment_of_inertia();}

	///\brief Set the spin frequency of the given zone.
	void set_spin_frequency(double spin_frequency)
	{__spin_frequency=spin_frequency;}

	///\brief The dimensionless tidal power or one of its derivatives.
	double tidal_power(
			///What to return
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,

			///If a spin-orbit lock is in effect and the time-lag is
			///discontinuous near zero forcing frequency, two possible values
			///can be calculated, assuming that the zone spin frequency
			///approaches the lock from below (default) or from above.
			bool above=false) const
	{return __power[2*deriv+(above? 1 : 0)];}

	///\brief The dimensionless tidal torque along x.
	///
	///See tidal_power() for a description of the arguments.
	double tidal_torque_x(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			bool above=false) const
	{return __torque_x[2*deriv+(above? 1 : 0)];}

	///\brief The dimensionless torque along y.
	///
	///See tidal_power() for a description of the arguments.
	double tidal_torque_y(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			bool above=false) const
	{return __torque_y[2*deriv+(above? 1 : 0)];}

	///\brief The dimensionless tidal torque along z.
	///
	///See tidal_power() for a description of the arguments.
	double tidal_torque_z(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			bool above=false) const
	{return __torque_z[2*deriv+(above? 1 : 0)];}

	///\brief Reads the eccentricity expansion coefficients of \f$p_{m,s}\f$.
	///
	///The given file should have been generated by
	///tabulate_eccentricity_expansion_coefficients.py.
	static void read_eccentricity_expansion(const std::string &fname)
	{__pms.read(fname);}

	///\brief Current outer radius of the zone or its derivative.
	///
	///The outermost zone's outer radius is considered to be the radius of
	///the body.
	virtual double outer_radius(
			///What to return:
			/// - 0 The boundary in \f$R_\odot\f$
			/// - 1 The rate of change of the boundary in \f$R_\odot/Gyr\f$
			/// - 2 The second derivative in \f$R_\odot/Gyr^2\f$
			int deriv_order=0) const =0;

	///\brief Mass coordinate of the zone's outer ouboundary or its
	///derivative.
	///
	///The outermost zone's boundary is considered to be the mass of the
	///body and should be constant.
	virtual double outer_mass(
			///What to return:
			/// - 0 The boundary in \f$M_\odot\f$
			/// - 1 The rate of change of the boundary in \f$M_\odot/Gyr\f$
			/// - 2 The second derivative in \f$M_\odot/Gyr^2\f$
			int deriv_order=0) const =0;

};

///Transforms a vector betwen the coordinates systems of two zones.
Eigen::Vector3D zone_to_zone_transform(
		///The zone whose coordinate system the vectors are currently in.
		const DissipatingZone &from_zone,

		///The zone whose coordinate system we want to transform the vectors
		///to.
		const DissipationgZone &to_zone, 
		
		///The vector to transform.
		const Eigen::Vector3d &vector

		///Derivatives with respect to inclination and periapsis can be
		///computed (assuming vector does not depend on these), in addition
		///to the regular transform. It is an error to request another
		///derivative.
		Dissipation::Derivative deriv=Dissipation::NO_DERIV,

		///If deriv is not NO_DERIV, derivatives can be computed with respect
		///to quantities of the from_zone (if this argument is true) or the
		///to_zone (if false).
		bool with_respect_to_from=false);

#endif
