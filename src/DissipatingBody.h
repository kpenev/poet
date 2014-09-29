#ifndef __DISSIPATING_BODY_H
#define __DISSIPATING_BODY_H

#include "DissipatingZone.h"
#include "AstronomicalConstants.h"
#include <valarray>
#include <cassert>

/**\file 
 *
 * \brief Declares the DissipatingBody class.
 *
 * \ingroup StellarSystem_group
 */

///\brief A base class for any body contributing to tidal dissipation.
///
///All torques are in units of
/// \f$M_\odotR_\odot^2 \mathrm{rad}\mathrm{day}^{-1}\mathrm{Gyr}^{-1}\f$
///and tidal power is in units of 
/// \f$M_\odotR_\odot^2 \mathrm{rad}^2\mathrm{day}^{-2}\mathrm{Gyr}^{-1}\f$
///
///Age derivatives are per Gyr
///
///Masses and radii are in solar units.
///
///Angular velocities are in rad/day.
///
///\ingroup StellarSystem_group
class DissipatingBody {
private:
	std::valarray< std::valarray<Eigen::Vector3d> > 
		///\brief The rate of angular momentum transfer between two
		///neighboring zones due to zone boundaries moving.
		///
		///The outer index is the index of the outer zone and each entry
		///consists of two sub-entries giving the rate at which each of the
		///two zones (outer first, inner second) involved gains angular 
		///momentum due to the moving boundary in each zone's coordinate 
		///system.
		__angular_momentum_transfer,

		///\brief Tidal torques on each zone for spin frequency approaching
		///potential lock from above.
		///
		///See tidal_torque() for more details.
		__tidal_torques_above,

		///\brief Tidal torques on each zone for spin frequency approaching
		///potential lock from below.
		///
		///See tidal_torque() for more details.
		__tidal_torques_below;

	///\brief Total tidal power in this body (and derivatives), for spin
	///approaching a potential lock from above.
	std::valarray<double> __tidal_power_above

		///\brief Total tidal power in this body (and derivatives), for spin
		///approaching a potential lock from below.
		__tidal_power_below;

	///\brief Angular momentum transfer between two neighboring zones due to
	///moving zone boundaries.
	///
	///The argument of periapsis and the inclination must be defined for both
	///zones involved.
	void angular_momentum_transfer(
			///The outer zone.
			const DissipatingZone &outer_zone,

			///The inner zone.
			const DissipatingZone &inner_zone,

			///The angular momentum change for the outer zone in the outer
			///zone's coordinate system.
			Eigen::Vector3D &outer_angmom_gain,

			///The angular momentum change for the inner zone in the inner
			///zone's coordinate system.
			Eigen::Vector3D &inner_angmom_gain,

			///Derivatives with respect to inclination and periapsis can be
			///computed, in addition to the actual transfer. It is an error
			///to request another derivative.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,

			///If deriv is not NO_DERIV, derivatives can be computed with 
			///respect to quantities of the outer zone (if this argument is 
			///true) or the inner zone (if false).
			bool with_respect_to_outer=false);

	///\brief Rate of angular momentum transfer (or its derivatives) to a
	///zone due to its top boundary moving.
	///
	///The result is in the coordinate system of the zone.
	///
	///The set_orbit() method must already have been called and all zones
	///must have their inclinations and arguments of periapsis set.
	Eigen::Vector3d angular_momentum_transfer_from_top(
			///The index of the zone to calculate angular momentum gain for.
			///Must not be zero.
			unsigned zone_index,

			///Whether to return the quantity or one of its derivatives.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,

			///If deriv is a zone specific quantity this argument determines
			///if derivative with respect to the quantity of the zone above
			///(true) or this zone (false) should be calculated.
			bool with_respect_to_outer=false);

	///\brief Rate of angular momentum transfer (or its derivatives) to a
	///zone due to its bottom boundary moving.
	///
	///The result is in the coordinate system of the zone.
	///
	///The set_orbit() method must already have been called and all zones
	///must have their inclinations and arguments of periapsis set.
	Eigen::Vector3d angular_momentum_transfer_from_bottom(
			///The index of the zone to calculate angular momentum gain for.
			///Must not be less than number_zones()-1.
			unsigned zone_index,

			///Whether to return the quantity or one of its derivatives.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,

			///If deriv is a zone specific quantity this argument determines
			///if derivative with respect to the quantity of the zone below
			///(true) or this zone (false) should be calculated.
			bool with_respect_to_inner=false);

	///\brief Rate of angular momentum transfer (or its derivatives) to a
	///zone due to moving zone boundaries.
	///
	///The result is in the coordinate system of the zone.
	///
	///The set_orbit() method must already have been called and all zones
	///must have their inclinations and arguments of periapsis set.
	Eigen::Vector3d angular_momentum_transfer_to_zone(
			///The index of the zone to calculate angular momentum gain for.
			unsigned zone_index,

			///Whether to return the quantity or one of its derivatives.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,

			///See matching argument of external_torque() for description.
			int deriv_zone=0);

public:
	///\brief Defines the orbit this body is in.
	///
	///The inclinations and arguments of periapsis must be already set for
	///all zones.
	void set_orbit(
			///The mass of the second body in the system.
			double companion_mass,

			///The semimajor axis of the orbit in \f$R_\odot\f$.
			double semimajor,

			///The eccentricity of the orbit
			double eccentricity);

	///\brief External torque acting on a single zone (last 
	///calculate_torques_power()).
	///
	///This includes torques coupling it to other zones (due to differential
	///rotation and/or due to mass transfer) and wind loss for the
	///surface zone.
	///
	///For each zone the torques are in a coordinate system with
	/// \f$\hat{x}\f$ along the ascending node of the orbit in the zone's
	///equatorial plane, and \f$\hat{z}\f$ along the zone's angular momentum.
	const Eigen::Vector3D &external_torque(
			///The index of the zone whose torque is needed.
			unsigned zone_index,

			///Whether to return the quantity or one of its derivatives.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			
			///Since external torques depend on neighboring zones, this
			///parameter is used to distinguish those (it is ignored if :
			///deriv is Dissipation::NO_DERIV or a quantity which is not
			///zone specific):
			/// - -1 Return the derivative with respect to the quantity for
			///   the zone above.
			/// - 0 Return the derivative with respect to the quantity for
			///   this
			///   zone.
			/// - 1 Return the derivative with respect to the quantity for
			///   the zone below.
			int deriv_zone=0);

	///\brief Tidal torque acting on the given zone (last 
	///calculate_torques_power()).
	///
	///The coordinate system is the same as for external_torque().
	const Eigen::Vector3D &tidal_torque(unsigned zone_index,

			///If a zone is currently in a lock, this argument decides
			///whetherthe torque for spin frequency above (true) or below
			///(false) of the lock should be returned.
			bool above)
	{
#ifdef DEBUG
		assert(zone_index<number_zones());
#endif
		return __tidal_torques[zone_index];
	}

	///\brief Energy dissipation rate due to tides in this body (last 
	///calculate_torques_power()).
	double tidal_power() const {return __tidal_power;}

	///The number of zones the body consists of.
	virtual unsigned number_zones() const =0;

	///A modifiable reference to one of the body's zones.
	virtual DissipatingZone &zone(
			///The index of the zone within the body. Sequential zones are
			///ordered from outside to inside (0 is the surface zone,
			///number_zones()-1 is the core).
			unsigned zone_index)=0;

	///\brief Coupling torque for two neighboring zones in the coordinate
	///system of the top zone.
	///
	///Units:
	/// \f$M_\odot R_\odot^2\mathrm{rad}\mathrm{day}^{-1}mathrm{Gyr}^{-1}\f$.
	virtual Eigen::Vector3d angular_momentum_coupling(
			///The index of the outer of the two zones whose coupling is
			///needed.
			unsigned top_zone_index,
			
			///The derivative of the coupling torque required.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			
			///For derivatives with respect to zone specific quantities,
			///this determines which zone's quantity to differentiate with
			///respect to (top zone if true, bottom zone if false).
			bool with_respect_to_top=false) const =0;

	///\brief Rate of angular momentum loss by the top zone of the body and
	///its derivatives.
	///
	///The spin frequency of the topmost zone of the body must already be
	///set.
	virtual double angular_momentum_loss(
			///The derivative of the angular momentum loss required.
			Dissipation::Derivative deriv) const =0;

	///\brief The current radius or its derivative with age of the body.
	double radius(
			///The order of the derivative to return.
			int deriv_order=0)
	{return zone(0).outer_radius(deriv_order);}

	///The mass of the body (constant with age).
	double mass() {return zone(0).outer_mass(Dissipation::NO_DERIV);}
	
};

#endif
