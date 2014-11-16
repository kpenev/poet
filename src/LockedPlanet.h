#ifndef __LOCKED_PLANET_H
#define __LOCKED_PLANET_H

#include "AstronomicalConstants.h"

/**\file
 * 
 * \brief Declares a class for planets that are always locked to the orbit.
 *
 * \ingroup StellarSystem_group
 */

///\brief The only zone of a LockedPlanet.
class LockedPlanetZone : virtual public DissipatingZone {
private:
	///The mass of the planet.
	double __mass,

		   ///The radius of the planet.
		   __radius;
public:
	LockedPlanetZone(double mass, double radius) :
		__mass(mass*AstroConst::jupiter_mass/AstroConst::solar_mass),
		__radius(radius*AstroConst::jupiter_radius/AstroConst::solar_radius)
	{}

	///See DissipatingZone::modified_phase_lag(), very large constant value.
	virtual double modified_phase_lag(
			///orbital_frequency_multiplier
			int ,
			int ,
			double ,
			Dissipation::Derivative ,
			double &) const
	{return 0;}

	///See DissipatingZone::love_coefficient(), always zero.
	double love_coefficient(int, int, Dissipation::Derivative) const
	{return 0;}

	///Tiny value (\f$10^{-6} M R^2\f$).
	double moment_of_inertia(
			///What to return:
			/// - 0 The moment of inertia in \f$M_\odot R_\odot^2\f$
			/// - 1 The rate of change of the moment of inertia in 
			///     \f$M_\odot R_\odot^2/Gyr\f$
			/// - 2 The second derivative in \f$M_\odot R_\odot^2/Gyr^2\f$
			int deriv_order=0) const 
	{return (deriv_order==0 ? 1e-6*__mass*std::pow(__radius, 2) : 0);}

	///The radius of the planet.
	double outer_radius(int deriv_order=0) const
	{return (deriv_order==0 ? __radius : 0);}

	///The mass of the planet.
	double outer_mass(int deriv_order=0) const
	{return (deriv_order==0 ? __mass : 0);}

	///No dissipation.
	bool dissipative() const {return false;}

	///\brief Calls the usual DissipatingZone::configure but with zero
	///inclination and periapsis.
	void configure(
			///The age to set the zone to.
			double age,

			///The angular velocity of the orbit in rad/day.
			double orbital_frequency,

			///The eccentricity of the orbit
			double eccentricity,
			
			///The absolute value of the angular momentum of the orbit.
			double orbital_angmom,

			///The angular momentum/frequency of the spin of the zone if the
			///zone is not locked (ignored it if is).
			double spin_angmom,
			
			///The inclination of the zone relative to the orbit.
			double ,
			
			///The argument of periapsis of the orbit in the equatorial
			///planet of the zone.
			double ,

			///Is spin_angmom angular momentum of freuqency? 
			bool spin_is_frequency)
	{DissipatingZone::configure(age, orbital_frequency, eccentricity,
								orbital_angmom, spin_angmom, 0, 0,
								spin_is_frequency);}
};

/**\brief Single zone non-evolving planets with huge dissipation, so they
 * always remain locked to the disk.
 */
class LockedPlanet : virtual public DissipatingBody {
private:
	///The only zone of the planet.
	LockedPlanetZone __zone;
public:
	///Create a planet with a constant mass and radius.
	LockedPlanet(double mass, double radius) : __zone(mass, radius) {};

	///The number of zones the body consists of.
	unsigned number_zones() const {return 1;}

	///Returns the only zone.
	const DissipatingZone &zone(unsigned) const {return __zone;}

	///Returns the only zone.
	DissipatingZone &zone(unsigned) {return __zone;}

	///Should never be called.
	Eigen::Vector3d angular_momentum_coupling(unsigned,
			Dissipation::Derivative =Dissipation::NO_DERIV,
			bool =false) const
	{assert(false);}

	///Always zero.
	double angular_momentum_loss(
			Dissipation::Derivative =Dissipation::NO_DERIV) const
	{return 0;}
};

#endif
