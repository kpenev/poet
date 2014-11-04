#ifndef __LOCKED_PLANET_H
#define __LOCKED_PLANET_H

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
		__mass(mass), __radius(radius) {}

	///See DissipatingZone::modified_phase_lag(), very large constant value.
	virtual double modified_phase_lag(
			///orbital_frequency_multiplier
			int ,
			int spin_frequency_multiplier,
			double forcing_frequency,
			Dissipation::Derivative deriv,
			double &above_lock_value) const
	{
		if(deriv!=Dissipation::NO_DERIV) return 0;
		if(forcing_frequency>0) return 1;
		else if(forcing_frequency<0) return -1;
		else  {
			above_lock_value=(spin_frequency_multiplier>0 ? -1 : 1);
			return -above_lock_value;
		}
	}

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
