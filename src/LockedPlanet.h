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
	bool dissipative() {return false;}
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

#ifdef DEBUG
	void configure(
			///The age to set the body to.
			double age,

			///The mass of the second body in the system.
			double companion_mass,

			///The semimajor axis of the orbit in \f$R_\odot\f$.
			double semimajor,

			///The eccentricity of the orbit
			double eccentricity,

			///The spin angular momenta of the non-locked zones of the body
			///(outermost zone to innermost).
			const double *spin_angmom,

			///The inclinations of the zones of the body (same order as 
			///spin_angmom). If NULL, all inclinations are assumed zero.
			const double *inclination=NULL,
			
			///The arguments of periapsis of the zones of the bodies (same
			///order as spin_angmom). If NULL, all periapses are assumed
			///zero.
			const double *periapsis=NULL,

			///If true, the outermost zone's spin is assumed locked to a 
			///disk and spin_angmom is assumed to start from the next zone.
			bool locked_surface=false,
			
			///If true, the outermost zone's inclination is assumed to be
			///zero and the inclination argument is assumed to start from the
			///next zone.
			bool zero_outer_inclination=false,
			
			///If true, the outermost zone's periapsis is assumed to be
			///zero and the inclination argument is assumed to start from the
			///next zone.
			bool zero_outer_periapsis=false)
	{
		DissipatingBody::configure(age, companion_mass, semimajor,
								   eccentricity, spin_angmom, inclination,
								   periapsis, locked_surface,
								   zero_outer_inclination,
								   zero_outer_periapsis);
	}
#endif

};

#endif
