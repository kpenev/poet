/**\file
 *
 * \brief Defines the Planet class.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __PLANET_H
#define __PLANET_H

#include "Functions.h"
#include "StellarDissipation.h"
#include "AstronomicalConstants.h"
#include <string>

///\brief Describes a planet around a star.
///
///\ingroup StellarSystem_group
class Planet : public DissipatingBody {
private:
	///The star this planet orbits.
	const Star &__star;

	///The semimajor axis (in AU) as a function of age (in Gyrs).
	InterpolatingFunctionALGLIB *__semimajor;

	
	double __mass,///< The mass of the planet in Solar masses.
		   __radius, ///< The radius of the planet in Solar radii.

		   ///\brief The age at which the planet will die.
		   ///
		   ///The end of the lifetime can be either due to the orbit
		   ///shrinking beyond the roche radius or due to the star leaving
		   ///the main sequence.
		   __lifetime,

		   ///The semimajor axis of the planet at the present time.
		   __current_semimajor,

		   __mstar, ///< stellar mass in physical units
	       __rroche, ///< tidal destruction radius in physical units
	       __mplanet, ///< planet mass in physical units
	       __rplanet; ///< planet radius in physical units

public:
	///Create a planet.
	Planet(
			///The star to place the planet around.
			const Star &star,

			///The mass of the new planet in Jupiter masses.
			double mass,
			
			///The radius of the new planet in Jupiter radii.
			double radius,

			///The planet-star orbital separation in AU.
			double semimajor);
	
	///Returns the mass of the planet in Solar masses.
	double mass() const {return __mass;}

	///Returns the radius of the planet in Solar radii.
	double radius() const {return __radius;}

	///Needed for DissipatingBody.
	double moment_of_inertia() const {return 0;}

	///\brief Returns semimajor axis (in AU) of the planet's orbit at the
	///given age.
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double semimajor(double age) const {return (*__semimajor)(age);}

	///\brief Returns the age derivative of the orbital semimajor axis.
	///
	///Units: AU/Gyr
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double semimajor_derivative(double age) const;

	///\brief Returns the orbital period of the planet at the given age.
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double orbital_period_age(double age) const 
	{return orbital_period_semimajor((*__semimajor)(age));}

	///\brief The orbital period that corresponds to the given semimajor
	///axis.
	///
	///Units: days
	double orbital_period_semimajor(double semimajor) const
	{return 2.0*M_PI/orbital_angular_velocity_semimajor(semimajor);}

	///\brief The angular velocity of the planet when the system has the
	///given age.
	///
	///Units: rad/day.
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double orbital_angular_velocity_age(double age) const
	{return orbital_angular_velocity_semimajor((*__semimajor)(age));}

	///\brief The orbital angular velocity that corresponds to the given
	///semimajor axis. 
	///
	///Units: rad/day
	double orbital_angular_velocity_semimajor(double semimajor) const;

	///\brief The derivative of the orbital angular velocity with 
	///respect to the semimajor axis.
	///
	///Units: \f$\frac{\mathrm{rad}}{\mathrm{day}\cdot\mathrm{AU}}\f$.
	double orbital_angular_velocity_semimajor_deriv(double semimajor) 
		const;

	///\brief The age of the planetary system (in Gyr) at which the planet is
	///no more.
	///
	///See the description of the lifetime member.
	double lifetime() const {return __lifetime;}

	///\brief Sets the evolution for the semimajor axis of the planet, and
	///its derivative.
	///
	///The semimajor axis values should be negative after the planet has
	///inspiralled into the star.
	///
	///\todo The implementation of this function is inconsistent with the new
	///scheme. In particular the requirement for negative semimajor axis
	///after death should be replaced by requiring NaN and the lifetime in
	///this case should be set to the duplicate age at which NaN first
	///appears.
	void set_semimajor_evolution(const std::valarray<double> &ages,
			const std::valarray<double> &semimajor_values,
			const std::valarray<double> &semimajor_derivatives);

	///\brief Returns the semimajor axis below which the planet is considered
	///destroyed.
	///
	///It is the larger of the surface of the star and the roche radius for
	///this particular planet.
	///
	///Units: AU.
	double minimum_semimajor(double age) const;

	///\brief The rate at which the semimajor axis is shrinking at the given
	///age. 
	///
	///Units: AU/Gyr
	///
	///The semimajor axis evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double tidal_decay(double age) const
	{return semimajor_derivative(age);}

	///\brief The orbital angular momentum corresponding to the given
	///semimajor axis.
	///
	///Units: \f$ M_\odot \cdot R_\odot^2\cdot\mathrm{rad}/\mathrm{day}\f$
	double orbital_angular_momentum(double semimajor) const;

	///\brief The rate of change of the orbital angular momentum.
	///
	///Units: \f$ \frac{M_\odot \cdot R_\odot^2\cdot\mathrm{rad}}
	/// {\mathrm{day}\cdot\mathrm{Gyr}}\f$
	double orbital_angmom_deriv(
			///Semimajor axis in AU.
			double semimajor, 

			///Age derivative of the samimajor axis in AU/Gyr.
			double semimajor_deriv) const;

	///\brief The rate of change of the orbital angular momentum.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2\cdot\mathrm{rad}}
	/// {\mathrm{day}\cdot\mathrm{Gyr}}\f$
	///
	///The semimajor axis evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double orbital_angmom_deriv(double age) const;

	///The semimajor axis of the planetary orbit at the present time in AU.
	double current_semimajor() const {return __current_semimajor;}

	///\brief A function defining the dissipation efficiency of the body.
	///
	///For planets this is identically zero.
	virtual double modified_phase_lag(int, double, const SpinOrbitLockInfo &,
			PhaseLag::Derivative =PhaseLag::NO_DERIV) const
	{return 0;}

	///\brief The spin frequency of the planet corresponding to the given
	///angular momentum. For now only a place holder.
	double spin_frequency(
			///The age of the planet in Gyr.
			double,

			///The angular momentum in \f$M_\odot R_\odot^2 rad/day\f$.
			double) const
	{throw Error::NotImplemented("Planet::spin_frequency()");}

	///\brief The partial angular momentum derivative of the spin frequency
	///of the planet for a given age.
	///
	///Units: \f$\left(M_\odot \cdot R_\odot^2\right)^{-1}\f$
	///
	///For now just a placeholder.
	double spin_frequency_angmom_deriv(
			///The age of the planet in Gyr.
			double,

			///The angular momentum in \f$M_\odot R_\odot^2 rad/day\f$.
			double) const
	{throw Error::NotImplemented("Planet::spin_frequency_angmom_deriv()");}
};


#endif
