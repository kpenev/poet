/**\file
 *
 * \brief Defines the Planet class.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __PLANET_H
#define __PLANET_H

#include "Functions.h"
#include "StellarQ.h"
#include "AstronomicalConstants.h"
#include <string>

///\brief Describes a planet around a star.
///
///\ingroup StellarSystem_group
class Planet {
private:
	///The star this planet orbits.
	Star* star;

	///The semimajor axis (in AU) as a function of age (in Gyrs).
	InterpolatingFunctionALGLIB *semimajor_axis;

	
	double mass,///< The mass of the planet in Jupiter masses.
		   radius, ///< The radius of the planet in Jupiter radii.

		   ///\brief The age at which the planet will die.
		   ///
		   ///The end of the lifetime can be either due to the orbit
		   ///shrinking beyond the roche radius or due to the star leaving
		   ///the main sequence.
		   lifetime,

		   ///The semimajor axis of the planet at the present time.
		   current_semimajor,

		   mstar, ///< stellar mass in physical units
	       rroche, ///< tidal destruction radius in physical units
	       mplanet, ///< planet mass in physical units
	       rplanet; ///< planet radius in physical units

public:
	///Create a planet.
	Planet(
			///The star to place the planet around.
			Star* star,

			///The mass of the new planet in Jupiter masses.
			double observed_mass,
			
			///The radius of the new planet in Jupiter radii.
			double observed_radius,

			///The planet-star orbital separation in AU.
			double observed_semimajor);

	///\brief Should not exist!!!
	///
	///essentially eliminate the planet from existence; set its mass to 1
	///Earth mass, its orbital distance to 1 AU
	void transform_into_earth();

	///Returns the mass of the planet in Jupiter masses.
	double get_mass() const {return mass;}

	///Returns the radius of the planet in Jupiter radii.
	double get_radius() const {return radius;}

	///\brief Returns semimajor axis (in AU) of the planet's orbit at the
	///given age.
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double get_semimajor(double age) const {return (*semimajor_axis)(age);}

	///\brief Returns the age derivative of the orbital semimajor axis.
	///
	///Units: AU/Gyr
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double get_semimajor_derivative(double age) const;

	///\brief Returns the orbital period of the planet at the given age.
	///
	///The orbital evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double orbital_period_age(double age) const 
	{return orbital_period_semimajor((*semimajor_axis)(age));}

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
	{return orbital_angular_velocity_semimajor((*semimajor_axis)(age));}

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
	double get_lifetime() const {return lifetime;}

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
	{return get_semimajor_derivative(age);}

	///\brief The rate at which the semimajor axis would be shrinking  due to
	///tidal dissipation.
	///
	///Units: AU/Gyr or \f$\mathrm{AU}^{6.5}/\mathrm{Gyr}\f$.
	///
	///If use_a6p5, returns the time derivative of \f$a^{6.5}\f$ in units of
	///AU instead (a better quantity to solve ODE for).
	double tidal_decay(
			///The age (in Gyr) of the star when the tidal decay should be
			///calculated.
			double age,
			
			///The semimajor axis (in AU) for which the tidal decay should be
			///calculated.
			double semimajor, 

			///The spin frequency (in rad/day) of the star for which tidal
			///decay should be calculated.
			double stellar_spin_frequency,
			
			///If true \f$ \frac{d}{dt} a^{6.5}\f$ is returned instead of
			/// \f$ \frac{d}{dt} (a) \f$.
			bool use_a6p5=false) const;

	///\brief The partial derivative of the orbital decay equation with
	///respect to the semimajor axis (^6.5).
	///
	///Units: \f$ \mathrm{Gyr}^{-1} \f$ regardless of use_a6p5.
	double tidal_decay_semimajor_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Semimajor axis in AU.
			double semimajor, 
			
			///Stellar spin frequency in rad/day
			double stellar_spin_frequency,
			
			///If true \f$\frac{\partial}{\partial a^{6.5}}
			/// \left(\frac{d}{dt} a^{6.5}\right) \f$ is returned instead of
			/// \f$ \frac{\partial}{\partial a} \left(\frac{d}{dt} a\right)
			/// \f$
			bool use_a6p5=false) const;

	///\brief The partial derivative of the orbital decay equation with
	///respect to the stellar spin frequency.
	///
	///Units: \f$\frac{\mathrm{AU}\cdot \mathrm{day}}
	/// {\mathrm{Gyr}\cdot\mathrm{rad}}\f$
	double tidal_decay_star_spin_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Semimajor axis in AU.
			double semimajor,

			///Stellar spin frequency in rad/day.
			double stellar_spin_frequency,
			
			///If true \f$ \frac{\partial}{\partial \Omega_*}
			/// \left(\frac{d}{dt} a^{6.5}\right)\f$ is returned instead of
			/// \f$ \frac{\partial}{\partial \Omega_*} \left(\frac{d}{dt} a
			/// \right)\f$
			bool use_a6p5=false) const;

	///\brief The partial derivative of the tidal torque equation with
	///respect to the stellar spin frequency.
	///
	///Units: \f$M_\odot\cdot R_\odot^2/\mathrm{Gyr}\f$
	double tidal_torque_star_spin_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Semimajor axis in AU
			double semimajor,

			///Age derivative of the samimajor axis in AU/Gyr.
			double semi_deriv,
			
			///Stellar spin frequency in rad/day.
			double stellar_spin_frequency) const;

	///\brief The partial derivative of the orbital decay equation with
	///respect to age.
	///
	///Units: \f$ \mathrm{AU}/\mathrm{Gyr}^2\f$ or
	/// \f$\mathrm{AU}^{6.5}/\mathrm{Gyr}^2\f$.
	double tidal_decay_age_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Semimajor axis in AU.
			double semimajor, 

			///Stellar spin frequency in rad/day
			double stellar_spin_frequency,

			///If true \f$ \frac{\partial}{\partial t} \left(\frac{d}{dt}
			/// a^{6.5}\right)\f$ is returned instead of
			/// \f$\frac{\partial}{\partial t} \left(\frac{d}{dt} a\right)\f$
			bool use_a6p5=false) const;

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

	///\brief The partial derivative with respect to age of the equation for
	///the rate of change of the orbital angular momentum.
	///
	///Units: \f$ \frac{M_\odot \cdot R_\odot^2\cdot\mathrm{rad}}
	/// {\mathrm{day}\cdot\mathrm{Gyr}^2}\f$
	double orbit_angmom_deriv_age_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Semimajor axis in AU.
			double semimajor,

			///Rate of change of the semimajor axis in AU/Gyr.
			double semimajor_deriv,
			
			///Stellar spin frequency in rad/day
			double stellar_spin_frequency) const;

	///\brief The partial derivative with respect to the semimajor axis of
	///the equation for the rate of change of the orbital angular momentum.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2\cdot\mathrm{rad}}
	/// {\mathrm{day}\cdot\mathrm{Gyr}\cdot\mathrm{AU}}\f$
	double orbital_angmom_deriv_semimajor_deriv(double semimajor, 
			double semimajor_deriv, double stellar_spin_frequency) const;

	///\brief The rate of change of the orbital angular momentum.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2\cdot\mathrm{rad}}
	/// {\mathrm{day}\cdot\mathrm{Gyr}}\f$
	///
	///The semimajor axis evolution must already be specified by calling the
	///set_semimajor_evolution method.
	double orbital_angmom_deriv(double age) const;

	///The semimajor axis of the planetary orbit at the present time in AU.
	double get_current_semimajor() const {return current_semimajor;}

	///\brief The full age derivative of \f$Q_*\f$ (due to \f$\Omega_*\f$
	///changing).
	///
	///Units: \f$\mathrm{Gyr}^{-1}\f$
	double get_tidal_Q_age_deriv(double age, double semi, double Lc) const;
};


#endif
