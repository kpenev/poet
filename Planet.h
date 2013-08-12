#ifndef __PLANET_H
#define __PLANET_H

#include "Functions.h"
#include "Star.h"
#include "AstronomicalConstants.h"
#include <string>

class Planet {
private:
	Star* star;

	///The semimajor axis (in AU) as a function of age (in Gyrs).
	InterpolatingFunctionALGLIB *semimajor_axis;

	
	double mass,///< The mass of the planet in Jupiter masses.
		   radius, ///< The radius of the planet in Jupiter radii.

		   ///The age at which the planet will die either due to its orbit
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
	///Create a planet with the given mass and radius and place it at the
	///given separation from its parent star.
	Planet(Star* star,
			double observed_mass, double observed_radius,
			double observed_semimajor);

	///essentially eliminate the planet from existence; set its mass to 1
	///Earth mass, its orbital distance to 1 AU
	void transform_into_earth();

	///Returns the mass of the planet in Jupiter masses.
	double get_mass() const {return mass;}

	///Returns the radius of the planet.
	double get_radius() const {return radius;}

	///Returns semimajor axis of the planet's orbit at the given age. The
	///orbital evolution must already be specified by calling the
	//set_semimajor_evolution method.
	double get_semimajor(double age) const {return (*semimajor_axis)(age);}

	///Returns the age derivative of the semimajor axis of the planet's 
	///orbit at the given age. The orbital evolution must already be 
	///specified by calling the set_semimajor_evolution method.
	double get_semimajor_derivative(double age) const;

	///Returns the orbital period of the planet at the given age. The
	///orbital evolution must already be specified by calling the
	//set_semimajor_evolution method.
	double orbital_period_age(double age) const 
	{return orbital_period_semimajor((*semimajor_axis)(age));}

	///Returns the orbital period that the planet would have, if it was
	///placed at the given semimajor axis. 
	double orbital_period_semimajor(double semimajor) const
	{return 2.0*M_PI/orbital_angular_velocity_semimajor(semimajor);}

	///Returns the angular velocity (in radians/day) of the planet when
	///the system is of the given age. The orbital evolution must 
	///already be specified by calling the
	///set_semimajor_evolution method.
	double orbital_angular_velocity_age(double age) const
	{return orbital_angular_velocity_semimajor((*semimajor_axis)(age));}

	///Returns the orbital angular velocity (in radians/day) that the 
	///planet would have, if it was placed at the given semimajor axis. 
	double orbital_angular_velocity_semimajor(double semimajor) const;

	///Returns the derivative of the orbital angular velocity with 
	///respect to the semimajor axis (in radians/day/AU) that the 
	///planet would have, if it was placed at the given semimajor axis. 
	double orbital_angular_velocity_semimajor_deriv(double semimajor) 
		const;

	///Returns the age of the planetary system at which the planet is
	///considered destroyed, see the description of the lifetime member.
	double get_lifetime() const {return lifetime;}

	///Specifies an evolution for the semimajor axis of the planet, and
	///its derivative at a set of ages. The semimajor axis values should
	///be negative after the planet has inspiralled into the star.
	void set_semimajor_evolution(const std::valarray<double> &ages,
			const std::valarray<double> &semimajor_values,
			const std::valarray<double> &semimajor_derivatives);

	///Returns the semimajor axis below which the planet is considered
	///destroyed (either the surface of the star or the roche radius,
	///whichever is larger).
	double minimum_semimajor(double age) const;

	///Returns the rate at which the semimajor axis is shrinking due to
	///the dissipation of tidal energy in the system's central star in
	///AU/Gyr. The semimajor axis evolution must already be specified by
	///calling the set_semimajor_evolution method.
	double tidal_decay(double age) const {return get_semimajor_derivative(age);}

	///Returns the rate at which the semimajor axis would be shrinking 
	///due to the dissipation of tidal energy in the system's central 
	///star in AU/Gyr if the orbit had the specified semimajor axis (in
	///AU) and the star was spinning at the given frequency 
	///(in radians/day). If use_a6p5, returns the time derivative
	///of a^6.5 in units of AU instead (a better quantity to solve ODE
	///for).
	double tidal_decay(double age, double semimajor, 
		double stellar_spin_frequency, bool use_a6p5=false) const;

	///Returns the derivative of the orbital decay with respect to the 
	///semimajor axis.
	double tidal_decay_semimajor_deriv(double age, double semimajor, 
		double stellar_spin_frequency, bool use_a6p5=false) const;

	///Returns the derivative of the orbital decay with respect to the 
	///stellar spin frequency.
	double tidal_decay_star_spin_deriv(double age, double semimajor,
					   double stellar_spin_frequency, bool use_a6p5=false) const;

	///Returns the derivative of the tidal torque with respect to the
	///stellar spin frequency.
	double tidal_torque_star_spin_deriv(double age, double semimajor,
			double semi_deriv, double stellar_spin_frequency) const;

	///Returns the age derivative of the orbital decay rate.
	double tidal_decay_age_deriv(double age, double semimajor, 
				     double stellar_spin_frequency,
				     bool use_a6p5=false) const;

	///Returns the angular momentum of the orbit if it has the given
	///semimajor axis.
	double orbital_angular_momentum(double semimajor) const;

	///Returns the rate at which the angular momentum of the orbit
	///changes if the semimajor axis has the given value and rate of
	///change. Units are Msun*Rsun^2/day/Gyr.
	double orbital_angmom_deriv(double semimajor, 
			double semimajor_deriv) const;

	///Returns the partial derivative with respect to the age
	///of the rate at which the angular momentum of the orbit changes if
	///the semimajor axis has the given value and rate of change.
	double orbit_angmom_deriv_age_deriv(double age, double semimajor,
			double semimajor_deriv, double stellar_spin_frequency) const;

	///Returns the partial derivative with respect to the current semimajor 
	///axis of the rate at which the angular momentum of the orbit changes if
	///the semimajor axis has the given value and rate of change. Units are
	///Msun*Rsun^2/(day*Gyr*AU)
	double orbital_angmom_deriv_semimajor_deriv(double semimajor, 
			double semimajor_deriv, double stellar_spin_frequency) const;

	///Returns the rate at which the angular momentum of the orbit is
	///changing for the given system age. The semimajor axis evolution 
	///must already be specified by calling the set_semimajor_evolution
	///method. Units are Msun*Rsun^2/day/Gyr.
	double orbital_angmom_deriv(double age) const;

	///Returns the semimajor axis of the planetary orbit at the present
	///time.
	double get_current_semimajor() const {return current_semimajor;}
	double get_tidal_Q_age_deriv(double age, double semi, double Lc) const;

};


#endif
