#ifndef __STAR_H
#define __STAR_H

#include "Functions.h"
#include "StellarZone.h"
#include "StellarEvolution.h"

#include <valarray>
#include <string>

enum WindSaturationState {NOT_SATURATED=-1, UNKNOWN, SATURATED};

class Star {
private:
	double mass,

		   ///The present age of the star
		   age,

		   ///The frequency range (in days^-1) over which 1/Q decays to zero
		   Q_transition_width,

		   ///The frequency to which the stellar convective zone is locked
		   ///while the disk is present (in days^-1)
		   disk_lock_frequency,
		   
		   ///The age (in Gyrs) when the circumstellar disk dissipates
		   disk_dissipation_age;

	const EvolvingStellarQuantity
	///The radius of the star in solar radii as a function of age in Gyr.
	*radius,

	///The luminosity of the star in solar luminosities as a function of age
	///in Gyr. If the luminosity is not defined this member should be NULL.
	*luminosity,

	///The moment of inertia of the convective envelope of the star in
	///solar masses times solar radii squared as a function of age in 
	///Gyrs.
	*conv_moment_of_inertia,

	///The moment of inertia of the entire star in solar
	///units as a function of age in Gyrs.
	*rad_moment_of_inertia,

	///The mass in the radiative core as a function of age in solar
	///masses.
	*rad_mass,

	///The radius of the radiative core in solar radii as a function of 
	///age in Gyr.
	*rad_radius;

	InterpolatingFunctionALGLIB
	///The angular momentum of the convective envelope of the star in
	///solar units per day as a function of age in Gyrs.
	*conv_angular_momentum,

	///The angular momentum of the radiative core of the star in solar
	///units per day as a function of age in Gyrs.
	*rad_angular_momentum;

	///The tidal quality factor of the star.
	double tidal_Q,

	///The age at which the star leaves the main sequence in Gyrs.
	lifetime,

	///The normalization constant (K) in the magnetic wind equation 
	///(Eq. 3).
	magnetic_wind_strength, 
	
	///The saturation frequency (radians/day) in the magnetic wind equation
	///(Eq. 3).
	magnetic_wind_saturation_freq,
	
	///The timescale on which the rotation of the core and the envelope
	///is coupled.
	core_env_coupling_timescale,

	///The age at which the core first forms
	core_formation,
	
	///The current convective zone angular momentum
	current_conv_angular_momentum,

	///The current radiative core angular momentum
	current_rad_angular_momentum;
public:
	///Create a star with the given properties.
	///current_conv_spin, current_rad_spin in units of radians/day
	///coupling_timescale in Gyr
	Star(double current_mass, double tidal_quality, 
			double wind_strength, double wind_saturation,
			double coupling_timescale, double dissipation_transition_width,
			double disk_lock_ang_vel, double disk_lock_time,
			const StellarEvolution &evolution,
			double current_age=NaN, double current_conv_spin=NaN,
			double current_rad_spin=NaN);
	
	///Same as above, but without unnecessary arguments
/*	Star(double current_mass, double tidal_quality,
			double wind_strength, double wind_saturation,
			double coupling_timescale,
			double dissipation_transition_width,
			const StellarEvolution &evolution);*/

	double get_test_val(double age);

	///Returns the current age of the star in Gyrs.
	double current_age() const;

	///Returns the mass of the star in solar masses.
	double get_mass() const;

	///Returns the radius of the star (in solar radii) for the given age
	///(in Gyrs).
	double get_radius(double age) const;

	///Returns the luminosity of the star in solar luminosities.
	double get_luminosity(double age) const;

	///Returns the age derivative of the log(radius) of the star (in Gyr^1)
	///for the given age (in Gyrs).
	double get_logradius_deriv(double age) const;

	///Returns the present radius of the star (in solar radii).
//	double get_radius() const;

	///Returns the transition width of the tidal quality factor (specified at
	///construction).
	double get_trans_width() const;

	///The age at which the radiative core first forms
	double core_formation_age() const;

	///Returns the frequency to which the stellar convective zone is locked
	///by the disk (in day^-1)
	double get_disk_lock_frequency() const;

	///Sets the frequency to which the stellar convective zone is locked
	///by the disk (in day^-1)
	void set_disk_lock_frequency(double value);

	///Returns the age at which the circumstellar disk dissipates, releasing
	///the convective zone rotation.
	double get_disk_dissipation_age() const;

	///Returns the radius of the radiative zone
	double get_rad_radius(double age) const;

	///Returns the derivative of the radiative zone radius
	double get_rad_radius_deriv(double age, unsigned order=1) const;

	///Returns the mass of the radiative core
	double get_rad_mass(double age) const;

	///Returns the derivative of the radiative core mass
	double get_rad_mass_deriv(double age) const;

	///Returns the moment of inertia (in units of
	///Msun*Rsun^2) of a particular zone of the star or the entire star 
	///depending on the zone argument at the given age (in Gyrs).
	double moment_of_inertia(double age, StellarZone zone) const;

	///Returns the age derivative of the moment of inertia (in units of
	///Msun*Rsun^2/Gyr) of a particular zone of the star or the entire star 
	///depending on the zone argument at the given age (in Gyrs).
	double moment_of_inertia_deriv(double age, StellarZone zone,
			int order=1) const;

	///Returns the angular momentum (in units of
	///Msun*Rsun^2*radians/day) of a particular zone of the star or the 
	///entire star depending on the zone argument at the given age 
	///(in Gyrs).
	double get_angular_momentum(double age, StellarZone zone) const;

	///Returns the derivative of the angular momentum of the specified zone.
	double angular_momentum_deriv(double age, StellarZone zone, int order);

	///Returns the spin frequency of the specified zone of the star in
	///radians per day for the given age.
	double spin_frequency(double age, StellarZone zone) const;

	///Returns the spin frequency of the specified zone of the star in
	///radians per day for the given age if the angular momentum of the
	///desired zone has the prescribed value. The angular momentum should
	///be in units of Msun Rsun^2/day
	double spin_frequency(double age, StellarZone zone,
			double angular_momentum) const;

	///Returns the spin frequency age derivative of the specified zone 
	///of the star in radians per day for the given age if the angular 
	///momentum of the desired zone has the prescribed value. The angular
	///momentum should be in units of Msun Rsun^2/day
	double spin_frequency_age_deriv(double age, StellarZone zone,
			double angular_momentum) const;

	///Returns the spin frequency derivative with respect to the angular
	///momentum of the specified zone of the star in radians per day for
	///the given age if the angular momentum of the desired zone has the
	///prescribed value. The angular momentum should be in units of 
	///Msun Rsun^2/day
	double spin_frequency_angmom_deriv(double age, StellarZone zone,
			double angular_momentum) const;

	///Returns the spin period of the specified zone of the star in days
	///at the given age.
	double spin_period(double age, StellarZone zone) const;

	///Prescribes an angular momentum evolution for the specified zone of
	//the star as a set of values at specified ages.
	void set_angular_momentum_evolution(
			const std::valarray<double> &ages, 
			const std::valarray<double> &L_values,
			const std::valarray<double> &L_derivatives,
			StellarZone zone);

	///Returns the mass of the specified zone (in solar masses) at the
	///specified age.
	double get_zone_mass(double age, StellarZone zone) const;

	///Retruns the timescale for core-envelope coupling in Gyrs.
	double get_core_env_coupling_timescale() const
	{return core_env_coupling_timescale;}

	///Retuns the strength of the magnetic wind (the K constant) in units of
	///Msun/Rsun^2*day^2/Gyr.
	double get_wind_strength() const {return magnetic_wind_strength;}

	///Retruns the saturation frequency of the magnetic wind in rad/day.
	double get_wind_saturation_frequency() const
	{return magnetic_wind_saturation_freq;}

	///Returns the torque on the stellar envelope due to the stellar 
	///magnetic wind, if the convective zone is spinning at the given
	///frequency (in radians/day). The units of the torque are
	///Msun*Rsun^2*radians/(day*Gyr)
	double wind_torque(double age, double conv_frequency,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///Returns the derivative of the torque on the stellar envelope due 
	///to the stellar magnetic wind with respect to the convective spin 
	///frequency, if the convective zone is spinning at the given
	///frequency (in radians/day). The units of the torque are
	///Msun*Rsun^2*radians/(day*Gyr)
	double wind_torque_freq_deriv(double age, double conv_frequency,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	/*
	///Returns the derivative of the torque on the stellar envelope due 
	///to the stellar magnetic wind with respect to stellar age, if the 
	///convective zone is spinning at the given frequency 
	///(in radians/day). The units of the torque are 
	///Msun*Rsun^2*radians/(day*Gyr)
	double wind_torque_age_deriv(double age, double conv_frequency) 
		const;
	*/

	///Returns the derivative of the torque on the stellar envelope due
	///to the stellar magnetic wind with respect to stellar age, if the
	///convective zone has the given angular momentum (in Msun*Rsun^2*rad/day)
	///The units of the torque are Msun*Rsun^2*radians/(day*Gyr). If 
	///const_angular_momentum is false, takes the age derivative assuming a
	///constant frequency.
	double wind_torque_age_deriv(double age, double angular_momentum,
			bool const_angular_momentum=true,
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///Returns the torque on the stellar envelope due to the stellar 
	///magnetic wind at the given age (the angular momentum evolution
	///for the convective zone should already be specified by calling the
	//set_angular_momentum_evolutio method with zone=convective). The 
	///units of the torque are Msun*Rsun^2*radians/(day*Gyr)
	double wind_torque(double age) const;

	///Returns the torque on the stellar envelope due to the
	///core-envelope coupling if the angular momenta of the core and the
	///envelope are as specified. The units of the torque are 
	///Msun*Rsun^2*radians/(day*Gyr)
	double differential_rotation_torque_angmom(double age, double Lconv, 
			double Lrad) const;

	///Returns the derivative (with respect to the angular momentum of 
	///the zone specified by the with_respect to parameter or time if 
	///that parameter is omitted) of the torque on the stellar envelope 
	///due to the core-envelope coupling if the angular momenta of the 
	///core and the envelope are as specified. The units of the torque 
	///are Msun*Rsun^2*radians/(day*Gyr)
	double differential_rotation_torque_deriv(double age, double Lconv, 
			double Lrad, StellarZone with_respect_to=total) const;

	///Returns the torque on the stellar envelope due to the
	///core-envelope coupling for the given amount of differential
	///rotation. The units of the torque are Msun*Rsun^2*radians/(day*Gyr)
	double differential_rotation_torque(double age, 
			double differential_rotation_amount, 
			double conv_frequency) const;

	///Returns the derivative of torque on the stellar envelope due to 
	///the core-envelope coupling with respect to whatever the 
	///derivatives of the differential rotation and convective frequency
	///are. The units of the  torque are Msun*Rsun^2*radians/(day*Gyr)
	double differential_rotation_torque_deriv(double age, 
			double differential_rotation_amount, 
			double differential_rotation_deriv,
			double conv_frequency,
			double conv_frequency_deriv, bool with_respect_to_age=false)
		const;

	///Returns the torque on the stellar envelope due to the
	///core-envelope coupling at the given age. Both the core and
	///envelope angular momenta evolutions must be already specified. 
	///The units of the torque are Msun*Rsun^2*radians/(day*Gyr)
	double differential_rotation_torque(double age) const;

	///Returns the amount of differential rotation between the convective
	///envelope and the radiative core at the given age that corresponds
	///to the specified angular momenta of the two zones.
	double differential_rotation(double age, double Lconv, double Lrad) 
		const;

	///Returns the derivative of the amount of differential rotation 
	///between the convective envelope and the radiative core (with 
	///respect to the angular momentum of the zone specified by the 
	///with_respect_to parameter or age if that parameter is omitted) at
	///the given age that corresponds to the specified angular momenta 
	///of the two zones.
	double differential_rotation_deriv(double age, double Lconv, 
		double Lrad, StellarZone with_respect_to=total) const;

	///Returns the amount of differential rotation between the convective
	///envelope and the radiative core at the given age. The angular
	///momenta evolution of the two stellar zones must already be
	//specified.
	double differential_rotation(double age) const;

	///Return the rate at which moment of inertia is transferred from the
	///envolope to the core due to the convective-radiative boundary 
	///evolving with time
	double core_inertia_gain(double age) const;

	///Return the age derivative of the rate at which moment of inertia 
	///is transferred from the envolope to the core due to the 
	///convective-radiative boundary evolving with time
	double core_inertia_gain_deriv(double age) const;

	///Returns the tidal quality factor of the star for tides having the
	///specified frequency (in radians/day).
	double get_tidal_Q(double tidal_frequency) const;

	///Returns the frequency derivative of tidal quality factor of the star 
	///for tides having the specified frequency (in radians/day).
	double get_tidal_Q_deriv(double tidal_frequency) const;

	///Returns the age at which the star leaves the main sequence in
	///Gyrs.
	double get_lifetime() const;

	///Returns the angular momentum of the core at the time of
	///observation
	double get_current_angular_momentum(StellarZone zone) const;
};

#endif
