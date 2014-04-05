#ifndef __STAR_H
#define __STAR_H

#include "Functions.h"
#include "StellarZone.h"
#include "StellarEvolution.h"

#include <valarray>
#include <string>

/**\file
 *
 * \brief Defines the StarBase class.
 * 
 * \ingroup StellarSystem_group
 */

///Tags identifying the state of the stellar wind.
enum WindSaturationState {
	NOT_SATURATED=-1,	///< The wind is not saturated.
	UNKNOWN,			///< The state of the wind is not known.
	SATURATED			///< The wind is saturated.
};

///More civilized output for WindSaturationState variables.
std::ostream &operator<<(std::ostream &os,
		const WindSaturationState &wind_state);

///\brief Describes a star hosting a planet apart from \f$Q^*\f$.
///
///\ingroup StellarSystem_group
class StarBase {
private:
	///The mass of the star in \f$M_\odot\f$.
	double mass,

		   ///The present age of the star in Gyr.
		   age,

		   ///\brief The frequency of the stellar convective zone while the
		   ///disk is present.
		   ///
		   ///Units: rad/day
		   disk_lock_frequency,
		   
		   ///The age (in Gyrs) when the circumstellar disk dissipates.
		   disk_dissipation_age;

	///\brief Is this a low mass star according to the stellar evolution with
	///which it was constructed?
	bool low_mass;

	const EvolvingStellarQuantity
	///The radius of the star in \f$R_\odot\f$ as a function of age in Gyr.
	*radius,

	///\brief The luminosity of the star in \f$L_\odot\f$ as a function of
	///age in Gyr.
	///
	///Since the luminosity does not enter in the equations solved it might
	///not be defined. If that is the case, this member should be NULL.
	*luminosity,

	///\brief The moment of inertia of the stellar convective envelope in
	///units of \f$M_\odot \cdot R_\odot^2\f$ as a function of age in Gyrs.
	*conv_moment_of_inertia,

	///\brief The moment of inertia of the radiative core in units of
	/// \f$M_\odot \cdot R_\odot^2\f$ as a function of age in Gyrs.
	*rad_moment_of_inertia,

	///The mass in the radiative core as a function of age in \f$M_\odot\f$
	*rad_mass,

	///\brief The radius of the radiative core in \f$R_\odot\f$ as a function
	///of age in Gyr.
	*rad_radius;

	InterpolatingFunctionALGLIB
	///\brief The angular momentum of the convective envelope of the star in
	///units units of
	/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$ as a
	///function of age in Gyrs.
	*conv_angular_momentum,

	///\brief The angular momentum of the radiative core of the star in units
	///of \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$ as a
	///function of age in Gyrs.
	*rad_angular_momentum;

	///The age at which the star leaves the main sequence in Gyrs.
	double lifetime,

		   ///\brief The normalization constant (K) in the magnetic wind
		   ///equation in
		   /// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
		   /// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
		   magnetic_wind_strength, 

		   ///\brief The saturation frequency (rad/day) in the magnetic wind
		   ///equation.
		   magnetic_wind_saturation_freq,

		   ///\brief The timescale on which the rotation of the core and the
		   ///envelope are coupled in Gyr.
		   core_env_coupling_timescale,

		   ///\brief The age at which the core first forms in Gyr.
		   core_formation,

		   ///\brief  The current convective zone angular momentum in
		   /// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
		   current_conv_angular_momentum,

		   ///\brief The current radiative core angular momentum in 
		   /// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
		   current_rad_angular_momentum;
public:
	///Create a star with the given properties.
	StarBase(
			///Mass of the star
			double current_mass,

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
			double current_rad_spin=NaN);
	
	///The current age of the star in Gyrs.
	double current_age() const;

	///The mass of the star in solar masses.
	double get_mass() const;

	///True iff the star is considered a low mass star by stellar evolution.
	bool is_low_mass() const {return low_mass;}

	///The radius of the star (in \f$R_\odot\f$) for the given age (in Gyr).
	double get_radius(double age) const;

	///The luminosity of the star in \f$L_\odot\f$ at the given age in (Gyr).
	double get_luminosity(double age) const;

	///\brief The age derivative of the log(radius) of the star (in
	/// \f$\mathrm{Gyr}^{-1}\f$) for the given age (in Gyr).
	double get_logradius_deriv(double age) const;

	//The present radius of the star (in \f$R_\odot\f$).
//	double get_radius() const;

	///The age at which the radiative core first forms in Gyr.
	double core_formation_age() const;

	///\brief The frequency to which the stellar convective zone is locked by
	///the disk (in rad/day).
	double get_disk_lock_frequency() const;

	///\brief Sets the frequency to which the stellar convective zone is
	///locked by the disk (in rad/day)
	void set_disk_lock_frequency(double value);

	///The age at which the circumstellar disk dissipates.
	double get_disk_dissipation_age() const;

	///The radius of the radiative zone.
	double get_rad_radius(double age) const;

	///The derivative of the radiative zone radius.
	double get_rad_radius_deriv(double age, unsigned order=1) const;

	///The mass of the radiative core.
	double get_rad_mass(double age) const;

	///The derivative of the radiative core mass.
	double get_rad_mass_deriv(double age) const;

	///\brief The moment of inertia of a stellar zone at the given age
	///(in Gyr).
	///
	///Units: \f$M_\odot \cdot R_\odot^2\f$
	double moment_of_inertia(double age, StellarZone zone) const;

	///\brief The age derivative of the moment of inertia of a stellar zone
	///at the given age (in Gyr).
	///
	///Units: \f$M_\odot \cdot R_\odot^2/\mathrm{Gyr}^\mathrm{order}\f$
	double moment_of_inertia_deriv(
			///Stellar age in Gyr.
			double age,
			
			///The zone for which we need the moment of inertia derivative.
			StellarZone zone,

			///The order of the derivative desired (if >2 result is zero).
			int order=1) const;

	///\brief The angular momentum of a stellar zone at the given age 
	///(in Gyrs).
	///
	///Units: \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
	double get_angular_momentum(double age, StellarZone zone) const;

	///\brief The derivative of the angular momentum of the specified zone.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double angular_momentum_deriv(double age, StellarZone zone, int order);

	///\brief The spin frequency of the specified zone of the star in
	///rad/day for the given age in Gyr.
	double spin_frequency(double age, StellarZone zone) const;

	///\brief The spin frequency of a stellar zone that corresponds to a
	///given angular momentum.
	///
	///Units: rad/day
	double spin_frequency(
			///The age at which the spin frequency is required in Gyr.
			double age,
			
			///The zone whose frequency is required.
			StellarZone zone,
			
			///The angular momentum of the zone in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum) const;

	///\brief The partial age derivative of the spin frequency of a stellar
	///zone for a specified angular momentum.
	///
	///Units: \f$\frac{\mathrm{rad}}{\mathrm{day}\cdot\mathrm{Gyr}}\f$
	double spin_frequency_age_deriv(
			///The age at which the spin frequency derivative is required in
			///Gyr.
			double age,
			
			///The zone whose frequency derivative is required.
			StellarZone zone,

			///The angular momentum of the zone in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum) const;


	///\brief The partial angular momentum derivative of the spin frequency
	///of a stellar zone for a given age.
	///
	///Units: \f$\left(M_\odot \cdot R_\odot^2\right)^{-1}\f$
	double spin_frequency_angmom_deriv(
			///The age at which the spin frequency derivative is required in
			///Gyr.
			double age,
			
			///The zone whose frequency derivative is required.
			StellarZone zone,

			///The angular momentum of the zone in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum) const;

	///\brief The spin period of the specified zone of the star at the given
	///age.
	///
	///Units: days
	double spin_period(
			///The age at which the spin period is required in Gyr.		
			double age,

			///The zone whose spin period is required.
			StellarZone zone) const;

	///\brief Prescribes the angular momentum evolution of the specified
	///stellar zone.
	void set_angular_momentum_evolution(
			///The ages at which the angular momentum is known in Gyr.
			const std::valarray<double> &ages, 

			///The angular momentum values at each of the ages in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			const std::valarray<double> &L_values,

			///The full age derivatives of the angular momentum at each of
			///the ages in \f$\frac{M_\odot \cdot R_\odot^2 \cdot
			/// \mathrm{rad}}{\mathrm{day}\cdot \mathrm{Gyr}}\f$
			const std::valarray<double> &L_derivatives,

			///The stellar zone whose angular momentum evolution is being
			///defined.
			StellarZone zone);

	///\brief Returns the mass of the specified zone (in \f$M_\odot\f$) at
	///the specified age in Gyr.
	double get_zone_mass(double age, StellarZone zone) const;

	///Retruns the timescale for core-envelope coupling in Gyrs.
	double get_core_env_coupling_timescale() const
	{return core_env_coupling_timescale;}

	///\brief Retuns the strength of the magnetic wind (the K constant).
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	double get_wind_strength() const {return magnetic_wind_strength;}

	///Retruns the saturation frequency of the magnetic wind in rad/day.
	double get_wind_saturation_frequency() const
	{return magnetic_wind_saturation_freq;}

	///\brief Returns the torque on the stellar envelope due to the stellar 
	///magnetic wind.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double wind_torque(
			///The age of the star in Gyr.
			double age,
			
			///The convective zone spin frequency in rad/day.
			double conv_frequency,

			///The saturation state of the wind to assume. Could be UNKNOWN
			///(default), in which case it is determined by comparing the
			///convective spin frequency to the saturation frequency.
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///\brief Partial derivative of the wind torque with respect to the
	///convective zone spin frequency.
	///
	///Units: \f$M_\odot \cdot R_\odot^2/\mathrm{Gyr}\f$
	double wind_torque_freq_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Spin frequency of the convective zone in rad/day.
			double conv_frequency,

			///The saturation state of the wind to assume. Could be UNKNOWN
			///(default), in which case it is determined by comparing the
			///convective spin frequency to the saturation frequency.
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

	///\brief Partial derivative of the wind torque with respect to the
	///stellar age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	double wind_torque_age_deriv(
			///Stellar age in Gyr.
			double age,

			///Convective zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum,

			///If true the partial derivative is taken with the convective
			///zone angular momentum held constant, otherwise the convective
			///zone spin frequency is held constant.
			bool const_angular_momentum=true,

			///The saturation state of the wind to assume. Could be UNKNOWN
			///(default), in which case it is determined by comparing the
			///convective spin frequency to the saturation frequency.
			WindSaturationState assume_wind_saturation=UNKNOWN) const;

	///\brief The wind torque at the given age in Gyrs.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///The angular momentum evolution for the convective zone should already
	///be specified by calling the set_angular_momentum_evolutio method with
	///zone=convective
	double wind_torque(double age) const;

	///\brief The torque on the stellar envelope due to the core-envelope
	///coupling.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double differential_rotation_torque_angmom(
			///Stellar age in Gyr.
			double age,
			
			///Angular momentum of the convective envelope in units of:
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv, 

			///Angular momentum of the radiative core in units of:
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lrad) const;

	///Returns the derivative (with respect to the angular momentum of 
	///the zone specified by the with_respect to parameter or time if 
	///that parameter is omitted) of the torque on the stellar envelope 
	///due to the core-envelope coupling if the angular momenta of the 
	///core and the envelope are as specified. The units of the torque 
	///are Msun*Rsun^2*radians/(day*Gyr)

	///\brief The partial derivative of the differential rotation torque.
	///
	///Units: \f$Gyr^{-1}\f$ or 
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	///
	///Depending on the with_respect_to argument, the derivative is with
	///respect to one of the angular momenta or age, with the other two
	///variables held constant.
	double differential_rotation_torque_deriv(
			///The stellar age in Gyr.
			double age,
			
			///The convective zone angular momentum.
			double Lconv, 

			///The radiative zone angular momentum.
			double Lrad,
			
			///Identifies either the zone whose angular momentum we are
			///taking the pratial derivaite with respect to, if it is
			//(convective or radiative) or that it should be taken with
			///respect to the stellar age if with_respect_to is total.
			StellarZone with_respect_to=total) const;

	///\brief The differential rotation torque on the stellar envelope for
	///the given amount of differential rotation.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double differential_rotation_torque(
			///The stellar age in Gyr.
			double age, 

			///The amount of differential rotation in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double differential_rotation_amount, 

			///The spin frequency of the convective zone in rad/day.
			double conv_frequency) const;

	///\brief The partial derivative of the differential rotation torque on
	///the stellar envelope with respect to whatever the derivatives of the
	///differential rotation and convective frequency are.
	///
	///The units of the torque are \f$\frac{M_\odot \cdot R_\odot^2 \cdot
	/// \mathrm{rad}}{\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double differential_rotation_torque_deriv(
			///The stellar ge in Gyr.
			double age, 

			///The amount of differential rotation in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double differential_rotation_amount, 

			///The derivative of the differential rotation with respect to
			///whatever we want the torque differentiated against.
			double differential_rotation_deriv,

			///The spin frequency of the convective envelope in rad/day.
			double conv_frequency,

			///The derivative of the spin frequency of the convective
			///envelope with respect to whatever we want the torque 
			///differentiated against.
			double conv_frequency_deriv,
			
			///Set to true if the derivatives of the differential rotation
			///and the convective zone frequency are with respect to age.
			bool with_respect_to_age=false)
		const;

	///\brief The torque on the stellar envelope due to the core-envelope
	///coupling.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///Both the core and envelope angular momenta evolutions must be already
	///specified.
	double differential_rotation_torque(double age) const;

	///\brief The amount of differential rotation between the convective
	///envelope and the radiative core.
	///
	///Units: \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
	//
	//at the given age that corresponds
	///to the specified angular momenta of the two zones.
	double differential_rotation(
			///Stellar age in Gyr.
			double age,
			
			///Convective envelope angular momentum in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv,

			///Radiative core angular momentum in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$		
			double Lrad) 
		const;

	///\brief The partial derivative of the differential rotation between
	///the convective envelope and the radiative core.
	///
	///Units: dimensionless or \f$\frac{M_\odot \cdot R_\odot^2 \cdot
	/// \mathrm{rad}}{\mathrm{day}\cdot \mathrm{Gyr}}\f$
	///
	///Depending on the with_respect_to argument, the derivative is with
	///respect to one of the angular momenta or age, with the other two
	///variables held constant.
	double differential_rotation_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Convective zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv, 

			///Radiative zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lrad,
			
			///Identifies either the zone whose angular momentum we are
			///taking the pratial derivaite with respect to, if it is
			//(convective or radiative) or that it should be taken with
			///respect to the stellar age if with_respect_to is total.
			StellarZone with_respect_to=total) const;

	///\brief The amount of differential rotation between the convective
	///envelope and the radiative core.
	///
	///The angular momenta evolution of the two stellar zones must already be
	///specified.
	double differential_rotation(
			///Stellar age in Gyr.
			double age) const;

	///\brief The rate of moment of inertia transfer from the envolope to the
	///core due to the convective-radiative boundary evolving with time.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double core_inertia_gain(double age) const;

	///\brief The age derivative of the rate of moment of inertia transfer
	///from the envolope to the core due to the convective-radiative boundary
	///evolving with time.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	double core_inertia_gain_deriv(double age) const;

	///\brief The tidal quality factor of the star for tides having the
	///specified frequency (in rad/day).
	virtual double get_tidal_Q(double tidal_frequency) const =0;

	///\brief The frequency derivative of tidal quality factor of the star
	///for tides having the specified frequency (in rad/day).
	virtual double get_tidal_Q_deriv(double tidal_frequency) const =0;

	///The age at which the star leaves the main sequence in Gyr.
	double get_lifetime() const;

	///The angular momentum of a stellar zone at the present time.
	double get_current_angular_momentum(StellarZone zone) const;

	///Delete the allocated interpolated stellar properties.
	~StarBase();
};

#endif
