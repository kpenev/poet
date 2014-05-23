#ifndef __STAR_H
#define __STAR_H

#include "DissipatingBody.h"
#include "Functions.h"
#include "StellarZone.h"
#include "StellarEvolution.h"

#include <valarray>
#include <string>
#include <complex>

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
class StarBase : public DissipatingBody {
private:
	///The mass of the star in \f$M_\odot\f$.
	double __mass,

		   ///The present age of the star in Gyr.
		   __age,

		   ///\brief The frequency of the stellar convective zone while the
		   ///disk is present.
		   ///
		   ///Units: rad/day
		   __disk_lock_frequency,
		   
		   ///The age (in Gyrs) when the circumstellar disk dissipates.
		   __disk_dissipation_age;

	///\brief Is this a low mass star according to the stellar evolution with
	///which it was constructed?
	bool __low_mass;

	const EvolvingStellarQuantity
		///\brief The radius of the star in \f$R_\odot\f$ as a function of
		///age in Gyr.
		*__radius,

		///\brief The luminosity of the star in \f$L_\odot\f$ as a function
		///of age in Gyr.
		///
		///Since the luminosity does not enter in the equations solved it
		///might not be defined. If that is the case, this member should be
		///NULL.
		*__luminosity,

		///\brief The moment of inertia of the stellar convective envelope in
		///units of \f$M_\odot \cdot R_\odot^2\f$ as a function of age in
		///Gyrs.
		*__conv_moment_of_inertia,

		///\brief The moment of inertia of the radiative core in units of
		/// \f$M_\odot \cdot R_\odot^2\f$ as a function of age in Gyrs.
		*__rad_moment_of_inertia,

		///\brief The mass in the radiative core as a function of age in
		///\f$M_\odot\f$
		*__rad_mass,

		///\brief The radius of the radiative core in \f$R_\odot\f$ as a
		///function of age in Gyr.
		*__rad_radius;

	InterpolatingFunctionALGLIB
		///\brief The angular momentum of the convective envelope of the star
		///in units units of
		/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$ as a
		///function of age in Gyrs.
		*__conv_angular_momentum,

		///\brief The angular momentum of the radiative core of the star in
		///units of
		/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$ as a
		///function of age in Gyrs.
		*__rad_angular_momentum;

	///The age at which the star leaves the main sequence in Gyrs.
	double __lifetime,

		   ///\brief The normalization constant (K) in the magnetic wind
		   ///equation in
		   /// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
		   /// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
		   __magnetic_wind_strength, 

		   ///\brief The saturation frequency (rad/day) in the magnetic wind
		   ///equation.
		   __magnetic_wind_saturation_freq,

		   ///\brief The timescale on which the rotation of the core and the
		   ///envelope are coupled in Gyr.
		   __core_env_coupling_timescale,

		   ///\brief The age at which the core first forms in Gyr.
		   __core_formation,

		   ///\brief  The current convective zone angular momentum in
		   /// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
		   __current_conv_angular_momentum,

		   ///\brief The current radiative core angular momentum in 
		   /// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
		   __current_rad_angular_momentum;
public:
	///Create a star with the given properties.
	StarBase(
			///Mass of the star
			double mass,

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
			double age=NaN,
			
			///The present surface rotation rate of the star in rad/day
			///(optional).
			double conv_spin=NaN,

			///The present core rotation rate of the star in rad/day
			///optional).
			double rad_spin=NaN);

	using DissipatingBody::angular_momentum;

	///\brief Returns a reference to this object properly set-up as a
	///DissipatingBody.
	StarBase &operator()(
			///The age to set the star to.
			double age,

			///The convective zone angular momentum in 
			/// \f$M_\odot R_\odot^2 rad/day\f$.
			double Lconv,

			///The inclination in rad.
			double theta)
	{__age=age; angular_momentum(Lconv); inclination(theta); return *this;}
	
	///The current age of the star in Gyrs.
	double age() const {return __age;}

	///The mass of the star in solar masses.
	double mass() const {return __mass;}

	///True iff the star is considered a low mass star by stellar evolution.
	bool is_low_mass() const {return __low_mass;}

	///The radius of the star (in \f$R_\odot\f$) for the given age (in Gyr).
	double radius(double age) const {return (*__radius)(age);}

	///The radius of the star (in \f$R_\odot\f$) for the current age.
	double radius() const {return radius(__age);}

	///The luminosity of the star in \f$L_\odot\f$ at the given age in (Gyr).
	double luminosity(double age) const;

	///The luminosity of the star in \f$L_\odot\f$ for the current age.
	double luminosity() const {return luminosity(__age);}

	///The age at which the radiative core first forms in Gyr.
	double core_formation_age() const {return __core_formation;}

	///\brief The frequency to which the stellar convective zone is locked by
	///the disk (in rad/day).
	double disk_lock_frequency() const {return __disk_lock_frequency;}

	///\brief Sets the frequency to which the stellar convective zone is
	///locked by the disk (in rad/day)
	void disk_lock_frequency(double value) {__disk_lock_frequency=value;}

	///The age at which the circumstellar disk dissipates.
	double disk_dissipation_age() const {return __disk_dissipation_age;}

	///The radius of the radiative zone for the given age.
	double rad_radius(double age) const;

	///The radius of the radiative zone at the current age.
	double rad_radius() const {return rad_radius(__age);}

	///The mass of the radiative core at the given age.
	double rad_mass(double age) const;

	///The mass of the radiative core at the current age.
	double rad_mass() const {return rad_mass(__age);}

	///\brief The moment of inertia of a stellar zone at the given age
	///(in Gyr).
	///
	///Units: \f$M_\odot \cdot R_\odot^2\f$
	double moment_of_inertia(double age, StellarZone zone) const;

	///\brief The moment of inertia of a stellar zone at the current age.
	///
	///Units: \f$M_\odot \cdot R_\odot^2\f$
	double moment_of_inertia(StellarZone zone) const
	{return moment_of_inertia(__age, zone);}

	///\brief The moment of inertia of the convective zone at the current
	///age.
	///
	///Units: \f$M_\odot \cdot R_\odot^2\f$
	///
	///Required by DissipatingBody.
	double moment_of_inertia() const
	{return moment_of_inertia(__age, convective);}

	///\brief The age derivative of the log(radius) of the star (in
	/// \f$\mathrm{Gyr}^{-1}\f$) for the given age (in Gyr).
	double logradius_deriv(double age) const;

	///\brief The age derivative of the log(radius) of the star (in
	/// \f$\mathrm{Gyr}^{-1}\f$) for the current age.
	double logradius_deriv() const {return logradius_deriv(__age);}

	///The derivative of the radiative zone radius at the given age.
	double rad_radius_deriv(double age, unsigned order=1) const;

	///The derivative of the radiative zone radius at the current age.
	double rad_radius_deriv(unsigned order=1) const
	{return rad_radius_deriv(__age, order);}

	///The derivative of the radiative core mass at the given age.
	double rad_mass_deriv(double age) const;

	///The derivative of the radiative core mass at the current age.
	double rad_mass_deriv() const {return rad_mass_deriv(__age);}

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

	///\brief The age derivative of the moment of inertia of a stellar zone
	///at the current age (in Gyr).
	///
	///Units: \f$M_\odot \cdot R_\odot^2/\mathrm{Gyr}^\mathrm{order}\f$
	double moment_of_inertia_deriv(
			///The zone for which we need the moment of inertia derivative.
			StellarZone zone,

			///The order of the derivative desired (if >2 result is zero).
			int order=1) const
	{return moment_of_inertia_deriv(__age, zone, order);}

	///\brief The angular momentum of a stellar zone at the given age 
	///(in Gyrs).
	///
	///Units: \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
	double angular_momentum(double age, StellarZone zone) const;

	///\brief The angular momentum of a stellar zone at the current age 
	///
	///Units: \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
	double angular_momentum(StellarZone zone) const
	{return angular_momentum(__age, zone);}

	///\brief The derivative of the angular momentum of the specified zone at
	///the given age in Gyr.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double angular_momentum_deriv(double age, StellarZone zone, int order);

	///\brief The derivative of the angular momentum of the specified zone at
	///the current age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double angular_momentum_deriv(StellarZone zone, int order)
	{return angular_momentum_deriv(__age, zone, order);}

	///\brief The spin frequency of the specified zone of the star in
	///rad/day for the given age in Gyr.
	double spin_frequency(double age, StellarZone zone) const;

	///\brief The spin frequency of the specified zone of the star in
	///rad/day for the current age.
	double spin_frequency(StellarZone zone) const
	{return spin_frequency(__age, zone);}

	///\brief The spin frequency of a stellar zone that corresponds to a
	///given angular momentum at the given age in Gyr.
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

	///\brief The spin frequency of a stellar zone that corresponds to a
	///given angular momentum at the current age.
	///
	///Units: rad/day
	double spin_frequency(
			///The zone whose frequency is required.
			StellarZone zone,
			
			///The angular momentum of the zone in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum) const
	{return spin_frequency(__age, zone, angular_momentum);}

	///\brief The partial age derivative of the spin frequency of a stellar
	///zone for a specified angular momentum at the given age.
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

	///\brief The partial age derivative of the spin frequency of a stellar
	///zone for a specified angular momentum at the current age.
	///
	///Units: \f$\frac{\mathrm{rad}}{\mathrm{day}\cdot\mathrm{Gyr}}\f$
	double spin_frequency_age_deriv(		
			///The zone whose frequency derivative is required.
			StellarZone zone,

			///The angular momentum of the zone in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum) const
	{return spin_frequency_age_deriv(__age, zone, angular_momentum);}

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
			///Expected but not used for uniformity.
			double) const
	{return 1.0/moment_of_inertia(age, zone);}

	///\brief The partial angular momentum derivative of the spin frequency
	///of a stellar zone for the current age.
	///
	///Units: \f$\left(M_\odot \cdot R_\odot^2\right)^{-1}\f$
	double spin_frequency_angmom_deriv(		
			///The zone whose frequency derivative is required.
			StellarZone zone,

			///The angular momentum of the zone in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double angular_momentum) const
	{return spin_frequency_angmom_deriv(__age, zone, angular_momentum);}

	///\brief The spin period of the specified zone of the star at the given
	///age.
	///
	///Units: days
	double spin_period(
			///The age at which the spin period is required in Gyr.		
			double age,

			///The zone whose spin period is required.
			StellarZone zone) const
	{return 2.0*M_PI/spin_frequency(age, zone);}

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
	double zone_mass(double age, StellarZone zone) const;

	///\brief Returns the mass of the specified zone (in \f$M_\odot\f$) at
	///the current age.
	double zone_mass(StellarZone zone) const {return zone_mass(__age, zone);}

	///Retruns the timescale for core-envelope coupling in Gyrs.
	double core_env_coupling_timescale() const
	{return __core_env_coupling_timescale;}

	///\brief Retuns the strength of the magnetic wind (the K constant).
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	double wind_strength() const {return __magnetic_wind_strength;}

	///Retruns the saturation frequency of the magnetic wind in rad/day.
	double wind_saturation_frequency() const
	{return __magnetic_wind_saturation_freq;}

	///\brief Returns the torque on the stellar envelope due to the stellar 
	///magnetic wind for the given age.
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
	///convective zone spin frequency at the given age.
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

	///\brief Partial derivative of the wind torque with respect to the
	///convective zone spin frequency at the current age.
	///
	///Units: \f$M_\odot \cdot R_\odot^2/\mathrm{Gyr}\f$
	double wind_torque_freq_deriv(
			///Spin frequency of the convective zone in rad/day.
			double conv_frequency,

			///The saturation state of the wind to assume. Could be UNKNOWN
			///(default), in which case it is determined by comparing the
			///convective spin frequency to the saturation frequency.
			WindSaturationState assume_wind_saturation=UNKNOWN) const
	{return wind_torque_freq_deriv(__age, conv_frequency,
			assume_wind_saturation);}

	/*
	///\brief Derivative of the wind torque on the stellar envelope with 
	///respect to stellar age at the given age.
	///
	///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$
	double wind_torque_age_deriv(double age, double conv_frequency) const;

	///\brief Derivative of the wind torque on the stellar envelope with 
	///respect to stellar age at the current age.
	///
	///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$
	double wind_torque_age_deriv(double conv_frequency) const
	{return wind_torque_age_deriv(__age, conv_frequency);}
	*/

	///\brief Partial derivative of the wind torque with respect to the
	///stellar age at the given age.
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

	///\brief Partial derivative of the wind torque with respect to the
	///stellar age at the current age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	double wind_torque_age_deriv(
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
			WindSaturationState assume_wind_saturation=UNKNOWN) const
	{return wind_torque_age_deriv(__age, angular_momentum,
			const_angular_momentum, assume_wind_saturation);}

	///\brief The wind torque at the given age in Gyrs.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///The angular momentum evolution for the convective zone should already
	///be specified by calling the set_angular_momentum_evolutio method with
	///zone=convective
	double wind_torque(double age) const
	{return wind_torque(age, spin_frequency(age, envelope));}

	///\brief The wind torque at the current age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///The angular momentum evolution for the convective zone should already
	///be specified by calling the set_angular_momentum_evolutio method with
	///zone=convective
	double wind_torque() const {return wind_torque(__age);}

	///\brief The torque on the stellar envelope due to the core-envelope
	///coupling at the given age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///The real part is the component along the convective angular momentum
	///and the imaginary part is 90 degrees away in the positive direction.
	std::complex<double> differential_rotation_torque_angmom(
			///Stellar age in Gyr.
			double age,
			
			///Angular momentum of the convective envelope in units of:
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv, 

			///Angular momentum of the radiative core in units of:
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			///The convention for the real and imaginary part is just like
			///for the return value.
			std::complex<double> Lrad) const;

	///\brief The torque on the stellar envelope due to the core-envelope
	///coupling at the current age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///The real part is the component along the convective angular momentum
	///and the imaginary part is 90 degrees away in the positive direction.
	std::complex<double> differential_rotation_torque_angmom(
			///Angular momentum of the convective envelope in units of:
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv, 

			///Angular momentum of the radiative core in units of:
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			///The convention for the real and imaginary part is just like
			///for the return value.
			std::complex<double> Lrad) const
	{return differential_rotation_torque_angmom(__age, Lconv, Lrad);}

	///\brief The partial derivative of the differential rotation torque at
	///the given age.
	///
	///Units: \f$Gyr^{-1}\f$ or 
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	///
	///Depending on the with_respect_to argument, the derivative is with
	///respect to one of the angular momenta or age, with the other two
	///variables held constant.
	std::complex<double> differential_rotation_torque_deriv(
			///The stellar age in Gyr.
			double age,
			
			///The convective zone angular momentum.
			double Lconv, 

			///The radiative zone angular momentum.
			std::complex<double> Lrad,
			
			///Identifies either the zone whose angular momentum we are
			///taking the pratial derivaite with respect to, if it is
			//(convective or radiative) or that it should be taken with
			///respect to the stellar age if with_respect_to is total.
			StellarZone with_respect_to=total) const;

	///\brief The partial derivative of the differential rotation torque at
	///the current age.
	///
	///See #differential_rotation_torque_deriv(double, double, double,
	///StellarZone) for details.
	std::complex<double> differential_rotation_torque_deriv(
			///The convective zone angular momentum.
			double Lconv, 

			///The radiative zone angular momentum.
			std::complex<double> Lrad,
			
			///Identifies either the zone whose angular momentum we are
			///taking the pratial derivaite with respect to, if it is
			//(convective or radiative) or that it should be taken with
			///respect to the stellar age if with_respect_to is total.
			StellarZone with_respect_to=total) const
	{return differential_rotation_torque_deriv(__age, Lconv, Lrad,
			with_respect_to);}

	///\brief The differential rotation torque on the stellar envelope for
	///the given amount of differential rotation at the given age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///The meaning of the real and complex parts is the same in
	///differential_rotation_torque_angmom.
	std::complex<double> differential_rotation_torque(
			///The stellar age in Gyr.
			double age, 

			///The amount of differential rotation in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			std::complex<double> differential_rotation_amount, 

			///The spin frequency of the convective zone in rad/day.
			double conv_frequency) const;

	///\brief The differential rotation torque on the stellar envelope for
	///the given amount of differential rotation at the current age.
	///
	///See #differential_rotation_torque(double, std::complex<double>,double)
	///for details
	std::complex<double> differential_rotation_torque(
			///The amount of differential rotation in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			std::complex<double> differential_rotation_amount, 

			///The spin frequency of the convective zone in rad/day.
			double conv_frequency) const
	{return differential_rotation_torque(__age, differential_rotation_amount,
			conv_frequency);}

	///\brief The partial derivative of the differential rotation torque on
	///the stellar envelope with respect to whatever the derivatives of the
	///differential rotation and convective frequency are.
	///
	///The units of the torque are \f$\frac{M_\odot \cdot R_\odot^2 \cdot
	/// \mathrm{rad}}{\mathrm{day} \cdot \mathrm{Gyr}}\f$
	std::complex<double> differential_rotation_torque_deriv(
			///The stellar ge in Gyr.
			double age, 

			///The amount of differential rotation in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			std::complex<double> differential_rotation_amount, 

			///The derivative of the differential rotation with respect to
			///whatever we want the torque differentiated against.
			std::complex<double> differential_rotation_deriv,

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

	///\brief The partial derivative of the differential rotation torque on
	///the stellar envelope with respect to whatever the derivatives of the
	///differential rotation and convective frequency are.
	///
	///The units of the torque are \f$\frac{M_\odot \cdot R_\odot^2 \cdot
	/// \mathrm{rad}}{\mathrm{day} \cdot \mathrm{Gyr}}\f$
	std::complex<double> differential_rotation_torque_deriv(
			///The amount of differential rotation in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			std::complex<double> differential_rotation_amount, 

			///The derivative of the differential rotation with respect to
			///whatever we want the torque differentiated against.
			std::complex<double> differential_rotation_deriv,

			///The spin frequency of the convective envelope in rad/day.
			double conv_frequency,

			///The derivative of the spin frequency of the convective
			///envelope with respect to whatever we want the torque 
			///differentiated against.
			double conv_frequency_deriv,
			
			///Set to true if the derivatives of the differential rotation
			///and the convective zone frequency are with respect to age.
			bool with_respect_to_age=false)
		const
	{return differential_rotation_torque_deriv(__age,
			differential_rotation_amount, differential_rotation_deriv,
			conv_frequency, conv_frequency_deriv, with_respect_to_age);}

	///\brief The torque on the stellar envelope due to the core-envelope
	///coupling for the given age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	///
	///Both the core and envelope angular momenta evolutions must be already
	///specified.
	std::complex<double> differential_rotation_torque(double age) const;

	///\brief The torque on the stellar envelope due to the core-envelope
	///coupling for the current age.
	///
	///See #differential_rotation_torque(double) for details.
	std::complex<double> differential_rotation_torque() const
	{return differential_rotation_torque(__age);}

	///\brief The amount of differential rotation between the convective
	///envelope and the radiative core for the given age.
	///
	///Units: \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
	///
	///The meaning of the real and imaginary parts is the same as for 
	///differential_rotation_torque_angmom
	std::complex<double> differential_rotation(
			///Stellar age in Gyr.
			double age,
			
			///Convective envelope angular momentum in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv,

			///Radiative core angular momentum in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$		
			std::complex<double> Lrad) const;

	///\brief The amount of differential rotation between the convective
	///envelope and the radiative core for the current age.
	///
	///See #differential_rotation(double, double, std::complex<double>) for
	///details.
	std::complex<double> differential_rotation(
			///Convective envelope angular momentum in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv,

			///Radiative core angular momentum in 
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$		
			std::complex<double> Lrad)  const
	{return differential_rotation(__age, Lconv, Lrad);}

	///\brief The partial derivative of the differential rotation between
	///the convective envelope and the radiative core at the given age.
	///
	///Units: dimensionless or \f$\frac{M_\odot \cdot R_\odot^2 \cdot
	/// \mathrm{rad}}{\mathrm{day}\cdot \mathrm{Gyr}}\f$
	///
	///Depending on the with_respect_to argument, the derivative is with
	///respect to one of the angular momenta or age, with the other two
	///variables held constant.
	std::complex<double> differential_rotation_deriv(
			///Stellar age in Gyr.
			double age,
			
			///Convective zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv, 

			///Radiative zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			std::complex<double> Lrad,
			
			///Identifies either the zone whose angular momentum we are
			///taking the pratial derivaite with respect to, if it is
			//(convective or radiative) or that it should be taken with
			///respect to the stellar age if with_respect_to is total.
			StellarZone with_respect_to=total) const;

	///\brief The partial derivative of the differential rotation between
	///the convective envelope and the radiative core at the current age.
	///
	///See #differential_rotation_deriv(double, double, double StellarZone)
	///for details.
	std::complex<double> differential_rotation_deriv(
			///Convective zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			double Lconv, 

			///Radiative zone angular momentum in
			/// \f$M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
			std::complex<double> Lrad,
			
			///Identifies either the zone whose angular momentum we are
			///taking the pratial derivaite with respect to, if it is
			//(convective or radiative) or that it should be taken with
			///respect to the stellar age if with_respect_to is total.
			StellarZone with_respect_to=total) const
	{return differential_rotation_deriv(__age, Lconv, Lrad,with_respect_to);}

	///\brief The amount of differential rotation between the convective
	///envelope and the radiative core at the given age.
	///
	///The angular momenta evolution of the two stellar zones must already be
	///specified.
	std::complex<double> differential_rotation(
			///Stellar age in Gyr.
			double age) const;

	///\brief The amount of differential rotation between the convective
	///envelope and the radiative core at the given age.
	///
	///The angular momenta evolution of the two stellar zones must already be
	///specified.
	std::complex<double> differential_rotation() const
	{return differential_rotation(__age);}

	///\brief The rate of moment of inertia transfer from the envolope to the
	///core due to the convective-radiative boundary evolving with time at
	///the given age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double core_inertia_gain(double age) const;

	///\brief The rate of moment of inertia transfer from the envolope to the
	///core due to the convective-radiative boundary evolving with time.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}}\f$
	double core_inertia_gain() const {return core_inertia_gain(__age);}

	///\brief The age derivative of the rate of moment of inertia transfer
	///from the envolope to the core due to the convective-radiative boundary
	///evolving with time at the given age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	double core_inertia_gain_deriv(double age) const;

	///\brief The age derivative of the rate of moment of inertia transfer
	///from the envolope to the core due to the convective-radiative boundary
	///evolving with time at the current age.
	///
	///Units: \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{rad}}
	/// {\mathrm{day} \cdot \mathrm{Gyr}^2}\f$
	double core_inertia_gain_deriv() const
	{return core_inertia_gain_deriv(__age);}

	///The age at which the star leaves the main sequence in Gyr.
	double lifetime() const {return __lifetime;}

	///Delete the allocated interpolated stellar properties.
	~StarBase();
};

#endif
