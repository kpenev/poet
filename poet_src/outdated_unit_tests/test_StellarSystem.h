/** \file
 *
 * \brief Defines the test suite and some related classes that exercises the
 * StellarSystem class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_STELLARSYSTEM_H
#define __TEST_STELLARSYSTEM_H

#include "Common.h"
#include "../StellarSystem.h"
#include <math.h>
#include <iostream>

///\brief Derives a bunch of quantities from the orbital parameters:
/// \f$a\f$, \f$L_{conv}\f$ and \f$L_{rad}\f$.
///
///All input quantities should be in the regular orbital evolution units and
///all derived quantities are in SI
///
///\ingroup UnitTests_group
class DeriveFromOrbit {
protected:
	///\brief Stellar convective zone moment of inertia in
	/// \f$\mathrm{kg}\cdot \mathrm{m}^2\f$.
	double __Iconv, 

		   ///\brief Stellar radiative zone angular momentum in
		   /// \f$\mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad/s}\f$.
		   __Lrad, 

		   __Mstar, ///< Stellar mass in kg
		   __Mplanet, ///< Planet mass in kg
		   __Rstar, ///< Stellar radius in m

		   ///\brief Wind strength constant in \f$ \mathrm{kg}\cdot
		   /// \mathrm{m}^2 \cdot\mathrm{s}/\mathrm{rad}^2 \f$
		   __windK, 


		   __wind_wsat, ///< Wind saturation frequency in rad/s
		   __coupling_timescale, ///< Core-envelope coupling timescale in s.

		   ///\brief Rate of change of the stellar convective zone angular
		   ///momentum in \f$\mathrm{kg}\cdot\mathrm{m}^2\mathrm{s}\f$
		   __dIconv_dt, 

		   ///\brief Stellar radiative zone moment of inertia in
		   /// \f$\mathrm{kg}\cdot \mathrm{m}^2\f$.
		   __Irad,

		   ///Stellar core radius in m.
		   __Rrad,

		   ///Stellar core mass in kg.
		   __Mrad,
		   
		   ///Rate of change of the stellar core mass in kg/s.
		   __Mrad_deriv,
		   
		   __a,///< Semimajor axis in m.

		   ///\brief Stellar convective zone angular momentum in 
		   /// \f$\mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad/s}\f$.
		   __Lconv,

		   __wconv, ///< Stellar convective zone spin frequency in rad/s.
		   __worb, ///< Orbital frequency in rad/s.

		   ///\brief Differential rotation amount in 
		   /// \f$\mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad/s}\f$.
		   __differential_rotation,

		   ///\brief Rate of transfer of angular momentum from the envelope
		   ///to the core due to the core growing in \f$\mathrm{kg}\cdot
		   /// \mathrm{m}^2 \cdot \mathrm{rad}/\mathrm{s}^2\f$.
		   __core_growth_torque,

		   ///The torque of the wind on the stellar convective zone in \f$
		   /// \mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad}/\mathrm{s}^2
		   /// \f$.
		   __wind_torque,
		   
		   ///The torque due to the coupling between the core and envelope on
		   ///the stellar convective zone in \f$\mathrm{kg}\cdot \mathrm{m}^2
		   /// \cdot \mathrm{rad}/\mathrm{s}^2\f$.
		   __coupling_torque;

	///\brief Calculates the quantities that are always calculated for all
	///sets of orbital parameters.
	void init_common(double age, double Lrad, const StellarSystem &system);
public:
	///In order to get a functioning object use operator()().
	DeriveFromOrbit() {}

	///Creates a fully functioning object.
	DeriveFromOrbit(double age, double a, double Lconv, double Lrad,
			const StellarSystem &system)
	{operator()(age, a, Lconv, Lrad, system);}

	///\brief For a default constructed object sets the orbit to derive the
	///quantities from.
	void operator()(double age, double a, double Lconv, double Lrad,
			const StellarSystem &system);

	///Semimajor axis in m.
	double a() const {return __a;}

	///\brief Stellar convective zone moment of inertia in
	/// \f$\mathrm{kg}\cdot \mathrm{m}^2\f$.
	double Iconv() const {return __Iconv;}

	///\brief Stellar radiative zone moment of inertia in
	/// \f$\mathrm{kg}\cdot \mathrm{m}^2\f$.
	double Irad() const {return __Irad;}

	///\brief Stellar convective zone angular momentum in 
	/// \f$\mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad/s}\f$.
	double Lconv() const {return __Lconv;}

	///\brief Stellar radiative zone angular momentum in
	/// \f$\mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad/s}\f$.
	double Lrad() const {return __Lrad;}

	///Stellar convective zone spin frequency in rad/s.
	double wconv() const {return __wconv;}

	///Orbital frequency in rad/s.
	double worb() const {return __worb;}

	///Stellar core radius in m.
	double Rrad() const {return __Rrad;}

	///Stellar mass in kg
	double Mstar() const {return __Mstar;}

	///Planet mass in kg
	double Mplanet() const {return __Mplanet;}

	///Stellar radius in m
	double Rstar() const {return __Rstar;}

	///\brief Differential rotation amount in 
	/// \f$\mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad/s}\f$.
	double differential_rotation() const {return __differential_rotation;}

	///Stellar core mass in kg.
	double Mrad() const {return __Mrad;}

	///Rate of change of the stellar core mass in kg/s.
	double Mrad_deriv() const {return __Mrad_deriv;}

	///\brief Rate of transfer of angular momentum from the envelope
	///to the core due to the core growing in \f$\mathrm{kg}\cdot
	/// \mathrm{m}^2 \cdot \mathrm{rad}/\mathrm{s}^2\f$.
	double core_growth_torque() const {return __core_growth_torque;}

	///\brief Wind strength constant in \f$ \mathrm{kg}\cdot
	/// \mathrm{m}^2 \cdot\mathrm{s}/\mathrm{rad}^2 \f$
	double windK() const {return __windK;}

	///Wind saturation frequency in rad/s
	double wind_wsat() const {return __wind_wsat;}

	///\brief The torque of the wind on the stellar convective zone in \f$
	/// \mathrm{kg}\cdot \mathrm{m}^2 \cdot \mathrm{rad}/\mathrm{s}^2
	/// \f$.
	double wind_torque() const {return __wind_torque;}

	///Core-envelope coupling timescale in s.
	double coupling_timescale() const {return __coupling_timescale;}

	///\brief The torque due to the coupling between the core and envelope on
	///the stellar convective zone in \f$\mathrm{kg}\cdot \mathrm{m}^2
	/// \cdot \mathrm{rad}/\mathrm{s}^2\f$.
	double coupling_torque() const {return __coupling_torque;}

	///\brief Rate of change of the stellar convective zone angular
	///momentum in \f$\mathrm{kg}\cdot\mathrm{m}^2\mathrm{s}\f$
	double dIconv_dt() const {return __dIconv_dt;}
};

///\brief Derives a bunch of quantities from the orbital parameters:
/// \f$a\f$ and \f$L_{rad}\f$ assuming a spin-orbit lock.
///
///All input quantities should be in the regular orbital evolution units and
///all derived quantities are in SI
///
///\ingroup UnitTests_group
class DeriveFromLockedOrbit : public DeriveFromOrbit {
public:
	///Creates a fully functioning object.
	DeriveFromLockedOrbit(double age, double a, double Lrad,
			const StellarSystem &system)
	{operator()(age, a, Lrad, system);}
	
	///\brief For a default constructed object sets the orbit to derive the
	///quantities from.
	void operator()(double age, double a, double Lrad,
			const StellarSystem &system);
};

///\brief Derives a bunch of quantities from the orbital parameters:
/// \f$L_{conv}\f$ and \f$L_{rad}\f$ for NO_PLANET evolution mode.
///
///All input quantities should be in the regular orbital evolution units and
///all derived quantities are in SI
///
///\ingroup UnitTests_group
class DeriveFromNoPlanet : public DeriveFromOrbit {
public:
	///Creates a fully functioning object.
	DeriveFromNoPlanet(double age, double Lconv, double Lrad,
			const StellarSystem &system)
	{operator()(age, Lconv, Lrad, system);}

	///\brief For a default constructed object sets the orbit to derive the
	///quantities from.
	void operator()(double age, double Lconv, double Lrad,
			const StellarSystem &system);
};

///\brief Derives a bunch of quantities from \f$L_{rad}\f$.
///
///All input quantities should be in the regular orbital evolution units and
///all derived quantities are in SI
///
///\ingroup UnitTests_group
class DeriveFromDiskLocked : public DeriveFromOrbit {
public:
	///Creates a fully functioning object.
	DeriveFromDiskLocked(double age, double Lrad,
			const StellarSystem &system)
	{operator()(age, Lrad, system);}

	///\brief For a default constructed object sets the orbit to derive the
	///quantities from.
	void operator()(double age, double Lrad, const StellarSystem &system);
};

///\brief The test suite that exercises the StellarSystem class.
///
///\ingroup UnitTests_group
class test_StellarSystem : public Test::Suite {
private:
	SystemData* ssdata;
	StarData* sdata;
	PlanetData* pdata;
	StellarSystem* system;

	///\brief Returns the expected differential equation for the given orbit
	///for the case of a planet around a star whose rotation is not locked to
	///the orbit.
	void predict_unlocked_orbit_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs, int dadt_sign=0);

	///\brief Returns the expected differential equation for the given orbit
	///for the case of a planet whose orbit is locked to the rotation of the
	///star.
	void predict_locked_orbit_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs);
	
	///\brief Returns the expected differential equation for the rotational
	///evolution of a star with no planet around it.
	void predict_no_planet_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs);

	///\brief Returns the torque on the core of the star when the convective
	///zone is locked to a disk.
	void predict_disk_locked_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs);
protected:
	///Create a random stellar system.
	void setup();

	///Cleanup the random stellar system created for the test.
	void tear_down();

public:
	///Create a test suite adding all the unit tests.
	test_StellarSystem();

	///A no-test
	void test_basic();

	///Tests the radius interpolation
	void test_interpolate();

	///\brief Tests the equation for the evolution of a^6.5, Lconv and Lrad
	///under the assumption that no lock to the orbit or the disk exists.
	void test_unlocked_orbit();

	///\brief Tests the equation for the evolution of a and Lrad under the
	///assumption that the orbit is locked to the stellar rotation
	void test_locked_orbit();

	///\brief Tests the equation for the evolution of Lconv and Lrad for a
	///single star with no planet around it.
	void test_no_planet_orbit();

	///\brief Tests the equation for the evolution of the core angular
	///momentum for a star whose convective zone is locked to a disk.
	void test_disk_locked_orbit();

	///Tests the jacobian of the unlocked orbit differential equation.
	void test_unlocked_jacobian();

	///Tests the jacobian of the locked orbit differential equation.
	void test_locked_jacobian();

	///\brief Tests the jacobian of the planet-free stellar rotation
	///differential equation.
	void test_no_planet_jacobian();

	///Tests the jacobian of the disk locked evolution
	void test_disk_locked_jacobian();

//	void test_equations();

	///\brief Estimates the part of the jacobian that corresponds to the
	///given index.
	void check_param_derivs(double age, double* params, int index,
			double* jacobian);

	///Estimates the full jacobian.
	void eval_jacobian(double age, double* params, double* jacobian,
			long double* ageDeriv);

	///Tests the jacobian for the FAST_PLANET/SLOW_PLANET evolution modes.
	void test_jacobian();
};

#endif /* TEST_STELLARSYSTEM_H_ */
