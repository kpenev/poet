/*
 * test_StellarSystem.h
 *
 *  Created on: Nov 20, 2012
 *      Author: stanley
 */

#ifndef __TEST_STELLARSYSTEM_H
#define __TEST_STELLARSYSTEM_H

#include "Common.h"
#include "../StellarSystem.h"
#include <math.h>
#include <iostream>

///Derives a bunch of quantities from the orbital parameters: a, Lconv, Lrad.
///All input quantities should be in the regular orbital evolution units and
///all derived quantities are in SI
class DeriveFromOrbit {
protected:
	//initialized by init_common
	double __Iconv, __Lrad, __Mstar, __Mplanet, __Rstar, __windK,
		   __wind_wsat, __coupling_timescale, __dIconv_dt, __Irad, __Rrad,
		   __Mrad, __Mrad_deriv,
		   
		   //Not initialized by init_common
		   __a, __Lconv, __wconv, __worb, __differential_rotation,
		   __core_growth_torque, __wind_torque, __coupling_torque;

	///Calculates the quantities that are always calculated for all sets of
	///orbital parameters.
	void init_common(double age, double Lrad, const StellarSystem &system);
public:
	DeriveFromOrbit() {}

	DeriveFromOrbit(double age, double a, double Lconv, double Lrad,
			const StellarSystem &system)
	{operator()(age, a, Lconv, Lrad, system);}

	void operator()(double age, double a, double Lconv, double Lrad,
			const StellarSystem &system);

	double a() const {return __a;}
	double Iconv() const {return __Iconv;}
	double Irad() const {return __Irad;}
	double Lconv() const {return __Lconv;}
	double Lrad() const {return __Lrad;}
	double wconv() const {return __wconv;}
	double worb() const {return __worb;}
	double Rrad() const {return __Rrad;}
	double Mstar() const {return __Mstar;}
	double Mplanet() const {return __Mplanet;}
	double Rstar() const {return __Rstar;}
	double differential_rotation() const {return __differential_rotation;}
	double Mrad() const {return __Mrad;}
	double Mrad_deriv() const {return __Mrad_deriv;}
	double core_growth_torque() const {return __core_growth_torque;}
	double windK() const {return __windK;}
	double wind_wsat() const {return __wind_wsat;}
	double wind_torque() const {return __wind_torque;}
	double coupling_timescale() const {return __coupling_timescale;}
	double coupling_torque() const {return __coupling_torque;}
	double dIconv_dt() const {return __dIconv_dt;}
};

class DeriveFromLockedOrbit : public DeriveFromOrbit {
public:
	DeriveFromLockedOrbit(double age, double a, double Lrad,
			const StellarSystem &system)
	{operator()(age, a, Lrad, system);}
	
	void operator()(double age, double a, double Lrad,
			const StellarSystem &system);
};

class DeriveFromNoPlanet : public DeriveFromOrbit {
public:
	DeriveFromNoPlanet(double age, double Lconv, double Lrad,
			const StellarSystem &system)
	{operator()(age, Lconv, Lrad, system);}

	void operator()(double age, double Lconv, double Lrad,
			const StellarSystem &system);
};

class DeriveFromDiskLocked : public DeriveFromOrbit {
public:
	DeriveFromDiskLocked(double age, double Lrad,
			const StellarSystem &system)
	{operator()(age, Lrad, system);}

	void operator()(double age, double Lrad, const StellarSystem &system);
};

class test_StellarSystem : public Test::Suite {
private:
	SystemData* ssdata;
	StarData* sdata;
	PlanetData* pdata;
	StellarSystem* system;

	///Returns the expected differential equation for the given orbit for
	///the case of a planet around a star whose rotation is not locked to the
	///orbit.
	void predict_unlocked_orbit_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs, int dadt_sign=0);

	///Returns the expected differential equation for the given orbit for the
	///case of a planet whose orbit is locked to the rotation of the star.
	void predict_locked_orbit_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs);
	
	///Returns the expected differential equation for the rotational evolution
	///of a star with no planet around it.
	void predict_no_planet_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs);

	///Returns the torque on the core of the star when the convective zone is
	///locked to a disk.
	void predict_disk_locked_diff_eq(const DeriveFromOrbit &quantities_SI,
			double *test_derivs);
protected:
	///Create a random stellar system.
	void setup();

	///Cleanup the random stellar system created for the test.
	void tear_down();

public:
	test_StellarSystem();

	///A no-test
	void test_basic();

	void test_interpolate();

	///Tests the equation for the evolution of a^6.5, Lconv and Lrad under
	///the assumption that no lock to the orbit or the disk exists.
	void test_unlocked_orbit();

	///Tests the equation for the evolution of a and Lrad under the 
	///assumption that the orbit is locked to the stellar rotation
	void test_locked_orbit();

	///Tests the equation for the evolution of Lconv and Lrad for a single
	///star with no planet around it.
	void test_no_planet_orbit();

	///Tests the equation for the evolution of the core angular momentum for
	///a star whose convective zone is locked to a disk.
	void test_disk_locked_orbit();

	///Tests the jacobian of the unlocked orbit differential equation.
	void test_unlocked_jacobian();

	///Tests the jacobian of the locked orbit differential equation.
	void test_locked_jacobian();

	///Tests the jacobian of the planet-free stellar rotation differential
	///equation.
	void test_no_planet_jacobian();

	///Tests the jacobian of the disk locked evolution
	void test_disk_locked_jacobian();

//	void test_equations();

	void check_param_derivs(double age, double* params, int index,
			double* jacobian);

	void eval_jacobian(double age, double* params, double* jacobian,
			long double* ageDeriv);

	void test_jacobian();
};

#endif /* TEST_STELLARSYSTEM_H_ */
