/**\file
 *
 * \brief Declares the test suite for the DissipatingZone class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_DISSIPATING_ZONE_H
#define __TEST_DISSIPATING_ZONE_H

#include "Common.h"
#include "../DissipatingZone.h"
#include "../OrbitalExpressions.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <map>
#include <utility>
#include <cassert>
#include <functional>

bool operator<(const std::pair<int, int> &v1, const std::pair<int, int> &v2)
{
	return v1.first<v2.first || (v1.first==v2.first && v1.second<v2.second);
}

class Lags : public std::map< std::pair<int, int>, double> {
public:
	///A reference to the requested element (creates it if it does not
	///exist).
	double &operator()(int m, int mp);

	///A copy of the requested element (error if it does not exist).
	double operator()(int m, int mp) const;
};

class TestingDissipatingZone : public DissipatingZone {
public:
	///Should place enough informatino to identify the zone to the given
	///stream.
	virtual void describe(std::ostream &os) const=0;
};

class ConstPhaseLagDissipatingZone : public TestingDissipatingZone {
private:
	Lags __lags;

	double __moment_of_inertia, __radius, __mass;
public:
	ConstPhaseLagDissipatingZone(const Lags &lags, double moment_of_inertia,
			double radius, double mass)
		: __lags(lags), __moment_of_inertia(moment_of_inertia),
		__radius(radius), __mass(mass) {}

	///\brief Should return true iff the given term is presently locked.
	virtual bool locked(int orbital_frequency_multiplier,
			int spin_frequency_multiplier) const {return false;}

	///\brief Should return the tidal phase lag time the love number for the
	///given tidal term (or one of its derivatives).
	///
	///In case the specified term is in a lock, it should return the phase
	///lag for the case of the spin frequency approaching the lock from
	///below. The lag for spin frequency approaching from above should be
	///written to above_lock_value. If the term is not locked 
	///leave above_lock_value untouched.
	virtual double modified_phase_lag(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,
			
			///The current orbital spin frequency in rad/day
			double orbital_frequency,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::Derivative deriv,

			///If the lag of a locked term is calculated this should be set
			///to the lag assuming the spin frequency is just above the lock.
			///Otherwise, leave untouched.
			double &above_lock_value) const
	{
		return (deriv!=Dissipation::NO_DERIV ? 0
				: __lags(spin_frequency_multiplier,
					orbital_frequency_multiplier));
	}

	double moment_of_inertia(int deriv) const
	{return (deriv ? 0 : __moment_of_inertia);}

	double love_coefficient(int, int, Dissipation::Derivative) const
	{return 0;}

	double outer_radius(int deriv) const {return (deriv ? 0 : __radius);}
	double outer_mass(int deriv) const {return (deriv ? 0 : __mass);}

	void describe(std::ostream &os) const;
};

/**\brief The test suite for the DissipatingZone class.
 *
 * \ingroup UnitTests_group
 */
class test_DissipatingZone : public Test::Suite {
private:
	///Did the currently run test fail?
	bool __current_failed;

	///\brief Returns the tidal dissipation torque along z according to
	///Lai 2012 Eq. (27).
	double torque_z_Lai(
			///The inclination angle in radians.
			double inclination,

			const Lags &lags);

	///\brief Returns the tidal dissipation torque along x according to 
	///Lai 2012  Eq. (35).
	double torque_x_Lai(
			///The inclination angle in radians.
			double inclination,

			const Lags &lags);

	///Returns the tidal dissipation power according to Lai 2012 Eq. (28).
	double power_Lai(
			///The inclination angle in radians.
			double inclination,

			const Lags &lags);

	///Performs a single check of the expected values of the tidal torques
	///and power.
	void single_test(
			///The zone to test. Modified because of set_orbit.
			TestingDissipatingZone &zone,

			///The orbital freuqency in rad/day
			double orbital_frequency,

			///The eccentricity of the orbit
			double eccentricity,

			///The mass of the primary.
			double m1,

			///The mass of the secondary.
			double m2,

			///The spin angular momentum of the zone.
			double spin_angmom,

			///The inclination to assume for the zone.
			double inclination,

			///The periapsis to assume for the zone.
			double periapsis,
			
			///The expected tidal torque in the z direction
			double torque_z,
			
			///The expected tidal torque in the x direction
			double torque_x,
			
			///The expected tidal power.
			double power);

protected:
	///No fixtures at this time
	void setup() {};

	///No fixtures at this time
	void tear_down() {};
public:
	///\brief Read eccentricity expansion coefficients from the given file
	///and add the tests.
	test_DissipatingZone(
			const std::string &eccentricity_expansion=
			"eccentricity_expansion_coef.txt");

	///\brief Tests the tidal torque and power against Eq. (27), (28) and
	///(35) in Lai 2012.
	void test_Lai();

	///\brief Tests the derivatives of the periapsis and inclination 
	///evolution rates.
	void test_inclination_periapsis_evol_deriv();
};

#endif
