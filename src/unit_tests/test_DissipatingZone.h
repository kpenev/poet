/**\file
 *
 * \brief Declares the test suite for the DissipatingZone class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_DISSIPATING_ZONE_H
#define __TEST_DISSIPATING_ZONE_H

#include "Common.h"
#include "ConstPhaseLagDissipatingZone.h"
#include "LaiExpressions.h"
#include "../DissipatingZone.h"
#include "../OrbitalExpressions.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <cassert>
#include <functional>
#include <ctime>

/**\brief The test suite for the DissipatingZone class.
 *
 * \ingroup UnitTests_group
 */
class test_DissipatingZone : public Test::Suite {
private:
	///Did the currently run test fail?
	bool __current_failed;

	///Performs a single check of the expected values of the tidal torques
	///and power.
	void single_test(
			///The eccentricity expansion order to use for the test (if
			///non-zero, eccentricity itself must be zero).
			int e_order,

			///The tidal lags to assume for the test.
			const Lags &lags,

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
