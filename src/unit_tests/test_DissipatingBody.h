/**\file 
 *
 * \brief Declares the test suite for exercising the DissipatingBody class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_DISSIPATING_BODY_H
#define __TEST_DISSIPATING_BODY_H

#include "TwoZoneBody.h"
#include "Common.h"
#include "ConstPhaseLagDissipatingZone.h"
#include "LaiExpressions.h"
#include "../DissipatingBody.h"
#include "../OrbitalExpressions.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctime>

/**\brief The test suite for the DissipatingBody class.
 *
 * \ingroup UnitTests_group
 */
class test_DissipatingBody : public Test::Suite {
private:
	///How many random configurations to initialize for tests.
	unsigned __ntests;

	///Generates a randomly configured body.
	TwoZoneBody *random_body(double &other_mass, double &a, Lags &lags_env,
							 Lags &lags_core, bool no_periapsis=false,
							 bool same_inclination=false) const;
protected:
	///No fixtures at this time
	void setup() {};

	///No fixtures at this time
	void tear_down() {};
public:
	///\brief Read eccentricity expansion coefficients from the given file
	///and add the tests.
	test_DissipatingBody(
			///How many random configurations to initialize for tests.
			unsigned ntests,

			///The name of the file containing eccentricity expansion
			///coefficiets.
			const std::string &eccentricity_expansion=
			"eccentricity_expansion_coef.txt");

	///\brief Tests the normalized torque and power expressions against eq.
	///(27), (28) and (35) in Lai 2012.
	void test_Lai_torque_power();

	///\brief Tests the energy and angular momentum rates of change
	///of the orbit due to tides for a 2-zone body with no periapsis
	///difference.
	void test_orbit_rates_same_periapsis();
};

#endif
