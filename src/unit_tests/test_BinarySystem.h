/**\file
 *
 * \brief Declares the test suite for exercising the BinarySystem class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_BINARY_SYSTEM_H
#define __TEST_BINARY_SYSTEM_H

#include "TwoZoneBody.h"
#include "Common.h"
#include "LaiExpressions.h"
#include "RandomDiskPlanetSystem.h"
#include "../BinarySystem.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctime>

/**\brief The test suite for the BinarySystem class.
 *
 * \ingroup UnitTests_group
 */
class test_BinarySystem : public Test::Suite {
private:
	///How many times to run random tests.
	unsigned __ntests;

	///\brief Tests if the orbit or differential equations returned by a 
	///system matches what is expected.
	void test_orbit_diff_eq(BinarySystem &system,
			const std::valarray<double> &expected,
			bool diff_eq=false);
protected:
	///No fixtures at this time
	void setup() {};

	///No fixtures at this time
	void tear_down() {};
public:
	///\brief Read eccentricity expansion coefficients from the given file
	///and add the tests.
	test_BinarySystem(
			///How many random configurations to initialize for tests.
			unsigned ntests,

			///The name of the file containing eccentricity expansion
			///coefficiets.
			const std::string &eccentricity_expansion=
			"eccentricity_expansion_coef.txt");

	///\brief Tests that the fill_orbit function works for single body with
	///locked surface rotation.
	void test_fill_orbit_locked_surface();

	///\brief Tests that the fill_orbit function works for single body with 
	///not-locked surface rotation.
	void test_fill_orbit_single();

	///\brief Tests that the fill_orbit function works for a binary with 
	///no locked zones.
	void test_fill_orbit_binary_no_locks();

	///\brief Tests that the fill_orbit function works for a binary with 
	///locked zones.
	void test_fill_orbit_binary_locks();

	///\brief Tests the differential equations for a surface locked single
	///body.
	void test_locked_surface_diff_eq();

	///\brief Tests the differential equations for a non-surface locked
	///single body with co-aligned zones.
	void test_single_aligned_diff_eq();

	///\brief Tests the differential equations for a non-surface locked
	///single body with all periapses=0.
	void test_single_zero_periapsis_diff_eq();

	///\brief Tests the differential equations for a binary with no locked 
	///zones, no eccentricity and all zones aligned with the orbit.
	void test_binary_no_locks_circular_aligned_diff_eq();

	///\brief Tests the differential equations for a banary with no locked
	///zones, no eccentricity but arbitrarily inclined zones with zero
	///periapsis.
	void test_binary_no_locks_circular_inclined_diff_eq();

	///Tests the differential equations for a banary with locked zones.
	void test_binary_locks_diff_eq();
};

#endif
