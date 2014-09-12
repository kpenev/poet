/**\file
 *
 * \brief Declares the test suite for the EccetricityExpansionCoefficients
 * class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_ECCENTRICITY_EXPANSION_COEFFICIENTS_H
#define __TEST_ECCENTRICITY_EXPANSION_COEFFICIENTS_H

#include "Common.h"
#include "../EccentricityExpansionCoefficients.h"
#include <iostream>
#include <vector>
#include <sstream>

/**\brief The test suite for the EccentricityExpansionCoefficients class.
 *
 * \ingroup UnitTests_group
 */
class test_EccentricityExpansionCoefficients : public Test::Suite {
private:
	///The series approximations to test.
	EccentricityExpansionCoefficients __pms;

	///Performs a single test against a known value.
	void single_test(
			///The first index (0 or +-2).
			int m, 

			///The second index.
			int s,

			///The value of the eccentricity to use.
			double e,

			///The maximum eccentricity order to include in the Taylor
			///series.
			unsigned max_e_power,
			
			///If true the result is differentiated w.r.t. to the
			///eccentricity.
			bool deriv,

			///The expected answer
			double answer);
protected:
	///No fixtures at this time
	void setup() {};

	///No fixtures at this time
	void tear_down() {};

public:
	///Read expansion coefficients from the given file and add the tests.
	test_EccentricityExpansionCoefficients(
			const std::string &fname="eccentricity_expansion_coef.txt");

	///\brief Compares the non-derivative values of expansion coefficients
	///against a known set.
	void test_values();

	///\brief Compares the derivatives of expansion coefficients against a
	///known set.
	void test_derivatives();
};

#endif
