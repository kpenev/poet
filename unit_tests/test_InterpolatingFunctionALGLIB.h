/**\file
 *
 * \brief Declares the test suite that exercises the
 * InterpolatingFuctionALGLIB class and the parent classes.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_INTERPOLATING_FUNCTION_ALGLIB
#define __TEST_INTERPOLATING_FUNCTION_ALGLIB

#include "Common.h"
#include "../Functions.h"
#include <sstream>

/**\brief The test suite that exercises the
 * InterpolatingFuctionALGLIB class and the parent classes.
 *
 * \ingroup UnitTests_group
 */
class test_InterpolatingFunctionALGLIB : public Test::Suite {
private:
	///\brief Returns an array of the values of f evaluated at the given x
	///positions.
	///
	///The only requirement for f is that it must provide operator().
	template<class FUNC_TYPE>
	std::valarray<double> get_y(const FUNC_TYPE &f, 
			const std::valarray<double> &x);

	///\brief Checks the intersections of the given interpolating function
	///with a number of y values.
	void check_solutions(const InterpolatingFunctionALGLIB &interpf,
			const std::valarray<double> &check_crossings,
			const std::valarray< std::valarray<double> > &crossing_locations,
			const std::string &test_description,
			double frac_tolerance=1e-12, double abs_tolerance=1e-12);

	///Checks the derivatives of the given interpolating function.
	void check_deriv(const InterpolatingFunctionALGLIB &interpf,
			const std::valarray<double> &x, 
			const std::valarray<double> &answers,
			const std::valarray<double> &first_deriv,
			const std::valarray<double> &second_deriv,
			const std::string &test_description,
			double frac_tolerance=1e-12, double abs_tolerance=1e-12);

	///Checks the results of an interpolation against the prescribed answers.
	void check_interp(const InterpolatingFunctionALGLIB &interpf,
			const std::valarray<double> &x, 
			const std::valarray<double> &answers,
			const std::string &test_description,
			double frac_tolerance=1e-12, double abs_tolerance=1e-12);

	///\brief Tests the interpolation of a constant function without
	///specifying derivative.
	void test_const_no_deriv();

	///\brief Tests the interpolation of a linear function without specifying
	///derivative
	void test_linear_no_deriv();

	///\brief Tests the interpolation of a constant function with
	///specified(=0) derivative
	void test_const_deriv();

	///\brief Tests the interpolation of a linear function with
	///specified(=slope) derivative
	void test_linear_deriv();

	///\brief Tests the interpolation of a quadratic function with specified
	///derivative
	void test_quadratic_deriv();

	///\brief Tests the interpolation of a cubic function with specified
	///derivative
	void test_cubic_deriv();

	///Tests the approximate interpolation with smoothing.
	void test_sin_smoothed();

protected:
	///No fixtures at this time.
	void setup() {};

	///No fixtures at this time.
	void tear_down() {};
public:
	///Add all tests to the test suite.
	test_InterpolatingFunctionALGLIB();

	///Do nothing.
	~test_InterpolatingFunctionALGLIB() {}
};

#endif
