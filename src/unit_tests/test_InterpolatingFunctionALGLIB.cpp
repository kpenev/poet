/**\file
 *
 * \brief Defines some of the methods of the test suite that exercises the
 * InterpolatingFuctionALGLIB class and the parent classes.
 *
 * \ingroup UnitTests_group
 */

#include "test_InterpolatingFunctionALGLIB.h"
#include <fstream>
#include <iomanip>

template<class FUNC_TYPE>
std::valarray<double> get_y(const FUNC_TYPE &f, 
		const std::valarray<double> &x)
{
	std::valarray<double> y(x.size());
	for(size_t i=0; i<x.size(); i++) y[i]=f(x[i]);
	return y;
}

void test_InterpolatingFunctionALGLIB::check_solutions(
		const InterpolatingFunctionALGLIB &interpf,
		const std::valarray<double> &check_crossings,
		const std::valarray< std::valarray<double> > &crossing_locations,
		const std::string &test_description,
		double frac_tolerance, double abs_tolerance)
{
	for(size_t y_i=0; y_i<check_crossings.size(); y_i++) {
		double y=check_crossings[y_i];
		if(crossing_locations[y_i].size()==1 && 
				std::isnan(crossing_locations[y_i][0])) {
			TEST_THROWS_MSG(interpf.crossings(y), std::domain_error,
					(test_description+
					 " no domain error exception thrown!").c_str())
			continue;
		}
		try {
			InterpSolutionIterator crossings=interpf.crossings(y);
			for(size_t c_i=0; c_i<crossing_locations[y_i].size(); c_i++){
				std::ostringstream test_msg;
				double expected=crossing_locations[y_i][c_i];
				test_msg << test_description << " crossing #" << c_i
					<< " with y=" << y << " expected at x=" << expected;
				TEST_ASSERT_MSG(!(crossings.out_of_range()), 
						(test_msg.str()+" not found!").c_str())
					test_msg << " found at x=" << *crossings << " (diff="
					<< expected-*crossings << ")!";
				TEST_ASSERT_MSG(check_diff(expected, *crossings, 
							frac_tolerance, abs_tolerance),
						test_msg.str().c_str())
					crossings++;
			}
			if(!crossings.out_of_range()) {
				size_t sol_ind=crossing_locations[y_i].size();
				std::ostringstream test_msg;
				test_msg << test_description << " with y=" << y 
					<< " extra solutions found:";
				for(;!crossings.out_of_range(); crossings++)
					test_msg << "#" << sol_ind++ << ": " << *crossings 
						<< ", ";
				TEST_ASSERT_MSG(false, test_msg.str().c_str())
			}
		} catch(...) {
			TEST_ASSERT_MSG(false, 
				(test_description+": unexpected exception thrown.").c_str())
		}
	}
}

void test_InterpolatingFunctionALGLIB::check_deriv(
		const InterpolatingFunctionALGLIB &interpf,
		const std::valarray<double> &x, 
		const std::valarray<double> &answers,
		const std::valarray<double> &first_deriv,
		const std::valarray<double> &second_deriv,
		const std::string &test_description,
		double frac_tolerance, double abs_tolerance)
{
	for(size_t i=0; i<x.size(); i++) {
		std::ostringstream test_msg;
		const CubicSplineDerivatives *deriv=interpf.deriv(x[i]);
		double ans=answers[i];
		test_msg << test_description << " Expected value: " << ans 
			<< " got: " << deriv->order(0) << std::endl;
		if(answers.size())
			TEST_ASSERT_MSG(check_diff(ans, deriv->order(0), frac_tolerance,
						abs_tolerance), test_msg.str().c_str())

		if(first_deriv.size()) {
			test_msg.clear();
			test_msg.str("");
			test_msg << test_description << " Expected first derivative: " 
				<< first_deriv[i] << " got: " << deriv->order(1) << std::endl;
			TEST_ASSERT_MSG(check_diff(first_deriv[i], deriv->order(1), 
					frac_tolerance, abs_tolerance), test_msg.str().c_str())
		}
		if(second_deriv.size()) {
			test_msg.clear();
			test_msg.str("");
			test_msg << test_description << " Expected second derivative: " 
				<< second_deriv[i] << " got: " << deriv->order(2) << std::endl;
			TEST_ASSERT_MSG(check_diff(second_deriv[i], deriv->order(2), 
					frac_tolerance, abs_tolerance), test_msg.str().c_str())
		}
		test_msg.clear();
		test_msg.str("");
		test_msg << test_description 
			<< " Expected third derivative=NaN got: " << deriv->order(3) 
			<< std::endl;
		TEST_ASSERT_MSG(std::isnan(deriv->order(3)), test_msg.str().c_str())
		test_msg.clear();
		test_msg.str("");
		test_msg << test_description << " Expected fourth derivative=0 got: "
			<< deriv->order(4) << std::endl;
		TEST_ASSERT_MSG(deriv->order(4)==0, test_msg.str().c_str())
	}
}

void test_InterpolatingFunctionALGLIB::check_interp(
		const InterpolatingFunctionALGLIB &interpf,
		const std::valarray<double> &x, 
		const std::valarray<double> &answers,
		const std::string &test_description,
		double frac_tolerance, double abs_tolerance)
{
	for(size_t i=0; i<x.size(); i++) {
		std::ostringstream test_msg;
		double ans=answers[i], interp_val=interpf(x[i]);
		test_msg << test_description << " Expected value: " << ans
			<< " got: " << interp_val << std::endl;
		TEST_ASSERT_MSG(check_diff(ans, interp_val, frac_tolerance, 
					abs_tolerance), test_msg.str().c_str())
	}
}

void test_InterpolatingFunctionALGLIB::test_const_no_deriv()
{
	std::valarray<double> x(5), y(0.0, x.size()), x_check(10),
		crossing_checks(5);
	std::valarray< std::valarray<double> > 
		no_crossing_locations(std::valarray<double>(),
				crossing_checks.size());
	for(size_t i=0; i<x.size(); i++) 
		x[i]=static_cast<double>(i)/(x.size()-1);
	for(size_t i=0; i<x_check.size(); i++)
		x_check[i]=static_cast<double>(i)/(x_check.size()-1);
	for(int two_y=-2; two_y<=2; two_y++) 
		crossing_checks[two_y+2]=0.5*two_y;
	no_crossing_locations[2].resize(1, NaN);
	InterpolatingFunctionALGLIB interpf(x,y);
	check_interp(interpf, x_check, 
			std::valarray<double>(0.0, x_check.size()), 
			"Checking the interpolation of a constant(=0) function (no "
			"derivatives specified).");
	check_deriv(interpf, x_check, 
			std::valarray<double>(0.0, x_check.size()), 
			std::valarray<double>(0.0, x_check.size()),
			std::valarray<double>(0.0, x_check.size()),
			"Checking the derivatives of a constant(=0) function (no "
			"derivatives specified).");
	check_solutions(interpf, crossing_checks, no_crossing_locations,
			"Checking the solutions of a constant(=0) function (no "
			"derivatives specified).");
	no_crossing_locations[2].resize(0);
	y=1.0;
	no_crossing_locations[4].resize(1, NaN);
	interpf=InterpolatingFunctionALGLIB(x,y);
	check_interp(interpf, x_check, 
			std::valarray<double>(1.0, x_check.size()), 
			"Checking the interpolation of a constant(=1) function (no "
			"derivatives specified).");
	check_deriv(interpf, x_check, 
			std::valarray<double>(1.0, x_check.size()), 
			std::valarray<double>(0.0, x_check.size()),
			std::valarray<double>(0.0, x_check.size()),
			"Checking the derivatives of a constant(=1) function (no "
			"derivatives specified).");
	check_solutions(interpf, crossing_checks, no_crossing_locations,
			"Checking the solutions of a constant(=1) function (no "
			"derivatives specified).");
	no_crossing_locations[4].resize(0);
	y=0.5;
	no_crossing_locations[3].resize(1, NaN);
	interpf=InterpolatingFunctionALGLIB(x,y);
	check_interp(interpf, x_check, 
			std::valarray<double>(0.5, x_check.size()), 
			"Checking the interpolation of a constant(=1/2) function (no "
			"derivatives specified).");
	check_deriv(interpf, x_check, 
			std::valarray<double>(0.5, x_check.size()), 
			std::valarray<double>(0.0, x_check.size()),
			std::valarray<double>(0.0, x_check.size()),
			"Checking the derivatives of a constant(=1/2) function (no "
			"derivatives specified).");
	check_solutions(interpf, crossing_checks, no_crossing_locations,
			"Checking the solutions of a constant(=1/2) function (no "
			"derivatives specified).");
}

void test_InterpolatingFunctionALGLIB::test_linear_no_deriv()
{
	for(int two_slope=-2; two_slope<=2; two_slope++) {
		double slope=0.5*two_slope;
		for(int two_zeropt=-2; two_zeropt<=2; two_zeropt++) {
			double zeropt=0.5*two_zeropt;
			std::valarray<double> x(5), y(x.size()), x_check(10), 
				y_check(x_check.size()), crossing_checks(5);
			std::valarray< std::valarray<double> > 
				crossing_locations(5);
			for(int two_y=-2; two_y<=2; two_y++) {
				int ind=two_y+2;
				double y=0.5*two_y;
				crossing_checks[ind]=y;
				double cross=(y-zeropt)/slope;
				if((cross>=0 && cross<=1) || std::isnan(cross)) 
					crossing_locations[ind].resize(1,cross);
			}
			for(size_t i=0; i<x.size(); i++) {
				x[i]=static_cast<double>(i)/(x.size()-1);
				y[i]=zeropt+slope*x[i];
			}
			for(size_t i=0; i<x_check.size(); i++) {
				x_check[i]=static_cast<double>(i)/(x_check.size()-1);
				y_check[i]=zeropt+slope*x_check[i];
			}
			std::ostringstream test_descr;
			test_descr << "Checking the interpolation of y="
				<< slope << "*x + " << zeropt 
				<< " with no derivatives specified"; 
			InterpolatingFunctionALGLIB interpf(x,y);
			check_interp(interpf, x_check, y_check,
					test_descr.str());
			test_descr.str("");
			test_descr << "Checking the derivatives of y="
				<< slope << "*x + " << zeropt 
				<< " with no derivatives specified"; 
			check_deriv(interpf, x_check, y_check,
					std::valarray<double>(slope, x_check.size()),
					std::valarray<double>(0.0, x_check.size()),
					test_descr.str());
			test_descr.str("");
			test_descr << "Checking the solutions of y="
				<< slope << "*x + " << zeropt 
				<< " with no derivatives specified"; 
			check_solutions(interpf, crossing_checks, crossing_locations,
					test_descr.str(), 1e-8, 1e-8);
		}
	}
}

void test_InterpolatingFunctionALGLIB::test_const_deriv()
{
	for(int two_value=-2; two_value<=2; two_value++) {
		double value=0.5*two_value;
		std::valarray<double> x(5), x_check(10), crossing_checks(5);
		std::valarray< std::valarray<double> > 
			no_crossing_locations(std::valarray<double>(),
					crossing_checks.size());
		for(int two_y=-2; two_y<=2; two_y++) 
			crossing_checks[two_y+2]=0.5*two_y;
		for(size_t i=0; i<x.size(); i++) 
			x[i]=static_cast<double>(i)/(x.size()-1);
		for(size_t i=0; i<x_check.size(); i++) 
			x_check[i]=static_cast<double>(i)/(x_check.size()-1);
		std::ostringstream test_descr;
		test_descr << "Checking the interpolation of a constant(="
			<< value << ") function with zero derivatives specified";
		InterpolatingFunctionALGLIB interpf(x, 
					std::valarray<double>(value, x.size()), 
					std::valarray<double>(0.0, x.size()));
		check_interp(interpf, x_check,
				std::valarray<double>(value, x_check.size()),
				test_descr.str());
		test_descr.str("");
		test_descr << "Checking the derivatives of a constant(="
			<< value << ") function with zero derivatives specified";
		check_deriv(interpf, x_check, 
				std::valarray<double>(value, x_check.size()),
				std::valarray<double>(0.0, x_check.size()),
				std::valarray<double>(0.0, x_check.size()),
				test_descr.str());
		test_descr.str("");
		test_descr << "Checking the solutions of a constant(="
			<< value << ") function with zero derivatives specified";
		no_crossing_locations[two_value+2].resize(1, NaN);
		check_solutions(interpf, crossing_checks, no_crossing_locations,
				test_descr.str(), 1e-8, 1e-8);
		no_crossing_locations[two_value+2].resize(0);
	}
}

void test_InterpolatingFunctionALGLIB::test_linear_deriv()
{
	for(int two_slope=-2; two_slope<=2; two_slope++) {
		double slope=0.5*two_slope;
		for(int two_zeropt=-2; two_zeropt<=2; two_zeropt++) {
			double zeropt=0.5*two_zeropt;
			std::valarray<double> x(5), y(x.size()), x_check(10), 
				y_check(x_check.size()), crossing_checks(5);
			std::valarray< std::valarray<double> > 
				crossing_locations(5);
			for(int two_y=-2; two_y<=2; two_y++) {
				int ind=two_y+2;
				double y=0.5*two_y;
				crossing_checks[ind]=y;
				double cross=(y-zeropt)/slope;
				if((cross>=0 && cross<=1) || std::isnan(cross))
					crossing_locations[ind].resize(1,cross);
			}
			for(size_t i=0; i<x.size(); i++) {
				x[i]=static_cast<double>(i)/(x.size()-1);
				y[i]=zeropt+slope*x[i];
			}
			for(size_t i=0; i<x_check.size(); i++) {
				x_check[i]=static_cast<double>(i)/(x_check.size()-1);
				y_check[i]=zeropt+slope*x_check[i];
			}
			std::ostringstream test_descr;
			test_descr << "Checking the interpolation of y="
				<< slope << "*x + " << zeropt 
				<< " with derivatives specified";
			InterpolatingFunctionALGLIB interpf(x, y,
						std::valarray<double>(slope, x.size()));
			check_interp(interpf, x_check, y_check,
					test_descr.str());
			test_descr.str("");
			test_descr << "Checking the derivatives of y="
				<< slope << "*x + " << zeropt 
				<< " with derivatives specified";
			check_deriv(interpf, x_check, y_check,
					std::valarray<double>(slope, x_check.size()),
					std::valarray<double>(0.0, x_check.size()),
					test_descr.str());
			test_descr.str("");
			test_descr << "Checking the solutions of y="
				<< slope << "*x + " << zeropt 
				<< " with derivatives specified";
			check_solutions(interpf, crossing_checks, crossing_locations,
					test_descr.str(), 1e-8, 1e-8);
		}
	}
}

void test_InterpolatingFunctionALGLIB::test_quadratic_deriv()
{
	for(int two_quadratic=-2; two_quadratic<=2; two_quadratic++) {
		double quadratic=0.5*two_quadratic;
		for(int two_linear=-2; two_linear<=2; two_linear++) {
			double linear=0.5*two_linear;
			for(int two_zeropt=-2; two_zeropt<=2; two_zeropt++) {
				double zeropt=0.5*two_zeropt;
				std::valarray<double> x(5), y(x.size()), yprime(x.size()), 
					x_check(10), y_check(x_check.size()), 
					yprime_check(x_check.size()), crossing_checks(5);
				std::valarray< std::valarray<double> > 
					crossing_locations(5);
				for(int two_y=-2; two_y<=2; two_y++) {
					int ind=two_y+2;
					double y=0.5*two_y;
					crossing_checks[ind]=y;
					double det=linear*linear-4.0*quadratic*(zeropt-y);
					if(det==0 || quadratic==0) { 
						double sol=(quadratic==0 
								? (y-zeropt)/linear
								: -linear/(2.0*quadratic));
						if((sol>-1e-8 && sol<1+1e-8) || std::isnan(sol))
							crossing_locations[ind].resize(1, sol);
					} else if(det>0) {
						double sol1=(-linear-std::sqrt(det))/(2.0*quadratic),
							   sol2=(-linear+std::sqrt(det))/(2.0*quadratic);
						if(quadratic<0) {
							double buff=sol1;
							sol1=sol2;
							sol2=buff;
						}
						if(sol1>-1e-8 && sol1<1+1e-8) {
							if(sol2>-1e-8 && sol2<1+1e-8) {
								crossing_locations[ind].resize(2);
								crossing_locations[ind][1]=sol2;
							}
							else crossing_locations[ind].resize(1);
							crossing_locations[ind][0]=sol1;
						} else if(sol2>-1e-8 && sol2<1+1e-8)
							crossing_locations[ind].resize(1, sol2);
					}
				}
				for(size_t i=0; i<x.size(); i++) {
					x[i]=static_cast<double>(i)/(x.size()-1);
					y[i]=zeropt+linear*x[i]+quadratic*x[i]*x[i];
					yprime[i]=linear+2.0*quadratic*x[i];
				}
				for(size_t i=0; i<x_check.size(); i++) {
					x_check[i]=static_cast<double>(i)/(x_check.size()-1);
					y_check[i]=zeropt+linear*x_check[i]+
						quadratic*x_check[i]*x_check[i];
					yprime_check[i]=linear+2.0*quadratic*x_check[i];
				}
				std::ostringstream test_descr;
				test_descr << "Checking the interpolation of y="
					<< quadratic << "*x^2 + " << linear << "*x + " 
					<< zeropt << " with derivatives specified";
				InterpolatingFunctionALGLIB interpf(x, y, yprime);
				check_interp(interpf, x_check, y_check, test_descr.str());
				test_descr.str("");
				test_descr << "Checking the derivatives of y="
					<< quadratic << "*x^2 + " << linear << "*x + " 
					<< zeropt << " with derivatives specified";
				check_deriv(interpf, x_check, y_check, yprime_check,
						std::valarray<double>(2.0*quadratic, x_check.size()),
						test_descr.str());
				test_descr.str("");
				test_descr << "Checking the solutions of y="
					<< quadratic << "*x^2 + " << linear << "*x + " 
					<< zeropt << " with derivatives specified";
				check_solutions(interpf, crossing_checks, crossing_locations,
						test_descr.str(), 1e-8, 1e-8);
			}
		}
	}
}

void test_InterpolatingFunctionALGLIB::test_cubic_deriv()
{
	for(int two_cubic=-2; two_cubic<=2; two_cubic++) {
		double cubic=0.5*two_cubic;
		for(int two_quadratic=-2; two_quadratic<=2; two_quadratic++) {
			double quadratic=0.5*two_quadratic;
			for(int two_linear=-2; two_linear<=2; two_linear++) {
				double linear=0.5*two_linear;
				for(int two_zeropt=-2; two_zeropt<=2; two_zeropt++) {
					double zeropt=0.5*two_zeropt;
					std::valarray<double> x(5), y(x.size()), 
						yprime(x.size()), 
						x_check(10), y_check(x_check.size()), 
						yprime_check(x_check.size()),
						ysecond_check(x_check.size());
					for(size_t i=0; i<x.size(); i++) {
						x[i]=static_cast<double>(i)/x.size();
						y[i]=zeropt + linear*x[i] + quadratic*x[i]*x[i] +
							cubic*x[i]*x[i]*x[i];
						yprime[i]=linear + 2.0*quadratic*x[i] + 
							3.0*cubic*x[i]*x[i];
					}
					for(size_t i=0; i<x_check.size(); i++) {
						x_check[i]=static_cast<double>(i)/x_check.size();
						y_check[i]=zeropt + linear*x_check[i] +
							quadratic*x_check[i]*x_check[i] + 
							cubic*x_check[i]*x_check[i]*x_check[i];
						yprime_check[i]=linear + 2.0*quadratic*x_check[i] +
							3.0*cubic*x_check[i]*x_check[i];
						ysecond_check[i]=2.0*quadratic+6.0*cubic*x_check[i];
					}
					std::ostringstream test_descr;
					test_descr << "Checking the interpolation of y="
						<< cubic << "*x^3 + " << quadratic << "*x^2 + " 
						<< linear << "*x + " << zeropt 
						<< " with derivatives specified";	
					InterpolatingFunctionALGLIB interpf(x, y, yprime);
					check_interp(interpf, x_check, y_check,test_descr.str());
					test_descr.str("");
					test_descr << "Checking the derivatives of y="
						<< cubic << "*x^3 + " << quadratic << "*x^2 + " 
						<< linear << "*x + " << zeropt 
						<< " with derivatives specified";
					check_deriv(interpf, x_check, y_check, yprime_check,
							ysecond_check, test_descr.str());
				}
			}
		}
	}
}

void test_InterpolatingFunctionALGLIB::test_sin_smoothed()
{
	unsigned frequency=3;
	const double noise=0.1;
	std::valarray<double> x(1000), y(x.size()),
		x_check(500), y_check(x_check.size()), yprime_check(x_check.size()),
		ysecond_check(x_check.size()), crossing_check(3);
	std::valarray< std::valarray<double> > 
		crossing_locations(std::valarray<double>(2*frequency),
				crossing_check.size()); 
	for(size_t i=0; i<x.size(); i++) {
		x[i]=static_cast<double>(i)/(x.size()-1);
		y[i]=std::sin(2.0*M_PI*frequency*x[i])+
			noise*(static_cast<double>(rand())/RAND_MAX-0.5);
	}
	for(size_t i=0; i<x_check.size(); i++) {
		x_check[i]=static_cast<double>(rand())/RAND_MAX;
		y_check[i]=std::sin(2.0*M_PI*frequency*x_check[i]);
		yprime_check[i]=2.0*M_PI*frequency*std::cos(2.0*M_PI*frequency*x_check[i]);
		ysecond_check[i]=-std::pow(2.0*M_PI*frequency, 2)*y_check[i];
	}
	for(size_t i=0; i<crossing_check.size(); i++) {
		double y=
			2.0*(static_cast<double>(i+1)/(crossing_check.size()+1)-0.5),
			   loc_base=std::asin(y);
		crossing_check[i]=y;
		size_t loc_ind=0;
		for(size_t n=0; n<frequency; n++)
			if(loc_base>=0) {
				crossing_locations[i][loc_ind++]=2.0*M_PI*n+loc_base;
				crossing_locations[i][loc_ind++]=(2.0*n+1.0)*M_PI-loc_base;
			} else {
				crossing_locations[i][loc_ind++]=(2.0*n+1.0)*M_PI-loc_base;
				crossing_locations[i][loc_ind++]=2.0*M_PI*(n+1)+loc_base;
			}
	}
	std::ostringstream test_descr;
	test_descr << "Checking the smoothed interpolation of y=sin(" 
		<< 2*frequency << "*pi*x) + " << noise << "*U(0,1)";
	InterpolatingFunctionALGLIB interpf(x, y, std::valarray<double>(), 6.5);
	check_interp(interpf, x_check, y_check, test_descr.str(), 0.1*noise, 
			0.1*noise);
	test_descr.str("");
	test_descr << "Checking the smoothed derivative of y=sin(" 
		<< 2*frequency << "*pi*x) + " << noise << "*U(0,1)";
	check_deriv(interpf, x_check, y_check, yprime_check, 
			std::valarray<double>(), test_descr.str(), 0.3*noise, 1.8*M_PI*noise);
	check_deriv(interpf, x_check, y_check, yprime_check, 
			ysecond_check, test_descr.str(), 2.0*noise, 72.0*M_PI*M_PI*noise);
}

test_InterpolatingFunctionALGLIB::test_InterpolatingFunctionALGLIB()
{
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_const_no_deriv)
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_linear_no_deriv)
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_const_deriv)
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_linear_deriv)
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_quadratic_deriv)
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_cubic_deriv)
	TEST_ADD(test_InterpolatingFunctionALGLIB::test_sin_smoothed)
}

#ifdef STANDALONE
int main()
{
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_InterpolatingFunctionALGLIB tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);  
}
#endif
