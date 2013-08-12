#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include "Common.h"
#include "Error.h"
//#include <solve_polynomial.h>
#include <cmath>
#include <limits>
#include <valarray>
#include <list>
#include <iterator>
#include <iostream>
#include "alglib/src/interpolation.h"

#ifndef NO_SERIALIZE
	#include <boost/serialization/base_object.hpp>
	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
	#include <boost/serialization/export.hpp>
#endif


///A class representing an iterator over a set of solutions to an 
///interpolating function.
class InterpSolutionIterator : 
	public std::iterator<std::input_iterator_tag, double, ptrdiff_t, 
		const double*, const double&> {
private:
	const alglib_impl::spline1dinterpolant *spline;
	int node_index;
	double min_diff, y;
	std::list<double>::const_iterator solution_iter;
	std::list<double> solutions;
	bool is_out_of_range;

	///Adds the solutions at the present node (or the closest node after
	///it if there are none) to the end of solutions, incrementing the 
	///node to one past the node from which solutions were added.
	void get_solutions();
public:
	///Default constructor of a non meaningful object
	InterpSolutionIterator() {};

	///Copy constructor
	InterpSolutionIterator(const InterpSolutionIterator &rhs);

	///Start iterating over the solution of the given spline.
	InterpSolutionIterator(const alglib::spline1dinterpolant &spline_var,
			double offset, double min_sol_distance=1e-8);

	///Go to the next solution.
	InterpSolutionIterator &operator++();

	///Go to the next solution.
	InterpSolutionIterator operator++(int);

	///Go to the previous solution.
	InterpSolutionIterator &operator--();

	///Go to the previous solution.
	InterpSolutionIterator operator--(int);

	///Returns the current solution
	const double &operator*() const;

	///Checks if this operator is at the same solution of the same
	///spline.
	bool operator==(const InterpSolutionIterator &rhs) const;

	///The opposite of operator==.
	bool operator!=(const InterpSolutionIterator &rhs) const;

	///Whether the currenty iterator is either before the first 
	///solution or past the last one.
	bool out_of_range() const; 
};

///A class representing a function which takes a single argument and returns
///a single value.
template<class InputDataType, class OutputDataType>
class OneArgumentFunction {
#ifndef NO_SERIALIZE
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {}
#endif
public:
	virtual OutputDataType operator()(InputDataType in_value) const =0;
	///The lower end of the range over which the function is defined
	virtual InputDataType range_high() const=0;
	///The upper end of the range over which the function is defined
	virtual InputDataType range_low() const=0;
	///Returns an iterator over the abscissas where the function takes 
	///the given y value.
	virtual InterpSolutionIterator crossings(double y=0) const =0;
};

///A class representing arbitrary order derivatives of a function.
class FunctionDerivatives {
public:
	///Returns the derivative of the given order of the function with 
	///respect to its argument.
	virtual double order(unsigned deriv_order=1) const =0;
	virtual ~FunctionDerivatives() {};
};

///A class for the derivatives of a cubic spline (=0 for order>2).
class CubicSplineDerivatives : public FunctionDerivatives {
private:
	double zeroth, ///< The value of the function
	       first, ///< The first derivative
	       second; ///< The second derivative
public:
	///Creates the derivative variable with the values of the function 
	///and first and second derivatives as specified. All higher order 
	///derivatives are zero.
	CubicSplineDerivatives(double func_value, double first_deriv, 
			double second_deriv);

	double order(unsigned deriv_order=1) const;
};

///A class representing a once differentiable function of a single argument
class OneArgumentDiffFunction : public OneArgumentFunction<double,double> {
private :
#ifndef NO_SERIALIZE
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & boost::serialization::base_object< OneArgumentFunction<double,double> >(*this);
	}
#endif
public:
	///Returns a pointer to the derivative of the function (allows aviding 
	///potentially expensive copy opertaions).
	virtual const FunctionDerivatives *deriv(double x) const=0;
};

///A class which returns a function which interpolates, with possible
///smoothing, between a set of points.
class InterpolatingFunctionALGLIB : public OneArgumentDiffFunction {
private:
#ifndef NO_SERIALIZE
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & boost::serialization::base_object<OneArgumentDiffFunction>(
				*this);
		ar & spline;
		ar & min_x & max_x;
	}
#endif

	///The interpolating function information necessary to evaluate at a
	///given point.
	alglib::spline1dinterpolant spline;

	///The abscissa range covered by the spline points
	double min_x, max_x;
public:
	InterpolatingFunctionALGLIB() {}; //needed by Boost serializer
	///Constuct an interpolating function based on the ALGLIB (smoothing)
	///cubic spline interpolation.
	InterpolatingFunctionALGLIB(const std::valarray<double> &x,
			const std::valarray<double> &y, 
			const std::valarray<double> &yprime=std::valarray<double>(),
			double smoothing=NaN, int degrees_of_freedom=-1);

	///Returns the value of the interpolating function at the given abscissa.
	inline double operator()(double x) const
	{return alglib::spline1dcalc(spline, x);}

	///Returns all the derivative of the interpolating function at the 
	///given abscissa in a newly allocated structure (must be destroyed 
	///when no longer needed).
	inline const CubicSplineDerivatives *deriv(double x) const
	{
		double v, dv, d2v; alglib::spline1ddiff(spline, x, v, dv, d2v); 
		return new CubicSplineDerivatives(v, dv, d2v);
	}

	inline double range_high() const {return max_x;}
	inline double range_low() const {return min_x;}

	///Returns an iterator over the abscissas where the function takes 
	///the given y value.
	InterpSolutionIterator crossings(double y=0) const;
};
//BOOST_CLASS_EXPORT_GUID(OneArgumentDiffFunction, "OneArgumentDiffFunction")
//

#ifndef NO_SERIALIZE
BOOST_CLASS_EXPORT_KEY(InterpolatingFunctionALGLIB)
#endif
//BOOST_CLASS_EXPORT_GUID(InterpolatingFunctionALGLIB,
	//	"InterpolatingFunctionALGLIB")
#endif
