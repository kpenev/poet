/**\file
 *
 * \brief A hierarchy of classes representing functions.
 * 
 * \ingroup Utilities_group
 */

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


///\brief An iterator over a set of solutions to an interpolating function.
///
///Iterates over all abscissas where the intepolating function takes a
///pre-determined value, merging solutions that are too close together.
///
///Works by adding solutions between later and later pairs of consecutive
///interpolation nodes whenever it runs out of found solutions and the next
///solution is requested.
class InterpSolutionIterator : 
	public std::iterator<std::input_iterator_tag, double, ptrdiff_t, 
		const double*, const double&> {
private:
	///The ALGLIB spline.
	const alglib_impl::spline1dinterpolant *spline;

	///The node up to which solutions have been reported
	int node_index;

	///\brief Controls when solutions are considered distinct.
	///
	///Solutions that differ by less than min_diff fraction of the
	///distance between two nodes are treated as a single solution.
	double min_diff,
	       
	       ///Iterate over abscissas when the interpolation=y.
	       y;

	///An iterator over the list of solutions found so far. 
	std::list<double>::const_iterator solution_iter;

	///The list of solutions found so far.
	std::list<double> solutions;

	///Whether we have gone past the last solution or before the first.
	bool is_out_of_range;

	///\brief Find another set of solutions.
	///
	///Adds the solutions between the the present node and the next (or
	///the closest node after it if there are none) to the end of
	///solutions, incrementing the node to one past the node from which
	///solutions were added.
	void get_solutions();
public:
	///Default constructor of a non meaningful object
	InterpSolutionIterator() {};

	///Copy constructor
	InterpSolutionIterator(const InterpSolutionIterator &rhs);

	///Start iterating over the solutions of the given spline.
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

	///\brief Checks if this iterator is at the same solution of the same
	///spline as rhs.
	bool operator==(const InterpSolutionIterator &rhs) const;

	///\brief The opposite of operator==.
	bool operator!=(const InterpSolutionIterator &rhs) const;

	///Whether we have gone past the last solution or before the first.
	bool out_of_range() const; 
};

///\brief The base class for functions which take a single argument and
///return a single value.
template<class InputDataType, class OutputDataType>
class OneArgumentFunction {
#ifndef NO_SERIALIZE
	///Needed for serialization to work.
	friend class boost::serialization::access;

	///Serialize this function.
	template<class Archive>
	void serialize(Archive &, const unsigned int) {}
#endif
public:
	///The value of the function at the given abscissa.
	virtual OutputDataType operator()(InputDataType in_value) const =0;

	///The lower end of the range over which the function is defined
	virtual InputDataType range_high() const=0;

	///The upper end of the range over which the function is defined
	virtual InputDataType range_low() const=0;

	///\brief An iterator over the abscissas where the function takes the
	///given y value.
	virtual InterpSolutionIterator crossings(double y=0) const =0;

	///Provide a virtual destructor for a virtual class.
	virtual ~OneArgumentFunction() {};
};

///A class representing arbitrary order derivatives of a function.
class FunctionDerivatives {
public:
	///\brief Derivative of the given order of the function with respect
	///to its argument.
	virtual double order(unsigned deriv_order=1) const =0;

	///Clean up.
	virtual ~FunctionDerivatives() {};
};

///A class for the derivatives of a cubic spline (=0 for order>2).
class CubicSplineDerivatives : public FunctionDerivatives {
private:
	double zeroth, ///< The value of the function
	       first, ///< The first derivative
	       second; ///< The second derivative
public:
	///\brief Constuct a spline derivative.
	///
	///Creates the derivative variable with the values of the function 
	///and first and second derivatives as specified. All higher order 
	///derivatives are zero.
	CubicSplineDerivatives(
			///The value of the function (zeroth derivative)
			double func_value,
			
			///The first derivative.
			double first_deriv, 

			///The second derivative.
			double second_deriv);

	///Returns the derivative of the given order (zero is allowed).
	double order(unsigned deriv_order=1) const;
};

///A class representing a once differentiable function of a single argument
class OneArgumentDiffFunction : public OneArgumentFunction<double,double> {
private :
#ifndef NO_SERIALIZE
	///Needed for serialization to work.
	friend class boost::serialization::access;

	///Serialize this function.
	//The second parameter is supposed to be version
	template<class Archive>
	void serialize(Archive & ar, const unsigned int) {
		ar & boost::serialization::base_object< OneArgumentFunction<double,double> >(*this);
	}
#endif
public:
	///\brief Returns a pointer to the derivative of the function.
	///
	///The use of a pointer allows avoiding potentially expensive copy
	///opertaions.
	virtual const FunctionDerivatives *deriv(double x) const=0;

	virtual ~OneArgumentDiffFunction() {}
};

///Function which interpolates, with possible smoothing, between points.
class InterpolatingFunctionALGLIB : public OneArgumentDiffFunction {
private:
#ifndef NO_SERIALIZE
	///Needed for serialization to work.
	friend class boost::serialization::access;

	///Serialize this function.
	//The second parameter should eventually be version
	template<class Archive>
	void serialize(Archive & ar, const unsigned int) {
		ar & boost::serialization::base_object<OneArgumentDiffFunction>(
				*this);
		ar & spline;
		ar & min_x & max_x;
	}
#endif

	///\brief The interpolating function information necessary to
	///evaluate at any given point.
	alglib::spline1dinterpolant spline;

	///The smallest abscissa covered by the spline points.
	double min_x,

	       ///The largest abscissa covered by the spline points.
	       max_x;
public:
	///Needed by the Boost serializer.
	InterpolatingFunctionALGLIB() {};

	///\brief Constuct an interpolating function.
	///
	///Based on the ALGLIB (smoothing) cubic spline interpolation.
	InterpolatingFunctionALGLIB(
			///The abscissas of tabulated points to fit.
			const std::valarray<double> &x,

			///The ordinates of tabulated points to fit.
			const std::valarray<double> &y, 

			///The values of the derivatives to impose on the
			///nodes.
			const std::valarray<double> &yprime=
				std::valarray<double>(),

			///How much smoothing to apply. Omit for no
			///smoothnig, i.e. the interpolating curve passes
			///through all the points).
			double smoothing=NaN,
			
			///How many degrees of freedom to use for smoothing
			///interpolation. Ignored for non-smoothing
			///interpolation. If omitted it is set to 3 times the
			///number of points being fitted.
			int degrees_of_freedom=-1);

	///Returns the value of the interpolating function at the given abscissa.
	inline double operator()(double x) const
	{return alglib::spline1dcalc(spline, x);}

	///\brief The derivatives of the interpolating function at the given
	///abscissa.
	///
	///Returns a newly allocated structure, which must be destroyed
	///when no longer needed.
	inline const CubicSplineDerivatives *deriv(double x) const
	{
		double v, dv, d2v; alglib::spline1ddiff(spline, x, v, dv, d2v); 
		return new CubicSplineDerivatives(v, dv, d2v);
	}

	///The maximum abscissa at which the function is defined.
	inline double range_high() const {return max_x;}

	///The minimum abscissa at which the function is defined.
	inline double range_low() const {return min_x;}

	///\brief Iterator over the abscissas where the function takes the
	///given y value.
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
