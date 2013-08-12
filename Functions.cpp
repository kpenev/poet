#include "Functions.h"
#include "Common.h"
#include <iostream>

#ifndef NO_SERIALIZE
BOOST_CLASS_EXPORT_IMPLEMENT(InterpolatingFunctionALGLIB)
#endif

const int max_default_degrees_of_freedom=1000;

///Adds the solutions at the present node (or the closest node after
///it if there are none) to the end of solutions, incrementing the 
///node to one past the node from which solutions were added.
void InterpSolutionIterator::get_solutions()
{
	bool found=false;
	while(node_index < spline->n - 1 && !found) {
		double *coef=(spline->c.ptr.p_double + 4*node_index);
		double x0=spline->x.ptr.p_double[node_index],
		       max_x=spline->x.ptr.p_double[node_index+1]-x0;
		std::valarray<double> cubic_sol=solve_cubic(
			coef[0]-y, coef[1], coef[2], coef[3]);
		double tolerance=min_diff*max_x;
		for(size_t sol_i=0; sol_i<cubic_sol.size(); sol_i++) { 
			double sol_value=cubic_sol[sol_i];
			if(sol_value>-tolerance && sol_value<max_x+tolerance
					&& (solutions.size()==0 || 
						solutions.back()<sol_value+x0-tolerance)) {
				solutions.push_back(sol_value+x0);
				found=true;
			}
		}
		node_index++;
	}
	is_out_of_range=!found;
}

///Copy constructor
InterpSolutionIterator::InterpSolutionIterator(
		const InterpSolutionIterator &rhs)
	: spline(rhs.spline), node_index(rhs.node_index), y(rhs.y),
	solution_iter(rhs.solution_iter), solutions(rhs.solutions), 
	is_out_of_range(false)
{}

///Start iterating over the solution of the given spline.
InterpSolutionIterator::InterpSolutionIterator(
		const alglib::spline1dinterpolant &spline_var, double offset,
		double min_sol_distance)
	: spline(spline_var.c_ptr()), node_index(0), min_diff(min_sol_distance),
	y(offset), is_out_of_range(false)
{
	get_solutions();
	solution_iter=solutions.begin();
}

///Go to the next solution.
InterpSolutionIterator &InterpSolutionIterator::operator++() 
{
	solution_iter++;
	if(solution_iter==solutions.end()) {
		solution_iter--;
		get_solutions();
		solution_iter++;
	}
	is_out_of_range=(solution_iter==solutions.end());
	return *this;
}

///Go to the next solution.
InterpSolutionIterator InterpSolutionIterator::operator++(int) 
{
	InterpSolutionIterator result(*this);
	operator++();
	return result;
}

///Go to the previous solution.
InterpSolutionIterator &InterpSolutionIterator::operator--()
{
	if(solution_iter!=solutions.begin()) solution_iter--;
	else is_out_of_range=true;
	return *this;
}

///Go to the previous solution.
InterpSolutionIterator InterpSolutionIterator::operator--(int)
{
	InterpSolutionIterator result(*this);
	operator--();
	return result;
}

///Returns the current solution
const double &InterpSolutionIterator::operator*() const
{
	return *solution_iter;
}

///Checks if this iterator is at the same solution of the same
///spline.
bool InterpSolutionIterator::operator==(const InterpSolutionIterator &rhs) 
	const
{
	return (spline==rhs.spline && node_index==rhs.node_index && 
			solution_iter==rhs.solution_iter);
}

///The opposite of operator==.
bool InterpSolutionIterator::operator!=(const InterpSolutionIterator &rhs) 
	const
{
	return !operator==(rhs);
}

///Whether the currenty iterator is either before the first 
///solution or past the last one.
bool InterpSolutionIterator::out_of_range() const
{
	return is_out_of_range;
}

///Creates the derivative variable with the values of the function 
///and first and second derivatives as specified. All higher order 
///derivatives are zero.
CubicSplineDerivatives::CubicSplineDerivatives(double func_value, 
		double first_deriv, double second_deriv) :
	zeroth(func_value), first(first_deriv), second(second_deriv)
{}

double CubicSplineDerivatives::order(unsigned deriv_order) const
{
	if(deriv_order==0) return zeroth;
	else if(deriv_order==1) return first;
	else if(deriv_order==2) return second;
	else if(deriv_order==3) return NaN;
	else return 0.0;
}

///Constuct an interpolating function based on the ALGLIB (smoothing)
///cubic spline interpolation.
InterpolatingFunctionALGLIB::InterpolatingFunctionALGLIB(
		///The abscissas of tabulated points from the curve to fit.
		const std::valarray<double> &x,
		///The ordinates of tabulated points from the curve to fit.
		const std::valarray<double> &y, 
		///The values of the derivatives to impose on the nodes.
		const std::valarray<double> &yprime,
		///How much smoothing to apply (omit for no smoothnig, i.e.
		///interpolating curve passes through all the points).
		double smoothing,
		///How many degrees of freedom to use, ignored for non-smoothing
		///interpolation, if omitted it is set to 3 times the number of
		///points being fitted.
		int degrees_of_freedom)
{
	if(degrees_of_freedom<0) degrees_of_freedom=std::min(
			max_default_degrees_of_freedom, static_cast<int>(3*x.size()));

	alglib::real_1d_array alglib_x, alglib_y, alglib_yprime;
	alglib_x.setcontent(x.size(), &x[0]);
	alglib_y.setcontent(x.size(), &y[0]);
	if(yprime.size()) alglib_yprime.setcontent(yprime.size(), &yprime[0]);
	min_x=x.min();
	max_x=x.max();

	if(std::isnan(smoothing)) {
		if(yprime.size())
			alglib::spline1dbuildhermite(alglib_x, alglib_y,
					alglib_yprime, spline);
		else 
			alglib::spline1dbuildcubic(alglib_x, alglib_y, spline);
	}
	else {
		if(yprime.size())
			throw Error::BadFunctionArguments(
				"Smoothing not supported when derivatives "
				"are specified in InterpolatingFunctionALGLIB "
				"constructor.");
		alglib::ae_int_t fit_info;
		alglib::spline1dfitreport fit_report;
		alglib::spline1dfitpenalized(alglib_x, alglib_y, 
				degrees_of_freedom, smoothing, fit_info, 
				spline, fit_report);
		if(fit_info<=0) throw Error::ALGLIB("Spline fit failed.");
	}
}

///Returns an iterator over the abscissas where the function takes 
///the given y value.
InterpSolutionIterator InterpolatingFunctionALGLIB::crossings(double y) 
	const
{
	return InterpSolutionIterator(spline, y);
}
