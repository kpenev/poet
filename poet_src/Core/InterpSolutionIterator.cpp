#include "InterpSolutionIterator.h"

namespace Core {
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

    InterpSolutionIterator::InterpSolutionIterator(
        const InterpSolutionIterator &rhs
    ) : 
        spline(rhs.spline),
        node_index(rhs.node_index), y(rhs.y),
        solution_iter(rhs.solution_iter),
        solutions(rhs.solutions), 
        is_out_of_range(false)
    {}

    InterpSolutionIterator::InterpSolutionIterator(
        const SerializableSpline1dInterpolant &spline_var,
        double offset,
        double min_sol_distance
    ) : 
        spline(spline_var.c_ptr()),
        node_index(0),
        min_diff(min_sol_distance),
        y(offset),
        is_out_of_range(false)
    {
        get_solutions();
        solution_iter=solutions.begin();
    }

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

    InterpSolutionIterator InterpSolutionIterator::operator++(int) 
    {
        InterpSolutionIterator result(*this);
        operator++();
        return result;
    }

    InterpSolutionIterator &InterpSolutionIterator::operator--()
    {
        if(solution_iter!=solutions.begin()) solution_iter--;
        else is_out_of_range=true;
        return *this;
    }

    InterpSolutionIterator InterpSolutionIterator::operator--(int)
    {
        InterpSolutionIterator result(*this);
        operator--();
        return result;
    }

    const double &InterpSolutionIterator::operator*() const
    {
        return *solution_iter;
    }

    bool InterpSolutionIterator::operator==(
        const InterpSolutionIterator &rhs
    ) const
    {
        return (spline==rhs.spline && node_index==rhs.node_index && 
                solution_iter==rhs.solution_iter);
    }

    bool InterpSolutionIterator::operator!=(
        const InterpSolutionIterator &rhs
    ) const
    {
        return !operator==(rhs);
    }

    bool InterpSolutionIterator::out_of_range() const
    {
        return is_out_of_range;
    }

}
