/**\file
 *
 * \brief Declarses a class for functions that are the inverse of some
 * analytical function.
 *
 * Inversion is performed numerically.
 *
 * \ingroup UnitTests_group
 */

#ifndef __INVERSE_FUNCTION_H
#define __INVERSE_FUNCTION_H

#include "../../Core/Functions.h"
#include "../../Core/InterpSolutionIterator.h"
#include <gsl/gsl_roots.h>

///\brief The invrse of an existing function.
///
///Uses the last returned result as a guess for the next solution, so this class
///workes best when invoked with consecutive arguments close to each other, like
///for example evaluating it at the set of ages of a POET solution in
///ascending/descending order.
///
///In other cases, use set_guess before each query.
class InverseFunction : public Core::OneArgumentDiffFunction {
private:
    ///The function being inverted.
    const Core::OneArgumentDiffFunction &__to_invert;

    ///\brief The relative tolerance in the value returned by the function at
    ///the best guess for the solution.
    double __tolerance;

    ///The GSL derivative-based solver.
    gsl_root_fdfsolver *__solver;

    ///The fdf argument used by the GSL solver.
    gsl_function_fdf __solver_fdf;

    ///The guess to use for the next solution.
    mutable double __guess;

    ///The value we are trying to match the function to.
    mutable double __target;

    ///GLS format function to invert.
    friend double gsl_f(double x, void *params);

    ///GLS format derivative of function to invert.
    friend double gsl_df(double x, void *params);

    ///GLS format function and its derivative to invert.
    friend void gsl_fdf(double x, void *params, double *f, double *df);

public:
    ///Invert the given function.
    InverseFunction(const OneArgumentDiffFunction &to_invert,
                    double initial_guess,
                    double tolerance = 1e-10);

    ///The value of the function at the given abscissa.
    double operator()(double x) const;

    ///The lower end of the range over which the function is defined is uknown.
    double range_high() const
    {
        throw Core::Error::Runtime(
            "Upper end of inverse function range is unknown."
        );
    }

    ///The upper end of the range over which the function is defined
    double range_low() const
    {
        throw Core::Error::Runtime(
            "Lower end of inverse function range is unknown."
        );
    }

    ///\brief An iterator over the abscissas where the function takes
    ///the given y value.
    Core::InterpSolutionIterator crossings(double = 0) const
    {
        throw Core::Error::Runtime(
            "Finding all solutinos of an inverse function not implemented."
        );
    };

    ///\brief Returns a pointer to the derivative of the function.
    ///
    ///The use of a pointer allows avoiding potentially expensive copy
    ///opertaions.
    const Core::FunctionDerivatives *deriv(double x) const;

    ~InverseFunction();
};

#endif
