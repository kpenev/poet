/**\file
 *
 * \brief Definitions of some of the methods of the InverseFunction class.
 *
 * \ingroup UnitTests_group
 */

#include "InverseFunction.h"

double gsl_f(double x, void *params)
{
    InverseFunction *func_object = reinterpret_cast<InverseFunction*>(params);
    return (func_object->__to_invert(x) - func_object->__target);
}

double gsl_df(double x, void *params)
{
    const Core::FunctionDerivatives *deriv = 
        reinterpret_cast<InverseFunction*>(params)->__to_invert.deriv(x);
    double result = deriv->order(1);
    delete deriv;
    return result;
}

void gsl_fdf(double x, void *params, double *f, double *df)
{
    InverseFunction *func_object = reinterpret_cast<InverseFunction*>(params);
    const Core::FunctionDerivatives *deriv = func_object->__to_invert.deriv(x);
    *f = deriv->order(0) - func_object->__target;
    *df = deriv->order(1);
    delete deriv;
}

InverseFunction::InverseFunction(const OneArgumentDiffFunction &to_invert,
                                 double initial_guess,
                                 double tolerance) :
    __to_invert(to_invert),
    __tolerance(tolerance),
    __solver(gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_steffenson)),
    __guess(initial_guess)
{
    __solver_fdf.f = gsl_f;
    __solver_fdf.df = gsl_df;
    __solver_fdf.fdf = gsl_fdf;
    __solver_fdf.params = reinterpret_cast<void*>(this);

    if(__solver == NULL)
        throw Core::Error::Runtime(
            "Failed to allocate GSL solver in inverse function."
        );

}

double InverseFunction::operator()(double x) const
{
    __target = x;
    if(
        gsl_root_fdfsolver_set(__solver,
                               const_cast<gsl_function_fdf*>(&__solver_fdf),
                               __guess)
    )
        throw Core::Error::Runtime(
            "Failed to initialize solver for inverse function."
        );
    double evaluation_error = Core::Inf, result = Core::NaN;
    while(evaluation_error > __tolerance) {
        if(gsl_root_fdfsolver_iterate(__solver))
            throw Core::Error::Runtime("Error iterating function inversion.");
        result = gsl_root_fdfsolver_root(__solver);
        evaluation_error = std::abs(__to_invert(result) - x);
    }
    __guess = result;
    assert(!std::isnan(result));
    return result;
}

const Core::FunctionDerivatives *InverseFunction::deriv(double x) const
{
    double value = operator()(x);
    const Core::FunctionDerivatives *deriv = __to_invert.deriv(value);
    double first_deriv = 1.0 / deriv->order(1);
    double second_deriv = deriv->order(2) * std::pow(first_deriv, 3);
    delete deriv;
    return new Core::CubicSplineDerivatives(value, first_deriv, second_deriv);
}

InverseFunction::~InverseFunction()
{
    gsl_root_fdfsolver_free(__solver);
}
