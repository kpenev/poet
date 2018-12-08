#include "testGravitationalPotential.h"
#include "testOrbitSolver.h"
#include "testTidalTorquePower.h"
#include "testDifferentialEquations.h"

#ifdef STANDALONE

#include <iostream>
#include <iomanip>

class SinFunction : public Core::OneArgumentDiffFunction {
private:
    double __frequency, __phase;
public:
    SinFunction(double frequency = 1.0, double phase = 0.0) :
        __frequency(frequency), __phase(phase) {}

    double operator()(double x) const
    {return std::sin(__frequency * x + __phase);}

    double range_high() const {return Core::Inf;}
    double range_low() const {return -Core::Inf;}

    ///\brief An iterator over the abscissas where the function takes
    ///the given y value.
    Core::InterpSolutionIterator crossings(double = 0) const
    {
        throw Core::Error::Runtime(
            "Finding all solutinos of Sin is not implemented."
        );
    };

    ///\brief Returns a pointer to the derivative of the function.
    ///
    ///The use of a pointer allows avoiding potentially expensive copy
    ///opertaions.
    const Core::FunctionDerivatives *deriv(double x) const
    {
        return new Core::CubicSplineDerivatives(
            operator()(x),
            __frequency * std::cos(__frequency * x + __phase),
            -std::pow(__frequency, 2) * operator()(x)
        );
    }
};


int main()
{
    std::cout.setf(std::ios_base::scientific);
    std::cout.precision(16);

    /*
    SinFunction sin_function;
    InverseFunction arcsin_function(sin_function, -1.0);

    std::cout << std::setw(25) << "x"
              << std::setw(25) << "f(x)"
              << std::setw(25) << "inverted"
              << std::setw(25) << "arcsin(x)"
              << std::endl;
    for(double x = -1.0; x <= 1.0; x+=1e-2) {
        double solution = arcsin_function(x),
               expected = std::asin(x);
        solution -= 2.0 * M_PI * std::floor(solution / (2.0 * M_PI));
        expected -= 2.0 * M_PI * std::floor(expected / (2.0 * M_PI));
        if(solution > 0.5 * M_PI && solution < 1.5 * M_PI)
            solution = 3.0 * M_PI - solution;
        std::cout << std::setw(25) << x
                  << std::setw(25) << sin_function(solution)
                  << std::setw(25) << solution
                  << std::setw(25) << expected
                  << std::endl;
    }
    return 0;

    Oblique10LinearQuantity q10(1.05 * M_PI, M_PI, 0.1 * M_PI);
    Oblique20LinearQuantity q20(1.05 * M_PI, M_PI, 0.1 * M_PI);
    std::cerr << std::setw(25) << "S"
              << std::setw(25) << "L10(S)"
              << std::setw(25) << "L20(S)"
              << std::endl;
    for(double s = 0.05 * M_PI; s <= 0.1 * M_PI; s+= 0.01 * M_PI)
        std::cerr << std::setw(25) << s
                  << std::setw(25) << q10(s)
                  << std::setw(25) << q20(s)
                  << std::endl;
    return 0;
    */

    Evolve::TidalPotentialTerms::read_eccentricity_expansion(
        "eccentricity_expansion_coef_O200.txt"
    );

	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	std::cerr.setf(std::ios_base::scientific);
	std::cerr.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);

    Test::Suite all_tests;
    all_tests.add(
        std::auto_ptr<Test::Suite>(new Evolve::test_GravitationalPotential)
    );
    all_tests.add(
        std::auto_ptr<Test::Suite>(new Evolve::test_TidalTorquePower)
    );
    all_tests.add(
        std::auto_ptr<Test::Suite>(new Evolve::test_DifferentialEquations)
    );
    all_tests.add(
        std::auto_ptr<Test::Suite>(new Evolve::test_OrbitSolver)
    );
    return (all_tests.run(output)
            ? EXIT_SUCCESS
            : EXIT_FAILURE);
}
#endif
