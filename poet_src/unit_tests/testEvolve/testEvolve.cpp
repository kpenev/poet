#include "testOrbitSolver.h"

#ifdef STANDALONE
int main()
{
    Evolve::DissipatingZone::read_eccentricity_expansion(
        "eccentricity_expansion_coef.txt"
    );

	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	std::cerr.setf(std::ios_base::scientific);
	std::cerr.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);
    Evolve::test_OrbitSolver tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);
    return 0;
}
#endif
