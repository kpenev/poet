#include "testGravitationalPotential.h"

namespace testGravitationalPotential {

    void print_orbit()
    {
        EccentricOrbit orbit(1.0, 0.1, M_PI, 0.8);
        std::cout << std::setw(25) << "phase"
                  << std::setw(25) << "x"
                  << std::setw(25) << "y"
                  << std::setw(25) << "z"
                  << std::endl;
        for(double phase = 0.0; phase < 4.0 * M_PI; phase += 0.001 * M_PI) {
            Eigen::Vector3d secondary_position = orbit.secondary_position(phase);
            std::cout << std::setw(25) << phase 
                      << std::setw(25) << secondary_position[0]
                      << std::setw(25) << secondary_position[1]
                      << std::setw(25) << secondary_position[2]
                      << std::endl;
        }
    }

    void print_tidal_potential()
    {
        double mprimary=1.0,
               msecondary=0.1,
               semimajor=M_PI,
               eccentricity=0.0,
               inclination=0.0,
               periapsis=0.0;
        TidalPotential exact_potential(mprimary,
                                       msecondary,
                                       semimajor,
                                       eccentricity,
                                       inclination,
                                       periapsis);
        TidalPotentialExpansion approx_potential(mprimary,
                                                 msecondary,
                                                 semimajor,
                                                 eccentricity,
                                                 inclination,
                                                 periapsis);
        double orbital_period = EccentricOrbit(mprimary,
                                               msecondary,
                                               semimajor,
                                               eccentricity).orbital_period();
        Eigen::Vector3d position(1.0e-3, 0.0, 0.0);
        std::ostringstream Uexact_label, Uapprox_label;
        Uexact_label << "Uexact("
                     << "x=" << position[0]
                     << ", y=" << position[1]
                     << ", z=" << position[2]
                     << ", t=time)";
        Uapprox_label << "Uapprox("
                      << "x=" << position[0]
                      << ", y=" << position[1]
                      << ", z=" << position[2]
                      << ", t=time)";
        std::cout << std::setw(35) << "time/Porb"
                  << std::setw(35) << Uexact_label.str()
                  << std::setw(35) << Uapprox_label.str()
                  << std::endl;
        for(
            double time = 0;
            time < 5.0 * orbital_period;
            time += 0.01 * orbital_period
        ) {
            std::cout << std::setw(35) << time/orbital_period
                      << std::setw(35) << exact_potential(position, time)
                      << std::setw(35) << approx_potential(position, time)
                      << std::endl;
        }
    }
}

int main(int argc, char **argv)
{
    Evolve::TidalPotentialTerms::read_eccentricity_expansion(
        "eccentricity_expansion_coef.txt"
    );
    std::cout.precision(16);
    std::cerr.precision(16);
    testGravitationalPotential::print_tidal_potential();
}
