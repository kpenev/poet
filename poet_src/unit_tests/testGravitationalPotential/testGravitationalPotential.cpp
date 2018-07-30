#include "testGravitationalPotential.h"

int main(int argc, char **argv)
{
    std::cout.precision(16);
    std::cerr.precision(16);
    testGravitationalPotential::EccentricOrbit orbit(1.0, 0.1, M_PI, 0.8);
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
