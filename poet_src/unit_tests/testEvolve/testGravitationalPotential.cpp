/**\file
 *
 * \brief Implement the non-inline methods of test_GravitationalPotential.
 *
 * \ingroup UnitTests_group
 */

#include "testGravitationalPotential.h"

namespace Evolve {

    void test_GravitationalPotential::print_orbit(
        const EccentricOrbit &orbit
    ) const
    {
        std::cout << std::setw(25) << "phase"
                  << std::setw(25) << "x"
                  << std::setw(25) << "y"
                  << std::setw(25) << "z"
                  << std::endl;
        for(double phase = 0.0; phase < 4.0 * M_PI; phase += 0.001 * M_PI) {
            Eigen::Matrix<long double, 3, 1> secondary_position =
                orbit.secondary_position(phase);
            std::cout << std::setw(25) << phase
                      << std::setw(25) << secondary_position[0]
                      << std::setw(25) << secondary_position[1]
                      << std::setw(25) << secondary_position[2]
                      << std::endl;
        }
    }

    void test_GravitationalPotential::print_tidal_potential(
        double primary_mass,
        double secondary_mass,
        double semimajor,
        double eccentricity,
        double inclination,
        double arg_of_periapsis,
        const Eigen::Matrix<long double, 3, 1> &position,
        unsigned expansion_order
    ) const
    {
        TidalPotential exact_potential(primary_mass,
                                       secondary_mass,
                                       semimajor,
                                       eccentricity,
                                       inclination,
                                       arg_of_periapsis);
        TidalPotentialExpansion approx_potential(primary_mass,
                                                 secondary_mass,
                                                 semimajor,
                                                 eccentricity,
                                                 inclination,
                                                 arg_of_periapsis);

        double orbital_period = exact_potential.orbit().orbital_period();
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
                      << std::setw(35) << approx_potential(position,
                                                           time,
                                                           expansion_order)
                      << std::endl;
        }
    }


    double test_GravitationalPotential::abs_precision(
        const Eigen::Matrix<long double, 3, 1> &position,
        const EccentricOrbit &orbit,
        double expansion_precision
    ) const
    {
        const double safety = 3.0;

        double scaled_distance = (
            position.norm()
            /
            (orbit.semimajor() * (1.0 - orbit.eccentricity()))
        );

        return (
            (
                Core::AstroConst::G
                *
                orbit.secondary_mass() * Core::AstroConst::solar_mass
                /
                (
                    orbit.semimajor() * Core::AstroConst::solar_radius
                    *
                    (1.0 - orbit.eccentricity())
                )
            )
            *
            std::pow(scaled_distance, 2)
            *
            (safety * scaled_distance + expansion_precision)
        );
    }

    void test_GravitationalPotential::test_single_point(
        TidalPotentialExpansion &approx_potential,
        const TidalPotential &exact_potential,
        const Eigen::Matrix<long double, 3, 1> &position
    )
    {
        double orbital_period = exact_potential.orbit().orbital_period(),
               abs_tolerance = abs_precision(
                   position,
                   exact_potential.orbit(),
                   approx_potential.get_expansion_precision()
               ),
               eccentricity = exact_potential.orbit().eccentricity();

        unsigned expansion_order =
            TidalPotentialTerms::required_expansion_order(eccentricity);

        std::ostringstream message_start;
        message_start.setf(std::ios_base::scientific);
        message_start.precision(16);
        message_start << "M = " << exact_potential.orbit().primary_mass()
                      << "; M' = " << exact_potential.orbit().secondary_mass()
                      << "; a = " << exact_potential.orbit().semimajor()
                      << "; e = " << eccentricity
                      << "; inclination = " << exact_potential.inclination()
                      << "; periapsis = " << exact_potential.arg_of_periapsis()
                      << "; U(x = " << position[0]
                      << ", y = " << position[1]
                      << ", z = " << position[2];

        for(
            double time = 0;
            time < 5.0 * orbital_period;
            time += 0.03 * M_PI * orbital_period
        ) {
            double expected = exact_potential(position, time);
            double got = approx_potential(position, time, expansion_order);

            std::ostringstream message;
            message.precision(16);
            message.setf(std::ios_base::scientific);
            message << message_start.str()
                    << ", t = " << time
                    << " = " << time / orbital_period
                    << " Porb): expected " << expected
                    << "; got " << got
                    << "; difference " << got - expected
                    << " > " << abs_tolerance;

#ifdef VERBOSE_DEBUG
            std::cerr << message.str() << std::endl;
#endif

            TEST_ASSERT_MSG(
                check_diff(got,
                           expected,
                           0.0,
                           abs_tolerance),
                message.str().c_str()
            );
        }
    }

    void test_GravitationalPotential::test_system(
        double primary_mass,
        double secondary_mass,
        double semimajor,
        double eccentricity,
        double inclination,
        double arg_of_periapsis
    )
    {
        TidalPotential exact_potential(primary_mass,
                                       secondary_mass,
                                       semimajor,
                                       eccentricity,
                                       inclination,
                                       arg_of_periapsis);
        TidalPotentialExpansion approx_potential(primary_mass,
                                                 secondary_mass,
                                                 semimajor,
                                                 eccentricity,
                                                 inclination,
                                                 arg_of_periapsis);

        long double test_offsets[]= {-0.01,
                                     -0.001,
                                     -1e-6,
                                     0.0,
                                     1e-6,
                                     0.001,
                                     0.01};
        unsigned num_offsets = sizeof(test_offsets) / sizeof(long double);
        typedef long double ldbl;
        for(unsigned z_index = 0; z_index < num_offsets; ++z_index)
            for(unsigned y_index = 0; y_index < num_offsets; ++y_index)
                for(unsigned x_index = 0; x_index < num_offsets; ++x_index) {
                    test_single_point(
                        approx_potential,
                        exact_potential,
                        Eigen::Matrix<ldbl, 3, 1> (test_offsets[x_index],
                                                   test_offsets[y_index],
                                                   test_offsets[z_index])
                    );
                }
    }

    void test_GravitationalPotential::test_expansion()
    {
        double test_angles[] = {-M_PI/2, -1.0, -0.1, 0.0, 0.1, 1.0, M_PI/2};
        unsigned num_angles = sizeof(test_angles) / sizeof(double);
        for(
            double e = 0.9;
            e <= 0.95;
            e += 0.1
        ) {
            std::cout << "Eccentricity : " << e << std::endl;

            for(
                unsigned inclination_i = 0;
                inclination_i < num_angles;
                ++inclination_i
            ) {
                for(
                    unsigned periapsis_i = 0;
                    periapsis_i < num_angles;
                    ++periapsis_i
                ) {
                    test_system(1.0,
                                0.1,
                                M_PI,
                                e,
                                test_angles[inclination_i],
                                2.0 * test_angles[periapsis_i]);
                }
            }
        }
    }

    test_GravitationalPotential::test_GravitationalPotential()
    {
        TEST_ADD(test_GravitationalPotential::test_expansion);
    }

} //End Evolve namespace.
