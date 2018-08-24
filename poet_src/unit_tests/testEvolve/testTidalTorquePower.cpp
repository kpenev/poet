/**\file
 * 
 * \brief Define the non-inline methods of test_TidalTorquePower and actually
 * run the tests if compiled in stand-alone mode.
 *
 * \ingroup UnitTests_group
 */

#include "testTidalTorquePower.h"

namespace Evolve {
    void test_TidalTorquePower::test_aligned_orbit()
    {
        const double abs_tolerance=1e-6;
        std::valarray<double> one_poly(1.0, 1);

        ConstPhaseLagZone zone(1.0, one_poly, one_poly, one_poly);

        zone.configure(true,//initialize
                       1.0,//age
                       1.0,//orbital frequency
                       0.0,//eccentricity
                       1.0,//orbital_angmom
                       0.0,//spin
                       0.0,//inclination
                       0.0,//periapsis
                       true);//spin_is_freq

        Planet::LockedPlanet planet(1.0, 1.0);
        BinarySystem dummy_system(planet, planet, "");

        zone.change_e_order(20, dummy_system, true, 0);

        std::cout << std::setw(25) << "e"
                  << std::setw(25) << "dE_dt"
                  << std::setw(25) << "dLx_dt"
                  << std::setw(25) << "dLy_dt"
                  << std::setw(25) << "dLz_dt"
                  << std::endl;
        for(double e = 0.0; e <= 0.0075; e += 0.0001) {
            zone.configure(false,//initialize
                           1.0,//age
                           1.0,//orbital frequency
                           e,//eccentricity
                           1.0,//orbital_angomm
                           0.0,//spin
                           0.0,//inclination
                           0.0,//periapsis
                           true);//spin_is_freq
            double factor = 1.2 * M_PI,
                   expected_tidal_power = factor * (1.0 + 14.25 * e * e),
                   expected_z_torque = factor * (1.0 + 7.5 * e * e),
                   got_tidal_power = zone.tidal_power(true),
                   got_x_torque = zone.tidal_torque_x(true),
                   got_y_torque = zone.tidal_torque_y(true),
                   got_z_torque = zone.tidal_torque_z(true);
            std::cout << std::setw(25) << e
                      << std::setw(25) << got_tidal_power
                      << std::setw(25) << got_x_torque
                      << std::setw(25) << got_y_torque
                      << std::setw(25) << got_z_torque
                      << std::endl;

            std::ostringstream message;
            message << "For e = " << e
                    << " tidal power = " << got_tidal_power
                    << ", expected: " << expected_tidal_power
                    << "; difference: " << (got_tidal_power
                                            -
                                            expected_tidal_power)
                    << " > " << abs_tolerance;
            TEST_ASSERT_MSG(
                check_diff(got_tidal_power,
                           expected_tidal_power,
                           0.0,
                           abs_tolerance),
                message.str().c_str()
            );
            TEST_ASSERT_MSG(
                check_diff(got_x_torque,
                           0.0,
                           0.0,
                           1e-15),
                "Non-zero x-torque for aligned zone!"
            );
            TEST_ASSERT_MSG(
                check_diff(got_y_torque,
                           0.0,
                           0.0,
                           1e-15),
                "Non-zero y-torque for aligned zone!"
            );
            message.str("");
            message << "For e = " << e
                    << " tidal Z torque = " << got_z_torque
                    << ", expected: " << expected_z_torque
                    << "; difference: " << (got_z_torque
                                            -
                                            expected_z_torque)
                    << " > " << abs_tolerance;
            TEST_ASSERT_MSG(
                check_diff(got_z_torque,
                           expected_z_torque,
                           0.0,
                           abs_tolerance),
                message.str().c_str()
            );
        }
    }

    test_TidalTorquePower::test_TidalTorquePower()
    {
        TEST_ADD(test_TidalTorquePower::test_aligned_orbit);
    }
}//End Evolve namespace.
