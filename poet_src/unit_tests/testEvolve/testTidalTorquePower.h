/**\file
 *
 * \brief Declare a unit tests class that check the calculations of tidal
 * torque and power.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_TIDAL_TORQUE_POWER_H
#define __TEST_TIDAL_TORQUE_POWER_H

#include "ConstPhaseLagZone.h"
#include "../../Planet/LockedPlanet.h"
#include "../../Evolve/BinarySystem.h"

#include "../shared/Common.h"

#include <iostream>
#include <iomanip>

namespace Evolve {

    /**\brief The test suite that compares the tidal torque and power for a
     * single zone against expectations.
     *
     * \ingroup UnitTests_group
     */
    class test_TidalTorquePower : public Test::Suite {
    private:
    public:
        ///Create the test suite.
        test_TidalTorquePower();

        ///\brief Test the dimensionless tidal torque and power for an
        ///equatorial orbit, as a function of eccentricity against Zahn77 (eq.
        ///3.6 and eq 3.8).
        void test_aligned_orbit();

        ///\brief Output a table showing the convergence as e-order is
        ///increased for aligned orbit and const phase lag.
        void test_convergence();
    };//End testTidalTorquePower class.


} //End Evolve namespace.

#endif
