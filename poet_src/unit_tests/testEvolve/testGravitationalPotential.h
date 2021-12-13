/**\file
 *
 * \brief Unit tests that check the expansion of the gravitational potential vs.
 * analytic expressions.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_GRAVITATIONAL_POTENTIAL_H
#define __TEST_GRAVITATIONAL_POTENTIAL_H

#include "EccentricOrbit.h"
#include "TidalPotential.h"
#include "TidalPotentialExpansion.h"

#include "../shared/Common.h"

#include <iostream>
#include <iomanip>
#include <sstream>

namespace Evolve {
    /**\brief The test suite that compares the tidal potential expansion to
     * the exact expression.
     *
     * \ingroup UnitTests_group
     */
    class test_GravitationalPotential : public Test::Suite {
    private:
        ///\brief Print to stdout the location of the secondary relative to the
        ///primary as a function of time for the given orbit.
        void print_orbit(const EccentricOrbit &orbit) const;

        ///\brief Print to stdout the value of the tidal potential without an
        ///approximation and using the expansion as a function of time.
        void print_tidal_potential(
            ///See same name argument to TidalPotential constructor.
            double primary_mass,

            ///See same name argument to TidalPotential constructor.
            double secondary_mass,

            ///See same name argument to TidalPotential constructor.
            double semimajor,

            ///See same name argument to TidalPotential constructor.
            double eccentricity,

            ///See same name argument to TidalPotential constructor.
            double inclination,

            ///See same name argument to TidalPotential constructor.
            double arg_of_periapsis,

            ///The position where to test for agreement between the two
            ///potentials, in the coordinate system expected by the potentials..
            const Eigen::Vector3d &position,

            ///The eccentricity expansion order to use.
            unsigned expansion_order
        ) const;

        ///The expected absolute precision in the potential expansion for the
        ///given orbit.
        double abs_precision_scale(
            ///The position where the potential expansion will be evaluated.
            const Eigen::Vector3d &position,

            ///The orbit for which the potential expansion is being tested.
            const EccentricOrbit &orbit
        ) const;

        ///\brief Test the expansion for a given system at a given position,
        ///sampling time.
        void test_single_point(
            ///The tidal potential expansion to test.
            TidalPotentialExpansion &approx_potential,

            ///The exact potential expected for the system.
            const TidalPotential &exact_potential,

            ///The position where to test for agreement between the two
            ///potentials, in the coordinate system expected by the potentials..
            const Eigen::Vector3d &position
        );

        ///Test the expansion for a given system, sampling positions and times.
        void test_system(
            ///See same name argument to TidalPotential constructor.
            double primary_mass,

            ///See same name argument to TidalPotential constructor.
            double secondary_mass,

            ///See same name argument to TidalPotential constructor.
            double semimajor,

            ///See same name argument to TidalPotential constructor.
            double eccentricity,

            ///See same name argument to TidalPotential constructor.
            double inclination,

            ///See same name argument to TidalPotential constructor.
            double arg_of_periapsis
        );

    protected:
        ///No fixtures at this time
        void setup() {};

        ///No fixtures at this time
        void tear_down() {};

    public:
        ///Create the test suite.
        test_GravitationalPotential();

        ///\brief Test the expansion of the potential for multiple system
        ///configurations, locations and times.
        void test_expansion();
    }; //End test_GravitationalPotential class.
} //End Evolve namespace.

#endif
