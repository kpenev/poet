/**\file
 *
 * \brief Unit tests that check the differential equations for eccentricity and
 * semimajor against analytic expressions.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_DIFFERENTIAL_EQUATIONS_H
#define __TEST_DIFFERENTIAL_EQUATIONS_H

#include "SingleTidalTermBody.h"
#include "../shared/Common.h"
#include "../shared/PolynomialEvolution.h"
#include "../shared/MakeStar.h"
#include "../../Planet/Planet.h"
#include "../../third_party_libs/alglib/alglib/src/interpolation.h"
#include <iomanip>

namespace Evolve {
    /**\brief The test suite that compares the differential equations for
     * eccentricity and semimajor axis to literature expansions.
     *
     * \ingroup UnitTests_group
     */
    class test_DifferentialEquations : public Test::Suite {
    private:
        ///\brief Output the predicted/expected semiamjor and eccentricity
        ///evolution rates to stdout.
        void output_rates(
            const alglib::real_1d_array &eccentricities,
            const alglib::real_1d_array &expected_semimajor_rate,
            const alglib::real_1d_array &predicted_semimajor_rate,
            const alglib::real_1d_array &expected_eccentricity_rate,
            const alglib::real_1d_array &predicted_eccentricity_rate
        ) const;


        ///\brief The Zahn (1977) coefficient for the semimajor evolution rate
        ///that corresponds to the given tidal term.
        double zahn77_semimajor_rate_coef(
            ///The multiple with which the orbital frequnecy enters in the tidal
            ///frequency expression for the selected term.
            int orbital_frequency_multiplier,

            ///The multiple with which the spin frequnecy enters in the tidal
            ///frequency expression for the selected term.
            int spin_frequency_multiplier,

            ///The eccentricity of the orbit.
            double eccentricity
        ) const;

        ///\brief Same as zahn77_semimajor_rate_coef() but for the eccentricity
        ///evolution.
        double zahn77_eccentricity_rate_coef(
            ///The multiple with which the orbital frequnecy enters in the tidal
            ///frequency expression for the selected term.
            int orbital_frequency_multiplier,

            ///The multiple with which the spin frequnecy enters in the tidal
            ///frequency expression for the selected term.
            int spin_frequency_multiplier
        ) const;

        ///\brief Check if two given dependencies on x agree up to a given order
        ///in x
        void check_agreement(
            ///The independent variable.
            const alglib::real_1d_array& x,

            ///The first dependence on x to compare.
            const alglib::real_1d_array& y1,

            ///The second dependence on x to compare.
            const alglib::real_1d_array& y2,

            ///The order up to which the two dependneces must agree.
            unsigned agreement_order,

            ///The maximum order of terms included in either dependence.
            unsigned max_order,

            ///A message describing what is being compared.
            const std::string &description
        );

    protected:
        ///No fixtures at this time
        void setup() {};

        ///No fixtures at this time
        void tear_down() {};
    public:
        ///Create the test suite.
        test_DifferentialEquations();

        ///\brief Test the a & e differential equations for aligned orbit.
        void test_aligned_equations();

        ///\brief Test the error estimate of the differential equations.
        void test_error_estimate();
    };// End test_DifferentialEquations class.

} //End Evolve namespace.

#endif
