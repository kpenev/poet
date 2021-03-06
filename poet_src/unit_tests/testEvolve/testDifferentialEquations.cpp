/**\file
 *
 * \brief Implement the non-inline methods of test_DifferentialEquations
 *
 * \ingroup UnitTests_group
 */

#include "testDifferentialEquations.h"

namespace Evolve {

    void test_DifferentialEquations::output_rates(
        const alglib::real_1d_array &eccentricities,
        const alglib::real_1d_array &expected_semimajor_rate,
        const alglib::real_1d_array &predicted_semimajor_rate,
        const alglib::real_1d_array &expected_eccentricity_rate,
        const alglib::real_1d_array &predicted_eccentricity_rate
    ) const
    {
        std::cout << std::setw(25) << "eccentricity"
                  << std::setw(25) << "semimajor_rate"
                  << std::setw(25) << "expected_semimajor_rate"
                  << std::setw(25) << "eccentricity_rate"
                  << std::setw(25) << "expected_ecc_rate"
                  << std::endl;
        for(
            unsigned e_index = 0;
            e_index < eccentricities.length();
            ++e_index
        ) {
                std::cout << std::setw(25) << eccentricities[e_index]
                          << std::setw(25)
                          << predicted_semimajor_rate[e_index]
                          << std::setw(25)
                          << expected_semimajor_rate[e_index]
                          << std::setw(25)
                          << predicted_eccentricity_rate[e_index]
                          << std::setw(25)
                          << expected_eccentricity_rate[e_index]
                          << std::endl;
        }


    }

    double test_DifferentialEquations::zahn77_semimajor_rate_coef(
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        double eccentricity
    ) const
    {

        double e2 = std::pow(eccentricity, 2);

        if(orbital_frequency_multiplier == 2 && spin_frequency_multiplier == 2)
            return 1.0 - 5.0 * e2;

        if(orbital_frequency_multiplier == 1 && spin_frequency_multiplier == 0)
            return 3.0 / 4.0 * e2;

        if(orbital_frequency_multiplier == 1 && spin_frequency_multiplier == 2)
            return 1.0 / 8.0 * e2;

        if(orbital_frequency_multiplier == 3 && spin_frequency_multiplier == 2)
            return 147.0 / 8.0 * e2;

        return 0.0;
    }

    double test_DifferentialEquations::zahn77_eccentricity_rate_coef(
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier
    ) const
    {

        if(orbital_frequency_multiplier == 2 && spin_frequency_multiplier == 2)
            return -1.0;

        if(orbital_frequency_multiplier == 1 && spin_frequency_multiplier == 0)
            return 1.5;

        if(orbital_frequency_multiplier == 1 && spin_frequency_multiplier == 2)
            return -1.0 / 4.0;

        if(orbital_frequency_multiplier == 3 && spin_frequency_multiplier == 2)
            return 49.0 / 4.0;

        return 0.0;
    }

    void test_DifferentialEquations::check_agreement(
        const alglib::real_1d_array& x,
        const alglib::real_1d_array& y1,
        const alglib::real_1d_array& y2,
        unsigned agreement_order,
        unsigned max_order,
        const std::string &description
    )
    {
        alglib::real_1d_array difference;
        difference.setlength(x.length());
        double min_difference = Core::Inf,
               max_difference = -Core::Inf,
               y_scale = 0;
        for(unsigned i = 0; i < x.length(); ++i) {
            difference[i] = y2[i] - y1[i];
            min_difference = std::min(difference[i], min_difference);
            max_difference = std::max(difference[i], max_difference);
            y_scale = std::max(std::abs(y1[i]), std::abs(y2[i]));
        }

        if(
            std::max(std::abs(max_difference), std::abs(min_difference))
            <
            1e-13 * y_scale
        )
            return;
        alglib::ae_int_t fit_info;
        alglib::barycentricinterpolant interp;
        alglib::polynomialfitreport fit_report;
        alglib::polynomialfit(x,
                              difference,
                              max_order + 1,
                              fit_info,
                              interp,
                              fit_report);

        double max_allowed_residual = 1e-12 * (max_difference - min_difference);

        std::ostringstream message;
        message << description
                << " too large residual from fit to difference: "
                << fit_report.maxerror
                << " > "
                << max_allowed_residual;

        TEST_ASSERT_MSG(fit_report.maxerror <= max_allowed_residual,
                        message.str().c_str());

        alglib::real_1d_array coefficients;

        alglib::polynomialbar2pow(interp, coefficients);

        double max_allowed_term = 0.0,
               max_abs_x = std::max(std::abs(x[0]), x[x.length() - 1]);
        for(unsigned power = 0; power < coefficients.length(); ++power)
            max_allowed_term = std::max(
                max_allowed_term,
                std::abs(coefficients[power]) * std::pow(max_abs_x, power)
            );
        max_allowed_term *= 1e-8;

        for(unsigned power = 0; power <= agreement_order; ++power) {
            double max_abs_term_value = (std::abs(coefficients[power])
                                         *
                                         std::pow(max_abs_x, power));
            message.str("");
            message << description
                    << " x^" << power << " coefficient too large: |"
                    << coefficients[power]
                    << " * "
                    << max_abs_x << "^" << power
                    << "| = "
                    << max_abs_term_value
                    << " > "
                    << max_allowed_term;

            TEST_ASSERT_MSG(max_abs_term_value <= max_allowed_term,
                            message.str().c_str());
        }
    }

    test_DifferentialEquations::test_DifferentialEquations()
    {
        TEST_ADD(test_DifferentialEquations::test_aligned_equations);
        TEST_ADD(test_DifferentialEquations::test_error_estimate);
    }

    void test_DifferentialEquations::test_aligned_equations()
    {
        const double MAX_ECCENTRICITY  = 0.5;
        const unsigned NUM_ECCENTRICITIES = 100;
        const double ECCENTRICITY_STEP = (
            MAX_ECCENTRICITY / (NUM_ECCENTRICITIES - 1)
        );

        alglib::real_1d_array eccentricities,
                              expected_semimajor_rate,
                              predicted_semimajor_rate,
                              expected_eccentricity_rate,
                              predicted_eccentricity_rate;
        eccentricities.setlength(NUM_ECCENTRICITIES);
        expected_semimajor_rate.setlength(NUM_ECCENTRICITIES);
        predicted_semimajor_rate.setlength(NUM_ECCENTRICITIES);
        expected_eccentricity_rate.setlength(NUM_ECCENTRICITIES);
        predicted_eccentricity_rate.setlength(NUM_ECCENTRICITIES);

        eccentricities[0] = 0.0;
        for(unsigned e_index = 1; e_index < NUM_ECCENTRICITIES; ++e_index)
            eccentricities[e_index] = (eccentricities[e_index - 1]
                                       +
                                       ECCENTRICITY_STEP);

        std::valarray<double> one_poly(1.0, 1);
        for(
            int orbital_frequency_multiplier = 0;
            orbital_frequency_multiplier <= 4;
            ++orbital_frequency_multiplier
        ) {
            for(
                int spin_frequency_multiplier = 0;
                spin_frequency_multiplier <= 4;
                ++spin_frequency_multiplier
            ) {
                SingleTidalTermBody primary(0.0,
                                            1.0,
                                            orbital_frequency_multiplier,
                                            spin_frequency_multiplier,
                                            1.0,
                                            one_poly,
                                            one_poly,
                                            one_poly);
                Planet::Planet secondary(1.0, 1.0);

                BinarySystem binary(primary, secondary);

                binary.change_e_order(2);

                std::valarray<double> parameters(0.0, 7);
                std::valarray<double> evolution_rates(parameters.size());
                parameters[0] = 10.0;
                parameters[1] = 0.0;
                binary.configure(true,
                                 0.0,
                                 &(parameters[0]),
                                 Core::BINARY);

                double common_rate_factor = (
                    -2.4 * M_PI
                    *
                    (primary.mass() + secondary.mass()) / primary.mass()
                    *
                    Core::AstroConst::G
                    *
                    Core::AstroConst::solar_mass * secondary.mass()
                    /
                    std::pow(Core::AstroConst::solar_radius, 3)
                ) * Core::AstroConst::Gyr;

                for(double semimajor = 2.0; semimajor <= 10.0; semimajor += 1.0) {

                    double mean_motion = std::sqrt(
                        Core::AstroConst::G * Core::AstroConst::solar_mass * (
                            primary.mass()
                            +
                            secondary.mass()
                        )
                        /
                        std::pow(semimajor * Core::AstroConst::solar_radius, 3)
                    );
                    double rate_factor = (
                        common_rate_factor
                        /
                        (mean_motion * std::pow(semimajor, 8))
                    );
                    for(
                        double primary_spin_frequency = 0.0;
                        primary_spin_frequency <= 5.0 * mean_motion;
                        primary_spin_frequency += 0.2 * mean_motion
                    ) {

                        for(
                            unsigned e_index = 0;
                            e_index < NUM_ECCENTRICITIES;
                            ++e_index
                        )
                        {
                            double age = 1.0,
                                   test_e = eccentricities[e_index];
                            parameters[0] = std::pow(semimajor, 6.5);
                            parameters[1] = test_e;
                            binary.differential_equations(age,
                                                          &(parameters[0]),
                                                          Core::BINARY,
                                                          &(evolution_rates[0]));
                            expected_semimajor_rate[e_index] = (
                                6.5 * std::pow(semimajor, 6.5)
                                *
                                rate_factor
                                *
                                zahn77_semimajor_rate_coef(
                                    orbital_frequency_multiplier,
                                    spin_frequency_multiplier,
                                    test_e
                                )
                            );
                            expected_eccentricity_rate[e_index] = (
                                test_e
                                *
                                rate_factor / 4.0
                                *
                                zahn77_eccentricity_rate_coef(
                                    orbital_frequency_multiplier,
                                    spin_frequency_multiplier
                                )
                            );
                            predicted_semimajor_rate[e_index] = evolution_rates[0];
                            predicted_eccentricity_rate[e_index] =
                                evolution_rates[1];
                        }
/*                        output_rates(eccentricities,
                                     expected_semimajor_rate,
                                     predicted_semimajor_rate,
                                     expected_eccentricity_rate,
                                     predicted_eccentricity_rate);*/
                        std::ostringstream message;
                        message << "orbital freq. mult. = "
                                << orbital_frequency_multiplier
                                << ", spin freq. mult. = "
                                << spin_frequency_multiplier
                                << " for a = " << semimajor
                                << ", W* = " << primary_spin_frequency;
                        check_agreement(
                            eccentricities,
                            expected_semimajor_rate,
                            predicted_semimajor_rate,
                            2,
                            4,
                            (
                                "Checking semimajor rate against Zahn (1977)"
                                +
                                message.str()
                            )
                        );
                        check_agreement(
                            eccentricities,
                            expected_eccentricity_rate,
                            predicted_eccentricity_rate,
                            2,
                            20,
                            (
                                "Checking eccentricity rate against Zahn (1977)"
                                +
                                message.str()
                            )
                        );
                    }
                }
            }
        }
    }

    void test_DifferentialEquations::test_error_estimate()
    {
        const double MAX_ECCENTRICITY = 0.9;
        const unsigned NUM_ECCENTRICITIES = 100;
        const double ECCENTRICITY_STEP = (
            MAX_ECCENTRICITY / (NUM_ECCENTRICITIES - 1)
        );
        const double a = 10.0;
        const double age = 1.0;

        StellarEvolution::MockStellarEvolution *no_evol =
            StellarEvolution::make_no_evolution();

        Star::InterpolatedEvolutionStar *star = make_const_lag_star(
            *no_evol,
            1.0,
            1.0,
            1.0,
            1.0
        );

        Planet::Planet planet(1.0, 1.0);

        BinarySystem binary(*star, planet);

        std::valarray<double> parameters(0.0, 10),
                              rough_rates(parameters.size()),
                              fine_rates(parameters.size()),
                              rate_errors(parameters.size());

        parameters[0] = a;

        binary.configure(true,
                         (MIN_AGE + MAX_AGE) / 2.0,
                         &(parameters[0]),
                         Core::BINARY);
        star->detect_saturation();

        std::cout << std::setw(25) << "e_order"
                  << std::setw(25) << "e"
                  << std::setw(25) << "rough_a_rate"
                  << std::setw(25) << "fine_a_rate"
                  << std::setw(25) << "a_rate_error"
                  << std::setw(25) << "rough_e_rate"
                  << std::setw(25) << "fine_e_rate"
                  << std::setw(25) << "e_rate_error"
                  << std::setw(25) << "rough_inc_rate"
                  << std::setw(25) << "fine_inc_rate"
                  << std::setw(25) << "inc_rate_error"
                  << std::setw(25) << "rough_spin_rate"
                  << std::setw(25) << "fine_spin_rate"
                  << std::setw(25) << "spin_rate_error"

                  << std::endl;


        for(unsigned e_order = 5; e_order <= 100; e_order+=5) {
            for(double e = 0.0; e <= MAX_ECCENTRICITY; e += ECCENTRICITY_STEP) {
                parameters[1] = e;

                star->zone(0).change_e_order(e_order, binary, true, 0);
                binary.differential_equations(age,
                                              &(parameters[0]),
                                              Core::BINARY,
                                              &(rough_rates[0]));
                binary.differential_equations(age,
                                              &(parameters[0]),
                                              Core::BINARY,
                                              &(rate_errors[0]),
                                              true);

                star->zone(0).change_e_order(2 * e_order, binary, true, 0);
                binary.differential_equations(age,
                                              &(parameters[0]),
                                              Core::BINARY,
                                              &(fine_rates[0]));

                std::cout << std::setw(25) << e_order
                          << std::setw(25) << e
                          << std::setw(25) << rough_rates[0]
                          << std::setw(25) << fine_rates[0]
                          << std::setw(25) << rate_errors[0]
                          << std::setw(25) << rough_rates[1]
                          << std::setw(25) << fine_rates[1]
                          << std::setw(25) << rate_errors[1]
                          << std::setw(25) << rough_rates[2]
                          << std::setw(25) << fine_rates[2]
                          << std::setw(25) << rate_errors[2]
                          << std::setw(25) << rough_rates[7]
                          << std::setw(25) << fine_rates[7]
                          << std::setw(25) << rate_errors[7]
                          << std::endl;
            }
        }

        delete star;
        delete no_evol;
    }

} //End Envolve namespace.
