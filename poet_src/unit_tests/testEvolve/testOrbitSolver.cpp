/**\file
 *
 * \brief Defines some of the methods of the test suite that exercises the
 * OrbitSolver class and the other clasess necessary to accomplish this.
 *
 * \ingroup UnitTests_group
 */

#include "testOrbitSolver.h"

namespace Evolve {

    const double TSTART = 2.0 * MIN_AGE;
    const double zero = 0.0;
    const double one = 1.0;

    const double Mjup_to_Msun = (Core::AstroConst::jupiter_mass
                                 /
                                 Core::AstroConst::solar_mass);

    const double Rjup_to_Rsun = (Core::AstroConst::jupiter_radius
                                 /
                                 Core::AstroConst::solar_radius);

    std::ostream &operator<<(std::ostream &os,
                             RealEvolutionQuantity q)
    {
        switch(q) {
            case SEMIMAJOR :
                os << "SEMIMAJOR"; break;
            case ECCENTRICITY:
                os << "ECCENTRICITY"; break;
            case CONV_INCLINATION:
                os << "CONV_INCLINATION"; break;
            case RAD_INCLINATION:
                os << "RAD_INCLINATION"; break;
            case CONV_PERIAPSIS:
                os << "CONV_PERIAPSIS"; break;
            case RAD_PERIAPSIS:
                os << "RAD_PERIAPSIS"; break;
            case CONV_ANGMOM:
                os << "CONV_ANGMOM"; break;
            case RAD_ANGMOM:
                os << "RAD_ANGMOM"; break;
            case AGE:
                os << "AGE"; break;
            default :
                assert(false);
        }
        return os;
    }


    StellarEvolution::PolynomialEvolutionQuantity nan_func(
        std::valarray<double>(Core::NaN, 2),
        TSTART,
        MAX_AGE
    );

    StellarEvolution::PolynomialEvolutionQuantity zero_func(
        std::valarray<double>(),
        TSTART,
        MAX_AGE
    );

    StellarEvolution::PolynomialEvolutionQuantity one_func(
        std::valarray<double>(1.0, 1),
        TSTART,
        MAX_AGE
    );

    StellarEvolution::PolynomialEvolutionQuantity two_func(
        std::valarray<double>(2.0, 1),
        TSTART,
        MAX_AGE
    );

    StellarEvolution::PolynomialEvolutionQuantity two_hundred_func(
        std::valarray<double>(200.0, 1),
        TSTART,
        MAX_AGE
    );

    void test_OrbitSolver::set_single_component_dissipation(
        double min_frequency,
        double max_frequency,
        double decay_scale,
        double phase_lag
    )
    {
        //phase lag(min_frequency - decay_scale)
        //=
        //phase lag(max_frequency + decay_scale)
        //=
        //suppression_factor * phase_lag
        const double suppression_factor = 0.01;

        std::vector<double> breaks(2);
        breaks[0] = min_frequency;
        breaks[1] = max_frequency;

        std::vector<double> powerlaw_indices(3);
        powerlaw_indices[0] = (
            std::log(suppression_factor)
            /
            std::log(1.0 - decay_scale / min_frequency)
        );
        powerlaw_indices[1] = 0.0;
        powerlaw_indices[2] = (
            std::log(suppression_factor)
            /
            std::log(1.0 + decay_scale / max_frequency)
        );

        BrokenPowerlawPhaseLagZone *zone;

        if(__star) {
            assert(__primary_planet == NULL);
            zone = &(__star->envelope());
        } else {
            assert(__primary_planet);
            zone = &(__primary_planet->zone());
        }

        zone->setup(breaks,
                   std::vector<double>(),
                   powerlaw_indices,
                   std::vector<double>(1, 0.0),
                   phase_lag);
    }

    void test_OrbitSolver::make_single_component_star(
        const StellarEvolution::Interpolator &evolution,
        double wind_strength,
        double wind_sat_freq,
        double coupling_timescale,
        double min_frequency,
        double max_frequency,
        double decay_scale,
        double phase_lag
    )
    {
        __star = new Star::InterpolatedEvolutionStar(1.0,//mass
                                                     0.0,//feh
                                                     wind_strength,
                                                     wind_sat_freq,
                                                     coupling_timescale,
                                                     evolution);
        set_single_component_dissipation(min_frequency,
                                         max_frequency,
                                         decay_scale,
                                         phase_lag);

    }

    void test_OrbitSolver::evolve(double wdisk,
                                  double tdisk,
                                  double initial_a,
                                  const double *initial_Lstar,
                                  double initial_incl,
                                  double secondary_mass,
                                  double tsecondary,
                                  double max_age,
                                  double secondary_radius,
                                  double precision,
                                  double max_step_factor)
    {
        Evolve::DissipatingBody *primary;
       if(__star)
           primary = __star;
       else
           primary = __primary_planet;

        secondary_mass *= Mjup_to_Msun;
        secondary_radius *= Rjup_to_Rsun;

        if(std::isnan(tsecondary)) tsecondary = tdisk;

        Planet::Planet secondary(secondary_mass,
                                 (secondary_mass ? secondary_radius : 0.0));
        secondary.configure(true, //init
                            tsecondary, //age
                            primary->mass(), //mass
                            initial_a, //semimajor
                            0.0, //eccentricity
                            &zero, //spin angmom
                            NULL, //inclination
                            NULL, //periapsis
                            false, //locked surface
                            true, //zero outer inclination
                            true);//zero outer periapsis

        __system = new Evolve::DiskBinarySystem(
            *primary,
            secondary,
            initial_a, //semimajor
            0.0, //eccentricity
            initial_incl, //inclination
            wdisk, //Wdisk
            tdisk, //disk dissipation age
            tsecondary //secondary formation age
        );
        double zeros[] = {0.0, 0.0};
        if(__star) {
            if(tdisk <= TSTART) {
                if(tsecondary <= TSTART) {
                    __system->configure(true, //init
                                        TSTART,
                                        initial_a, //semimajor
                                        0.0, //eccentricity
                                        initial_Lstar, //spin angmom
                                        zeros, //inclination
                                        zeros, //periapsis
                                        Core::BINARY);
                } else {
                    __system->configure(true, //init
                                        TSTART,
                                        Core::NaN, //semimajor
                                        Core::NaN, //eccentricity
                                        initial_Lstar, //spin angmom
                                        zeros, //inclination
                                        zeros, //periapsis
                                        Core::SINGLE);
                }
            } else {
                __system->configure(true, //init
                                    TSTART,
                                    Core::NaN, //semimajor
                                    Core::NaN, //eccentricity
                                    initial_Lstar, //spin angmom
                                    NULL, //inclination
                                    NULL, //periapsis
                                    Core::LOCKED_SURFACE_SPIN);
            }
            __star->detect_saturation();
        } else {
            double initial_inclinations[] = {initial_incl, 0.0, 0.0};
            __system->configure(true, //init
                                tsecondary,
                                initial_a, //semimajor
                                0.0, //eccentricity
                                initial_Lstar, //spin angmom
                                initial_inclinations, //inclination
                                zeros, //periapsis
                                Core::BINARY);
        }
        __solver = new Evolve::OrbitSolver(max_age, precision);
        (*__solver)(*__system,
                    (max_age - __system->age()) * max_step_factor,//time step
                    std::list<double>()); //no required ages
    }

    std::vector< const std::list<double> *> test_OrbitSolver::get_evolution()
        const
    {
        std::vector< const std::list<double> *>
            tabulated_real_quantities(NUM_REAL_QUANTITIES);

            tabulated_real_quantities[AGE] = &(__solver->evolution_ages());

            tabulated_real_quantities[SEMIMAJOR] =
                &(__system->semimajor_evolution());

            tabulated_real_quantities[ECCENTRICITY] =
                &(__system->eccentricity_evolution());

            Evolve::DissipatingZone *primary_envelope;

            if(__star) {
                tabulated_real_quantities[RAD_INCLINATION] =
                    &(__star->core().get_evolution_real(Evolve::INCLINATION));

                tabulated_real_quantities[RAD_PERIAPSIS] =
                    &(__star->core().get_evolution_real(Evolve::PERIAPSIS));

                tabulated_real_quantities[RAD_ANGMOM] = &(
                    __star->core().get_evolution_real(Evolve::ANGULAR_MOMENTUM)
                );
                primary_envelope = &(__star->envelope());
            } else {
                tabulated_real_quantities[RAD_INCLINATION] = NULL;
                tabulated_real_quantities[RAD_PERIAPSIS] = NULL;
                tabulated_real_quantities[RAD_ANGMOM] = NULL;
                primary_envelope = &(__primary_planet->zone());
            }

            tabulated_real_quantities[CONV_INCLINATION] = &(
                primary_envelope->get_evolution_real(Evolve::INCLINATION)
            );

            tabulated_real_quantities[CONV_PERIAPSIS] = &(
                primary_envelope->get_evolution_real(Evolve::PERIAPSIS)
            );

            tabulated_real_quantities[CONV_ANGMOM] = &(
                primary_envelope->get_evolution_real(
                    Evolve::ANGULAR_MOMENTUM
                )
            );
        return tabulated_real_quantities;
    }

    void test_OrbitSolver::test_solution(
        const std::vector< const std::list<double> * > &
            tabulated_real_quantities,
        std::vector<const Core::OneArgumentDiffFunction *>
            expected_real_quantities,
        const ExpectedEvolutionMode<Core::EvolModeType> &expected_evol_mode,
        const ExpectedEvolutionMode<bool> &expected_wind_mode,
        double min_age,
        double max_age,
        bool debug_mode
    )
    {
        const std::list<Core::EvolModeType> &tabulated_modes =
            __solver->mode_evolution();

        const std::list<bool> *tabulated_wind_sat =
            (__star ? &(__star->wind_saturation_evolution()) : NULL);

        unsigned num_ages = tabulated_real_quantities[AGE]->size();

        std::ostringstream msg_start;
        std::ostringstream msg;
        msg.precision(16);
        msg.setf(std::ios_base::scientific);
        msg_start.precision(16);
        msg_start.setf(std::ios_base::scientific);
        msg << msg_start.str()
            << num_ages
            << " tabulated ages, ";
        bool all_same_size = true;

        for(unsigned q = 0; q < AGE; ++q) {
            if(tabulated_real_quantities[q]) {
                msg << tabulated_real_quantities[q]->size()
                    << " tabulated "
                    << static_cast<RealEvolutionQuantity>(q)
                    << ", ";
                all_same_size = (
                    all_same_size
                    &&
                    tabulated_real_quantities[q]->size() == num_ages
                );
            }
        }

        msg << tabulated_modes.size() << " tabulated modes";

        all_same_size = (all_same_size
                         &&
                         tabulated_modes.size() == num_ages);
        if(__star) {
            msg << ", "
                << tabulated_wind_sat->size()
                << " tabulated wind saturations";
            all_same_size = (all_same_size
                             &&
                             tabulated_wind_sat->size() == num_ages);
        }


        if(debug_mode) std::cout << msg.str() << std::endl;
        TEST_ASSERT_MSG(all_same_size, msg.str().c_str());

        msg.str("");
        msg << "Expected age range: " << min_age << " to " << max_age
            << ", actual age range: " << tabulated_real_quantities[AGE]->front()
            << " to " << tabulated_real_quantities[AGE]->back();
        if(debug_mode) std::cout << msg.str() << std::endl;
        TEST_ASSERT_MSG(
            (
                tabulated_real_quantities[AGE]->front() == min_age
                &&
                tabulated_real_quantities[AGE]->back() == max_age
            ),
            msg.str().c_str()
        );

        std::vector< std::list<double>::const_iterator >
            real_tabulated_iter(AGE);
        for(unsigned q = 0; q < AGE; ++q)
            if(tabulated_real_quantities[q])
                real_tabulated_iter[q] = tabulated_real_quantities[q]->begin();
        std::list<Core::EvolModeType>::const_iterator
            tabulated_mode_iter = tabulated_modes.begin();
        std::list<bool>::const_iterator tabulated_wind_sat_iter;
        if(__star)
            tabulated_wind_sat_iter = tabulated_wind_sat->begin();
        double last_checked_age = Core::NaN;

        for(
            std::list<double>::const_iterator
                age_i = tabulated_real_quantities[AGE]->begin();
            age_i != tabulated_real_quantities[AGE]->end();
            ++age_i
        ) {
            std::vector<double> expected_real_values(AGE);
            for(unsigned q = 0; q < AGE; ++q) {
                expected_real_values[q] =
                    (*(expected_real_quantities[q]))(*age_i);
            }
            Core::EvolModeType expected_mode = expected_evol_mode(*age_i);
            bool expected_wind_sat = expected_wind_mode(*age_i);

            std::ostringstream age_msg_start;
            age_msg_start.precision(16);
            age_msg_start.setf(std::ios_base::scientific);
            age_msg_start << msg_start.str()
                          << "age = " << *age_i
                          << ", mode = " << *tabulated_mode_iter;
            if(__star) {
                age_msg_start << ", wind is ";
                if(!(*tabulated_wind_sat_iter)) age_msg_start << " not ";
                age_msg_start << "saturated";
            }
            for(unsigned q = 0; q < AGE; ++q)
                if(tabulated_real_quantities[q])
                    age_msg_start << ", "
                                  << static_cast<RealEvolutionQuantity>(q)
                                  << " = "
                                  << *real_tabulated_iter[q];

            msg.str("");
            msg << age_msg_start.str() << " age is out of range.";
            TEST_ASSERT_MSG(*age_i >= min_age && *age_i <= max_age,
                            msg.str().c_str());

            std::list<double>::const_iterator next_age_i = age_i;
            ++next_age_i;
            bool
                can_skip = (
                    next_age_i != tabulated_real_quantities[AGE]->end()
                    &&
                    std::abs(*next_age_i - *age_i) < 1e-5
                    &&
                    expected_evol_mode.near_break(*age_i)
                    &&
                    expected_evol_mode.near_break(*next_age_i)
                ),
                skipped = (*age_i == last_checked_age);

            msg.str("");
            msg << age_msg_start.str() << ": mode is not "
                << expected_mode << ", but " << *tabulated_mode_iter;

            if(debug_mode) std::cout << msg.str() << std::endl;
            if(!skipped && expected_mode != *tabulated_mode_iter) {
                if(can_skip) skipped = true;
                else TEST_ASSERT_MSG(false, msg.str().c_str());
            }

            if(__star) {
                msg.str("");
                msg << age_msg_start.str() << ": wind is ";
                if(!(*tabulated_wind_sat_iter)) msg << " not ";
                msg << "saturated, but should";
                if(!expected_wind_sat) msg << " not ";
                msg << "be.";
                if(debug_mode) std::cout << msg.str() << std::endl;
                if(!skipped && expected_wind_sat != *tabulated_wind_sat_iter) {
                    if(can_skip) skipped = true;
                    else TEST_ASSERT_MSG(false, msg.str().c_str());
                }
            }
/*            TEST_ASSERT_MSG(
                (
                    skipped
                    ||
                    expected_wind_sat == *tabulated_wind_sat_iter
                    ||
                    expected_wind_mode.near_break(*age_i)
                ),
                msg.str().c_str()
            );*/

            for(unsigned q = 0; q < AGE; ++q) {
                if(!tabulated_real_quantities[q]) continue;
                msg.str("");
                msg << age_msg_start.str() << ": "
                    << static_cast<RealEvolutionQuantity>(q)
                    << " is not "
                    << expected_real_values[q]
                    << " but "
                    << *real_tabulated_iter[q]
                    << ", difference = "
                    << (*real_tabulated_iter[q]) - expected_real_values[q];
                if(debug_mode) std::cout << msg.str() << std::endl;
                if(
                    !skipped
                    &&
                    !check_diff((*real_tabulated_iter[q]),
                                expected_real_values[q],
                                1e-5,
                                0.0)
                ) {
                    if(can_skip) {
                        skipped = true;
                        break;
                    }
                    TEST_ASSERT_MSG(false, msg.str().c_str());
                }
            }
            if(!skipped)
                last_checked_age = *age_i;
            else if(debug_mode)
                std::cerr << "Skipped checks for t = " << *age_i
                          << std::endl;
            for(unsigned q = 0; q < AGE; ++q)
                if(tabulated_real_quantities[q])
                    ++(real_tabulated_iter[q]);
            if(__star)  ++tabulated_wind_sat_iter;
            ++tabulated_mode_iter;
        }
    }

    void test_OrbitSolver::test_no_planet_scenario(
        const StellarEvolution::Interpolator &stellar_evol,
        double *initial_Lstar,
        double windK,
        double wind_sat_freq,
        double core_env_coupling_time,
        std::vector<const Core::OneArgumentDiffFunction *>
            &expected_real_quantities,
        const ExpectedEvolutionMode<bool> &expected_wind_mode,
        double max_age,
        bool debug_mode
    ) {
            ExpectedEvolutionMode<Core::EvolModeType> single_mode,
                                                      binary_mode;
            single_mode.add_break(TSTART, Core::SINGLE);
            binary_mode.add_break(TSTART, Core::BINARY);

            __star = make_const_lag_star(stellar_evol,
                                         windK,
                                         wind_sat_freq,
                                         core_env_coupling_time,
                                         1.0);//phase lag
            evolve(0.0,//Wdisk
                   0.0,//tdisk
                   1.0,//initial semimajor
                   initial_Lstar,
                   0.0,//initial inclination
                   1.0,//planet mass
                   max_age,//planet formation age.
                   max_age);//end evolution age
            expected_real_quantities[SEMIMAJOR] = &nan_func;
            expected_real_quantities[ECCENTRICITY] = &nan_func;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;

            test_solution(get_evolution(),
                          expected_real_quantities,
                          single_mode,
                          expected_wind_mode,
                          TSTART,
                          max_age,
                          debug_mode);

            if(initial_Lstar[0] == 0) return;

            for(double phase_lag = 0.0; phase_lag < 1.5; phase_lag += 1.0)
                for(
                    double mplanet = 0.0;
                    mplanet < 1.5 - phase_lag;
                    mplanet += 1.0
                ) {
                    delete __star;
                    delete __system;
                    delete __solver;

                    __star = make_const_lag_star(stellar_evol,
                                                 windK,
                                                 wind_sat_freq,
                                                 core_env_coupling_time,
                                                 phase_lag);
                    evolve(0.0,//Wdisk
                           0.0,//tdisk
                           1.0,//initial semimajor
                           initial_Lstar,
                           0.0,//initial inclination
                           mplanet,//planet mass
                           TSTART,//planet formation age.
                           max_age);//end evolution age

                    expected_real_quantities[SEMIMAJOR] = &one_func;
                    expected_real_quantities[ECCENTRICITY] = &zero_func;
                    test_solution(get_evolution(),
                                  expected_real_quantities,
                                  binary_mode,
                                  expected_wind_mode,
                                  TSTART,
                                  max_age,
                                  debug_mode);
                    delete __system;
                    delete __solver;

                    evolve(0.0,//Wdisk
                           0.0,//tdisk
                           200.0,//initial semimajor
                           initial_Lstar,
                           0.0,//initial inclination
                           mplanet,//planet mass
                           TSTART,//planet formation age.
                           max_age);//end evolution age

                    expected_real_quantities[SEMIMAJOR] = &two_hundred_func;
                    test_solution(get_evolution(),
                                  expected_real_quantities,
                                  binary_mode,
                                  expected_wind_mode,
                                  TSTART,
                                  max_age,
                                  debug_mode);

                }
            delete __star;
            __star = NULL;
            delete __system;
            __system = NULL;
            delete __solver;
            __solver = NULL;

    }

    std::vector<const Core::OneArgumentDiffFunction *>
        test_OrbitSolver::calculate_expected_unlocked_evolution(
            double phase_lag,
            double secondary_mass,
            bool decaying
        )
        {
            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            double msecondary_si = (secondary_mass
                                    *
                                    Core::AstroConst::jupiter_mass),
                   mprimary_si = Core::AstroConst::solar_mass,
                   alpha = (
                       2.4 * M_PI
                       *
                       std::sqrt(
                           Core::AstroConst::G * (msecondary_si + mprimary_si)
                           /
                           Core::AstroConst::solar_radius
                       )
                       *
                       msecondary_si / mprimary_si
                       *
                       phase_lag * Core::AstroConst::Gyr
                       /
                       Core::AstroConst::solar_radius
                   ),
                   a6p5_offset,
                   a0,
                   Lscale = (
                       -msecondary_si
                       /
                       std::pow(Core::AstroConst::solar_radius, 1.5)
                       *
                       std::sqrt(
                           Core::AstroConst::G
                           /
                           (msecondary_si + Core::AstroConst::solar_mass)
                       )
                       *
                       Core::AstroConst::day
                   );

            if(decaying) {
                a6p5_offset = std::pow(2.0, 6.5) + 6.5 * alpha;
                a0 = std::pow(a6p5_offset, 1.0 / 6.5);
            } else {
                a0 = 2.6;
                a6p5_offset = std::pow(a0, 6.5);
            }

            std::valarray<double> a6p5_poly_coef(2);
            a6p5_poly_coef[0] = a6p5_offset;
            a6p5_poly_coef[1] = (decaying ? -1.0 : 1.0) * 6.5 * alpha;
            StellarEvolution::PolynomialEvolutionQuantity *a6p5_evol = (
                new StellarEvolution::PolynomialEvolutionQuantity(
                    a6p5_poly_coef,
                    TSTART,
                    1.0
                )
            );
            __temp_functions.push_back(a6p5_evol);

            FunctionToPower *sqrta_evol = new FunctionToPower(a6p5_evol, 1.0/13.0);
            __temp_functions.push_back(sqrta_evol);


            ExponentialPlusFunc *Lconv_unscaled = new ExponentialPlusFunc(
                sqrta_evol,
                (decaying ? 0.0 : -1e5) - std::sqrt(a0),
                0
            );
            __temp_functions.push_back(Lconv_unscaled);

            expected_real_quantities[SEMIMAJOR] = new FunctionToPower(
                a6p5_evol,
                1.0 / 6.5
            );
            __temp_functions.push_back(expected_real_quantities[SEMIMAJOR]);

            expected_real_quantities[ECCENTRICITY] = &zero_func;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = new ScaledFunction(
                Lconv_unscaled,
                Lscale
            );
            __temp_functions.push_back(expected_real_quantities[CONV_ANGMOM]);
            expected_real_quantities[RAD_ANGMOM] = &zero_func;

            return expected_real_quantities;

        }

    std::vector<const Core::OneArgumentDiffFunction *>
        test_OrbitSolver::calculate_expected_disklocked_to_fast_to_locked(
            double lgQ,
            double tdisk,
            double async,
            double tsync,
            double tend,
            bool include_disk_lock
        )
        {
            const double Q = std::pow(10.0, lgQ),
                         alpha = (
                             -4.5
                             *
                             std::sqrt(
                                 Core::AstroConst::G
                                 /
                                 (
                                     Core::AstroConst::solar_radius
                                     *
                                     Core::AstroConst::solar_mass
                                 )
                             )
                             *
                             Core::AstroConst::jupiter_mass / Q
                             *
                             Core::AstroConst::Gyr
                             /
                             Core::AstroConst::solar_radius
                         ),
                         Lscale = (
                             Core::AstroConst::jupiter_mass
                             /
                             std::pow(Core::AstroConst::solar_radius, 1.5)
                             *
                             std::sqrt(Core::AstroConst::G
                                       /
                                       (
                                           Core::AstroConst::jupiter_mass
                                           +
                                           Core::AstroConst::solar_mass
                                       )
                             )
                             *
                             Core::AstroConst::day
                         ),
                         beta = (
                             std::sqrt(
                                 Core::AstroConst::G
                                 *
                                 (
                                     Core::AstroConst::solar_mass
                                     +
                                     Core::AstroConst::jupiter_mass
                                 )
                             )
                             *
                             Core::AstroConst::day
                             /
                             std::pow(Core::AstroConst::solar_radius, 1.5)
                         ),
                         a6p5_offset = (std::pow(async, 6.5)
                                        -
                                        6.5 * alpha * tsync),
                         a_formation=std::pow(
                             a6p5_offset + 6.5 * alpha * tdisk,
                             1.0 / 6.5
                         ),
                         Ic = (
                             Lscale
                             *
                             (std::sqrt(a_formation) - std::sqrt(async))
                             /
                             (beta * (std::pow(async, -1.5)
                                      -
                                      0.5 * std::pow(a_formation, -1.5)))
                         ),
                         wdisk = 0.5 * beta / std::pow(a_formation, 1.5),
                         wlocked = beta / std::pow(async, 1.5);

            std::valarray<double> a6p5_poly_coef(2);
            a6p5_poly_coef[0] = a6p5_offset;
            a6p5_poly_coef[1] = 6.5 * alpha;
            StellarEvolution::PolynomialEvolutionQuantity
                *a6p5_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    a6p5_poly_coef,
                    tdisk,
                    tsync
                ),
                *Lconv_disk = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(Ic * wdisk, 1),
                    TSTART,
                    tdisk
                ),
                *Lconv_locked = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(Ic * wlocked, 1),
                    tsync,
                    tend
                ),
                *nan_disk = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(Core::NaN, 2),
                    TSTART,
                    tdisk
                ),
                *a_locked = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(async, 1),
                    tsync,
                    tend
                ),
                *Lrad_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(),
                    TSTART,
                    tend
                ),
                *zero_e = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(),
                    tdisk,
                    tend
                );
            __temp_functions.push_back(a6p5_evol);
            __temp_functions.push_back(Lconv_disk);
            __temp_functions.push_back(Lconv_locked);
            __temp_functions.push_back(nan_disk);
            __temp_functions.push_back(a_locked);
            __temp_functions.push_back(Lrad_evol);
            __temp_functions.push_back(zero_e);

            FunctionToPower *a_fast = new FunctionToPower(a6p5_evol,
                                                          1.0 / 6.5),
                            *sqrta_evol = new FunctionToPower(a6p5_evol,
                                                              1.0 / 13.0);
            __temp_functions.push_back(a_fast);
            __temp_functions.push_back(sqrta_evol);

            ExponentialPlusFunc *Lconv_unscaled = new ExponentialPlusFunc(
                sqrta_evol,
                -Ic * wdisk / Lscale - std::sqrt(a_formation), 0
            );
            __temp_functions.push_back(Lconv_unscaled);

            ScaledFunction *Lconv_fast = new ScaledFunction(Lconv_unscaled,
                                                            -Lscale);
            __temp_functions.push_back(Lconv_fast);

            PiecewiseFunction *a_evol = new PiecewiseFunction,
                              *e_evol = new PiecewiseFunction,
                              *Lconv_evol = new PiecewiseFunction;
            __temp_functions.push_back(a_evol);
            __temp_functions.push_back(e_evol);
            __temp_functions.push_back(Lconv_evol);

            StellarEvolution::PolynomialEvolutionQuantity *zero_quantity;

            if(include_disk_lock) {
                a_evol->add_piece(nan_disk);
                e_evol->add_piece(nan_disk);
                Lconv_evol->add_piece(Lconv_disk);
                zero_quantity = &zero_func;
            } else {
                zero_quantity = zero_e;
            }

            a_evol->add_piece(a_fast);
            a_evol->add_piece(a_locked);

            e_evol->add_piece(zero_e);

            Lconv_evol->add_piece(Lconv_fast);
            Lconv_evol->add_piece(Lconv_locked);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            expected_real_quantities[SEMIMAJOR] = a_evol;
            expected_real_quantities[ECCENTRICITY] = e_evol;
            expected_real_quantities[CONV_INCLINATION] = zero_quantity;
            expected_real_quantities[RAD_INCLINATION] = zero_quantity;
            expected_real_quantities[CONV_PERIAPSIS] = zero_quantity;
            expected_real_quantities[RAD_PERIAPSIS] = zero_quantity;
            expected_real_quantities[CONV_ANGMOM] = Lconv_evol;
            expected_real_quantities[RAD_ANGMOM] = zero_quantity;

            return expected_real_quantities;
        }

    std::vector<const Core::OneArgumentDiffFunction *>
        test_OrbitSolver::calculate_expected_polar_1_0(double tdisk,
                                                       double wstar,
                                                       double worb)

        {
            double aorb = std::pow(
                         (
                             Core::AstroConst::G
                             *
                             (
                                 Core::AstroConst::solar_mass
                                 +
                                 Core::AstroConst::jupiter_mass
                             )
                         )
                         /
                         std::pow(
                             worb / Core::AstroConst::day,
                             2
                         ),
                         1.0 / 3.0
                     ) / Core::AstroConst::solar_radius;

            StellarEvolution::PolynomialEvolutionQuantity
                *disk_nan_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(Core::NaN, 1),
                    TSTART,
                    tdisk
                ),
                *fixed_a_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(aorb, 1),
                    tdisk,
                    MAX_AGE
                ),
                *fixed_e_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(),
                    tdisk,
                    MAX_AGE
                ),
                *disk_zero_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(),
                    TSTART,
                    tdisk
                ),
                *halfpi_evol = new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(M_PI / 2.0, 1),
                    tdisk,
                    MAX_AGE
                );
            __temp_functions.push_back(disk_nan_evol);
            __temp_functions.push_back(fixed_a_evol);
            __temp_functions.push_back(fixed_e_evol);
            __temp_functions.push_back(disk_zero_evol);
            __temp_functions.push_back(halfpi_evol);

            PiecewiseFunction *a_evol = new PiecewiseFunction,
                              *e_evol = new PiecewiseFunction,
                              *conv_incl_evol = new PiecewiseFunction,
                              *rad_incl_evol = new PiecewiseFunction;
            __temp_functions.push_back(a_evol);
            __temp_functions.push_back(e_evol);
            __temp_functions.push_back(conv_incl_evol);
            __temp_functions.push_back(rad_incl_evol);

            a_evol->add_piece(disk_nan_evol);
            a_evol->add_piece(fixed_a_evol);

            e_evol->add_piece(disk_nan_evol);
            e_evol->add_piece(fixed_e_evol);

            conv_incl_evol->add_piece(disk_zero_evol);
            conv_incl_evol->add_piece(halfpi_evol);

            rad_incl_evol->add_piece(disk_zero_evol);
            rad_incl_evol->add_piece(halfpi_evol);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            expected_real_quantities[SEMIMAJOR] = a_evol;
            expected_real_quantities[ECCENTRICITY] = e_evol;
            expected_real_quantities[CONV_INCLINATION] = conv_incl_evol;
            expected_real_quantities[RAD_INCLINATION] = rad_incl_evol;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] =
                new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(wstar, 1),
                    TSTART,
                    MAX_AGE
                );
            __temp_functions.push_back(
                expected_real_quantities[CONV_ANGMOM]
            );
            expected_real_quantities[RAD_ANGMOM] = &one_func;

            return expected_real_quantities;
        }

    std::vector<const Core::OneArgumentDiffFunction *>
        test_OrbitSolver::calculate_expected_polar_2_0(double tdisk,
                                                       double wstar,
                                                       double worb,
                                                       double phase_lag,
                                                       double &lconv_decay_rate,
                                                       double &semimajor)
        {
            semimajor = std::pow(
                    (
                        Core::AstroConst::G
                        *
                        (
                            Core::AstroConst::solar_mass
                            +
                            Core::AstroConst::jupiter_mass
                        )
                    )
                    /
                    std::pow(
                        worb / Core::AstroConst::day,
                        2
                    ),
                    1.0 / 3.0
                ) / Core::AstroConst::solar_radius;
            lconv_decay_rate = (
                0.3 * M_PI
                *
                (
                    Core::AstroConst::G
                    *
                    std::pow(
                        (
                            Core::AstroConst::jupiter_mass
                            /
                            std::pow(
                                semimajor * Core::AstroConst::solar_radius,
                                3
                            )
                        ),
                        2
                    )
                    *
                    std::pow(Core::AstroConst::solar_radius, 5)
                )
                /
                (
                    Core::AstroConst::solar_mass
                    *
                    std::pow(Core::AstroConst::solar_radius, 2)
                )
                *
                Core::AstroConst::day
                *
                Core::AstroConst::Gyr
                *
                phase_lag
            );
            double lorb = (
                (
                    Core::AstroConst::jupiter_mass
                    /
                    Core::AstroConst::solar_mass
                )
                /
                (
                    1.0
                    +
                    Core::AstroConst::jupiter_mass
                    /
                    Core::AstroConst::solar_mass
                )
                *
                semimajor * semimajor
                *
                worb
            );
            double lconv_offset = lconv_decay_rate * tdisk;

            std::valarray<double> spindown_coef(2);
            spindown_coef[0] = wstar + lconv_decay_rate * tdisk;
            spindown_coef[1] = - lconv_decay_rate;

            std::valarray<double> conv_cosinc_numer_coef(3),
                conv_cosinc_denom_coef(2),
                rad_cosinc_numer_coef(3),
                rad_cosinc_denom_coef(3);
            conv_cosinc_numer_coef[0] = - (
                (2.0 * wstar + lconv_offset) * lconv_offset
            );
            conv_cosinc_numer_coef[1] = 2.0 * wstar * lconv_decay_rate;
            conv_cosinc_numer_coef[2] = -std::pow(lconv_decay_rate, 2);

            conv_cosinc_denom_coef[0] = 2.0 * lorb * (wstar + lconv_offset);
            conv_cosinc_denom_coef[1] = -2.0 * lconv_decay_rate * lorb;

            conv_cosinc_numer_coef[0] += conv_cosinc_denom_coef[0];
            conv_cosinc_numer_coef[1] += conv_cosinc_denom_coef[1];

            rad_cosinc_numer_coef[0] = - 2.0 * lorb * lconv_decay_rate * tdisk;
            rad_cosinc_numer_coef[1] = 2.0 * lorb * lconv_decay_rate;

            rad_cosinc_denom_coef[0] = (2.0 * lorb * lorb
                                        -
                                        2.0 * wstar * lconv_decay_rate * tdisk
                                        -
                                        std::pow(lconv_decay_rate * tdisk, 2));
            rad_cosinc_denom_coef[1] =
                2.0 * lconv_decay_rate * (worb + lconv_decay_rate * tdisk);
            rad_cosinc_denom_coef[2] = - std::pow(lconv_decay_rate, 2);

            rad_cosinc_numer_coef += rad_cosinc_denom_coef;

            StellarEvolution::PolynomialEvolutionQuantity
                *disk_nan_evol =
                new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(Core::NaN, 1),
                    TSTART,
                    tdisk
                ),
                    *fixed_a_evol =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            std::valarray<double>(semimajor, 1),
                            tdisk,
                            MAX_AGE
                        ),
                    *fixed_e_evol
                        = new StellarEvolution::PolynomialEvolutionQuantity(
                            std::valarray<double>(),
                            tdisk,
                            MAX_AGE
                        ),
                    *disk_cosinc_evol =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            std::valarray<double>(2.0, 1),
                            TSTART,
                            tdisk
                        ),
                    *conv_cosinc_binary_numer =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            conv_cosinc_numer_coef,
                            tdisk, MAX_AGE
                        ),
                    *conv_cosinc_binary_denom =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            conv_cosinc_denom_coef,
                            tdisk,
                            MAX_AGE
                        ),
                    *rad_cosinc_binary_numer =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            rad_cosinc_numer_coef,
                            tdisk,
                            MAX_AGE
                        ),
                    *rad_cosinc_binary_denom =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            rad_cosinc_denom_coef,
                            tdisk,
                            MAX_AGE
                        ),
                    *lconv_disk =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            std::valarray<double>(wstar, 1),
                            TSTART,
                            tdisk
                        ),
                    *lconv_spindown =
                        new StellarEvolution::PolynomialEvolutionQuantity(
                            spindown_coef,
                            tdisk,
                            MAX_AGE
                        );
            __temp_functions.push_back(disk_nan_evol);
            __temp_functions.push_back(fixed_a_evol);
            __temp_functions.push_back(fixed_e_evol);
            __temp_functions.push_back(disk_cosinc_evol);
            __temp_functions.push_back(conv_cosinc_binary_numer);
            __temp_functions.push_back(conv_cosinc_binary_denom);
            __temp_functions.push_back(rad_cosinc_binary_numer);
            __temp_functions.push_back(rad_cosinc_binary_denom);
            __temp_functions.push_back(lconv_disk);
            __temp_functions.push_back(lconv_spindown);

            FunctionRatio *conv_cosinc_binary = new FunctionRatio(
                conv_cosinc_binary_numer,
                conv_cosinc_binary_denom
            );
            __temp_functions.push_back(conv_cosinc_binary);
            FunctionRatio *rad_cosinc_binary = new FunctionRatio(
                rad_cosinc_binary_numer,
                rad_cosinc_binary_denom
            );
            __temp_functions.push_back(rad_cosinc_binary);

            PiecewiseFunction *a_evol = new PiecewiseFunction,
                              *e_evol = new PiecewiseFunction,
                              *conv_cosinc_evol = new PiecewiseFunction,
                              *rad_cosinc_evol = new PiecewiseFunction,
                              *lconv_evol = new PiecewiseFunction;
            __temp_functions.push_back(a_evol);
            __temp_functions.push_back(e_evol);
            __temp_functions.push_back(conv_cosinc_evol);
            __temp_functions.push_back(rad_cosinc_evol);
            __temp_functions.push_back(lconv_evol);

            a_evol->add_piece(disk_nan_evol);
            a_evol->add_piece(fixed_a_evol);

            e_evol->add_piece(disk_nan_evol);
            e_evol->add_piece(fixed_e_evol);

            conv_cosinc_evol->add_piece(disk_cosinc_evol);
            conv_cosinc_evol->add_piece(conv_cosinc_binary);

            rad_cosinc_evol->add_piece(disk_cosinc_evol);
            rad_cosinc_evol->add_piece(rad_cosinc_binary);

            lconv_evol->add_piece(lconv_disk);
            lconv_evol->add_piece(lconv_spindown);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            expected_real_quantities[SEMIMAJOR] = a_evol;
            expected_real_quantities[ECCENTRICITY] = e_evol;
            expected_real_quantities[CONV_INCLINATION] = conv_cosinc_evol;
            expected_real_quantities[RAD_INCLINATION] = rad_cosinc_evol;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = lconv_evol;
            expected_real_quantities[RAD_ANGMOM] = &one_func;

            return expected_real_quantities;
        }

    std::vector<const Core::OneArgumentDiffFunction *>
        test_OrbitSolver::calculate_expected_oblique_m_0(
            unsigned m,
            double tdisk,
            double worb,
            double initial_inc,
            double initial_wstar,
            double phase_lag,
            double &min_wstar
        )
        {
            assert(m == 1 || m == 2);
            double aorb = std::pow(
                    (
                        Core::AstroConst::G
                        *
                        (
                            Core::AstroConst::solar_mass
                            +
                            Core::AstroConst::jupiter_mass
                        )
                    )
                    /
                    std::pow(
                        worb / Core::AstroConst::day,
                        2
                    ),
                    1.0 / 3.0
                ) / Core::AstroConst::solar_radius;
            double lorb = (Mjup_to_Msun
                           /
                           (1.0 + Mjup_to_Msun)
                           *
                           aorb * aorb
                           *
                           worb);

            double ltot = std::sqrt(
                std::pow(initial_wstar, 2)
                +
                2.0 * initial_wstar * lorb * std::cos(initial_inc)
                +
                std::pow(lorb, 2)
            );

            min_wstar = ltot - lorb;

            double linear_quantity_rate = (
                0.6 * M_PI
                *
                (
                    Core::AstroConst::G
                    *
                    std::pow(
                        (
                            Core::AstroConst::jupiter_mass
                            /
                            std::pow(
                                aorb * Core::AstroConst::solar_radius,
                                3
                            )
                        ),
                        2
                    )
                    *
                    std::pow(Core::AstroConst::solar_radius, 5)
                )
                /
                (
                    Core::AstroConst::solar_mass
                    *
                    std::pow(Core::AstroConst::solar_radius, 2)
                )
                *
                Core::AstroConst::day
                *
                Core::AstroConst::Gyr
                *
                phase_lag
            );

            StellarEvolution::PolynomialEvolutionQuantity *disk_nan_evol =
                new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(Core::NaN, 1),
                    TSTART,
                    tdisk
                );
            __temp_functions.push_back(disk_nan_evol);

            StellarEvolution::PolynomialEvolutionQuantity *fixed_a_evol =
                new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(aorb, 1),
                    tdisk,
                    MAX_AGE
                );
            __temp_functions.push_back(fixed_a_evol);

            StellarEvolution::PolynomialEvolutionQuantity *fixed_e_evol =
                new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(),
                    tdisk,
                    MAX_AGE
                );
            __temp_functions.push_back(fixed_e_evol);

            PiecewiseFunction *a_evol = new PiecewiseFunction,
                              *e_evol = new PiecewiseFunction;
            __temp_functions.push_back(a_evol);
            __temp_functions.push_back(e_evol);

            a_evol->add_piece(disk_nan_evol);
            a_evol->add_piece(fixed_a_evol);

            e_evol->add_piece(disk_nan_evol);
            e_evol->add_piece(fixed_e_evol);

            Core::OneArgumentDiffFunction *lconv_evol;
            if(m == 1)
               lconv_evol =
                   new InverseLinearLconvEvolution<Oblique10LinearQuantity>(
                       tdisk,
                       ltot,
                       lorb,
                       initial_wstar,
                       linear_quantity_rate
                   );
            else
                lconv_evol =
                    new InverseLinearLconvEvolution<Oblique20LinearQuantity>(
                        tdisk,
                        ltot,
                        lorb,
                        initial_wstar,
                        linear_quantity_rate / 2.0
                    );
            __temp_functions.push_back(lconv_evol);

            ConservedLEConvObliquityEvolution *conv_obliq_evol =
                new ConservedLEConvObliquityEvolution(*lconv_evol,
                                                      lorb,
                                                      ltot,
                                                      tdisk);
            __temp_functions.push_back(conv_obliq_evol);

            ConservedLERadObliquityEvolution *rad_obliq_evol =
                new ConservedLERadObliquityEvolution(*lconv_evol,
                                                     *conv_obliq_evol,
                                                     lorb,
                                                     tdisk);
            __temp_functions.push_back(rad_obliq_evol);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);
            expected_real_quantities[SEMIMAJOR] = a_evol;
            expected_real_quantities[ECCENTRICITY] = e_evol;
            expected_real_quantities[CONV_INCLINATION] = conv_obliq_evol;
            expected_real_quantities[RAD_INCLINATION] = rad_obliq_evol;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = lconv_evol;
            expected_real_quantities[RAD_ANGMOM] = &one_func;

            return expected_real_quantities;
        }

    void test_OrbitSolver::test_disk_locked_no_stellar_evolution()
    {
        try {
            StellarEvolution::MockStellarEvolution *no_evol =
                StellarEvolution::make_no_evolution();
            __star = make_const_lag_star(*no_evol,
                                         1.0,
                                         1.0,
                                         1.0);

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);

            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);
            expected_real_quantities[SEMIMAJOR] = &nan_func;
            expected_real_quantities[ECCENTRICITY] = &nan_func;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = &one_func;
            expected_real_quantities[RAD_ANGMOM] = &one_func;

            evolve(1.0, MAX_AGE, 1.0, &one);
            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          MAX_AGE);

            delete __solver;
            delete __system;

            evolve(1.0, MAX_AGE, 1.0, &zero);
            expected_real_quantities[RAD_ANGMOM] = new ExponentialPlusFunc(
                &one_func,
                -std::exp(TSTART/2.0),
                -0.5
            );
            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          MAX_AGE);

            delete expected_real_quantities[RAD_ANGMOM];
            delete no_evol;
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
                        ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
                        ex.what()).c_str());
        }
    }

    void test_OrbitSolver::test_disk_locked_with_stellar_evolution()
    {
        try {
            StellarEvolution::MockStellarEvolution *evol1 =
                StellarEvolution::make_linear_I_evolution();

            __star = make_const_lag_star(*evol1, 1.0, 1.0, 1.0);

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);

            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            std::valarray<double> temp_array=std::valarray<double>(1.0, 2);
            temp_array[0] = -1;
            StellarEvolution::PolynomialEvolutionQuantity *temp_poly =
                new StellarEvolution::PolynomialEvolutionQuantity(temp_array,
                                                                  TSTART,
                                                                  MAX_AGE);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);
            expected_real_quantities[SEMIMAJOR] = &nan_func;
            expected_real_quantities[ECCENTRICITY] = &nan_func;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] =
                new StellarEvolution::PolynomialEvolutionQuantity(
                    std::valarray<double>(1.0, 2),
                    TSTART,
                    MAX_AGE
                );
            expected_real_quantities[RAD_ANGMOM] =
                new ExponentialPlusFunc(temp_poly, 1, -0.5);

            double initial_Lrad = TSTART - 1 + std::exp(-TSTART / 2.0);
            evolve(1.0,
                   MAX_AGE,
                   1.0,
                   &initial_Lrad);
            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          MAX_AGE);

            delete expected_real_quantities[CONV_ANGMOM];
            delete expected_real_quantities[RAD_ANGMOM];
            delete temp_poly;
            delete __solver;
            delete __system;
            delete __star;
            delete evol1;

            std::valarray< std::valarray<double> > Ic_arr(
                std::valarray<double>(1.0, 1),
                2
            );
            Ic_arr[1][0]=-1.0/6.0;
            StellarEvolution::MockStellarEvolution evol2(
                -1,
                std::valarray< std::valarray<double> >(//R
                    std::valarray<double>(1.0, 1),
                    1
                ),
                Ic_arr,                                //Iconv
                std::valarray< std::valarray<double> >(//Irad
                    std::valarray<double>(1.0/6.0, 1),
                    2
                ),
                std::valarray< std::valarray<double> >(//Rcore
                    std::valarray<double>(0.5, 1),
                    1
                ),
                std::valarray< std::valarray<double> >(//Mcore
                    std::valarray<double>(1.0, 1),
                    2
                )
            );

            __star = make_const_lag_star(evol2, 1.0, 1.0, 1.0);

            std::valarray<double> Lc_coef(1.0, 2);
            Lc_coef[1]=-1.0/6.0;
            expected_real_quantities[SEMIMAJOR] = &nan_func;
            expected_real_quantities[ECCENTRICITY] = &nan_func;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] =
                new StellarEvolution::PolynomialEvolutionQuantity(Lc_coef,
                                                                  TSTART,
                                                                  MAX_AGE);
            expected_real_quantities[RAD_ANGMOM] =
                new StellarEvolution::PolynomialEvolutionQuantity(
                        std::valarray<double>(1.0/6.0, 2),
                        TSTART,
                        MAX_AGE
                );
            initial_Lrad = (TSTART / 2.0 + 1.0) / 6.0;
            evolve(1.0,
                   MAX_AGE,
                   1.0,
                   &initial_Lrad);
            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          MAX_AGE);

            delete expected_real_quantities[CONV_ANGMOM];
            delete expected_real_quantities[RAD_ANGMOM];
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()).c_str());
        }
    }

    void test_OrbitSolver::test_no_planet_evolution()
    {
        try {
            const double rt2 = std::sqrt(2.0);

            StellarEvolution::MockStellarEvolution
                *stellar_evol = StellarEvolution::make_no_evolution();
            double initial_Lstar[] = {0.0, 0.0};

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);
            expected_real_quantities[CONV_ANGMOM] = &zero_func;
            expected_real_quantities[RAD_ANGMOM] = &zero_func;
            ExpectedEvolutionMode<bool> unsat_wind_mode, sat_wind_mode;
            unsat_wind_mode.add_break(TSTART, false);
            sat_wind_mode.add_break(TSTART, true);

            test_no_planet_scenario(*stellar_evol,
                                    initial_Lstar,
                                    1.0,//Wind K
                                    2.0,//wind sat freq.
                                    Core::Inf,//core-env coupling timescale
                                    expected_real_quantities,
                                    unsat_wind_mode);

            delete stellar_evol;
            stellar_evol = StellarEvolution::make_linear_I_evolution();

            initial_Lstar[0] = 1.0;
            expected_real_quantities[CONV_ANGMOM] = &one_func;
            test_no_planet_scenario(*stellar_evol,
                                    initial_Lstar,
                                    0.0,//Wind K
                                    2.0,//wind sat freq.
                                    Core::Inf,//core-env coupling timescale
                                    expected_real_quantities,
                                    unsat_wind_mode);


            initial_Lstar[0] = 0.5 * (1.0 + std::exp(-TSTART));
            initial_Lstar[1] = 0.5 * (1.0 - std::exp(-TSTART));

            StellarEvolution::PolynomialEvolutionQuantity half_func(
                std::valarray<double>(0.5, 1),
                TSTART,
                MAX_AGE
            );
            expected_real_quantities[CONV_ANGMOM] = new ExponentialPlusFunc(
                &half_func,
                0.5,
                -1.0
            );
            expected_real_quantities[RAD_ANGMOM] = new ExponentialPlusFunc(
                &half_func,
                -0.5,
                -1.0
            );
            test_no_planet_scenario(*stellar_evol,
                                    initial_Lstar,
                                    0.0,//Wind K
                                    1.0,//wind sat freq.
                                    1.0,//core-env coupling timescale
                                    expected_real_quantities,
                                    unsat_wind_mode);

            delete expected_real_quantities[CONV_ANGMOM];
            delete expected_real_quantities[RAD_ANGMOM];

            initial_Lstar[0] = 1.0/(1.0+TSTART);
            initial_Lstar[1] = 1.0;

            double wind_sat_age = std::pow(2.0, 0.25) - 1;
            std::valarray<double> late_denom_coef(3);
            late_denom_coef[0] = 2.0 * rt2 - 2.0;
            late_denom_coef[1] = 4.0 * rt2;
            late_denom_coef[2] = 2.0 * rt2;
            StellarEvolution::PolynomialEvolutionQuantity
                one_func_early(std::valarray<double>(1.0, 1),
                               TSTART, wind_sat_age),
                one_plus_t(std::valarray<double>(1.0, 2), TSTART, MAX_AGE),
                late_denom2(late_denom_coef, wind_sat_age, MAX_AGE);
            FunctionToPower late_denom(&late_denom2, 0.5);
            FunctionRatio early_solution(&one_func_early, &one_plus_t),
                          late_solution(&one_plus_t, &late_denom);
            PiecewiseFunction full_solution;
            full_solution.add_piece(&early_solution);
            full_solution.add_piece(&late_solution);

            expected_real_quantities[CONV_ANGMOM] = &full_solution;
            expected_real_quantities[RAD_ANGMOM] = &one_func;
            ExpectedEvolutionMode<bool> changing_wind_mode;
            changing_wind_mode.add_break(TSTART, true);
            changing_wind_mode.add_break(wind_sat_age, false);

            test_no_planet_scenario(*stellar_evol,
                                    initial_Lstar,
                                    2.0,//Wind K
                                    1.0 / rt2, //wind sat freq.
                                    Core::Inf,//core-env coupling timescale
                                    expected_real_quantities,
                                    changing_wind_mode);

            double b1 = (std::sqrt(2) - 1) / (2.0 * rt2),
                   b2 = (std::sqrt(2) + 1) / (2.0 * rt2);
            initial_Lstar[0] = (b1 * std::exp(-TSTART / (2.0 + rt2))
                                +
                                b2 * std::exp(-TSTART / (2.0 - rt2)));
            initial_Lstar[1] = (
                0.5 / rt2 * std::exp(-1.0 / (2.0 + rt2) * (TSTART))
                -
                0.5 / rt2 * std::exp(-1.0 / (2.0 - rt2) * (TSTART))
            );

            ExponentialPlusFunc
                Lc1(&zero_func, b1, -1.0 / (2.0 + rt2)),
                Lc2(&zero_func, b2, -1.0 / (2.0 - rt2)),
                Lr1(&zero_func, 0.5 / rt2, -1.0 / (2.0 + rt2)),
                Lr2(&zero_func, -0.5 / rt2, -1.0 / (2.0 - rt2));

            expected_real_quantities[CONV_ANGMOM] = new FuncPlusFunc(&Lc1,
                                                                     &Lc2);
            expected_real_quantities[RAD_ANGMOM] = new FuncPlusFunc(&Lr1,
                                                                    &Lr2);
            delete stellar_evol;
            stellar_evol = StellarEvolution::make_no_evolution();
            test_no_planet_scenario(*stellar_evol,
                                    initial_Lstar,
                                    100.0,//Wind K
                                    0.1,//wind sat freq.
                                    1.0,//core-env coupling timescale
                                    expected_real_quantities,
                                    sat_wind_mode,
                                    2.0);
            delete expected_real_quantities[CONV_ANGMOM];
            delete expected_real_quantities[RAD_ANGMOM];

            delete stellar_evol;
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
                                    ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
                                    ex.what()).c_str());
        }
    }

    void test_OrbitSolver::test_unlocked_evolution()
    {
        try {
            ExpectedEvolutionMode<Core::EvolModeType> binary_mode;
            binary_mode.add_break(TSTART, Core::BINARY);
            ExpectedEvolutionMode<bool> unsat_wind_mode, sat_wind_mode;
            unsat_wind_mode.add_break(TSTART, false);
            sat_wind_mode.add_break(TSTART, true);

            StellarEvolution::MockStellarEvolution
                *no_evol = StellarEvolution::make_no_evolution(1.0);

            const double mplanet = 100;
            double lag = 1e-8 / mplanet;

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities = calculate_expected_unlocked_evolution(
                    lag,
                    mplanet
                );

            __star = make_const_lag_star(*no_evol,
                                         0.0,//wind K
                                         100.0,//wsat
                                         Core::Inf,//tcoup
                                         lag);

            while(true) {
                double initial_a = (
                    *expected_real_quantities[SEMIMAJOR]
                )(TSTART);

                std::valarray<double> initial_L(3);
                initial_L[0] = (*expected_real_quantities[CONV_ANGMOM])(
                    TSTART
                );
                initial_L[1] = 0.0;
                initial_L[2] = 0.0;

                evolve(0.0,//wdisk
                       TSTART,//tdisk
                       initial_a,
                       &(initial_L[0]),
                       0.0,//initial inclination
                       mplanet,//planet mass
                       Core::NaN,//tplanet
                       1.0,//stop evolution age
                       0.0001);//Rplanet
                test_solution(get_evolution(),
                              expected_real_quantities,
                              binary_mode,
                              unsat_wind_mode,
                              TSTART,
                              1.0);

                delete __system;
                delete __solver;

                expected_real_quantities =
                    calculate_expected_unlocked_evolution(
                        lag,
                        mplanet,
                        false
                    );

                initial_a = (*expected_real_quantities[SEMIMAJOR])(TSTART);
                initial_L[0] = (*expected_real_quantities[CONV_ANGMOM])(TSTART);

                evolve(0.0,//wdisk
                       TSTART,//tdisk
                       initial_a,
                       &(initial_L[0]),
                       0.0,//initial inclination
                       mplanet,//planet mass
                       Core::NaN,//tplanet
                       1.0,//stop evolution age
                       0.0001);//Rplanet
                test_solution(get_evolution(),
                              expected_real_quantities,
                              binary_mode,
                              sat_wind_mode,
                              TSTART,
                              1.0);

                if(__star) {
                    delete __star;
                    __star = NULL;

                    __primary_planet = new Planet::Planet(1.0, 1.0, 1.0);
                    __primary_planet->zone().setup(
                        std::vector<double>(),//Wtide breaks
                        std::vector<double>(),//W* breaks
                        std::vector<double>(1, 0.0),//Wtide pow.
                        std::vector<double>(1, 0.0),//W* pow.
                        lag
                    );
                } else
                    break;
            }

            delete no_evol;
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()).c_str());
        }
    }

    ///\brief The equation that should be solved in order to get the semimamjor
    ///axis at a gien time.
    ///
    ///The paramaters should be: alpha, beta, kappa, C, t
    double locked_unsat_eq(double a, void *params)
    {
        double *dbl_par = static_cast<double *>(params);
        double alpha = dbl_par[0],
               beta = dbl_par[1],
               kappa = dbl_par[2],
               c = dbl_par[3],
               t = dbl_par[4];
        return (
            0.1 * alpha * std::pow(a, 5)
            -
            0.5 * beta * std::pow(a, 3)
            +
            kappa*t
            -
            c
        );
    }

    ///\brief The equation that should be solved in order to get the semimamjor
    ///axis at a given time for the locked evolution when the wind is saturated.
    ///
    ///The paramaters should be: alpha, beta, kappa, C, t
    double locked_sat_eq(double a, void *params)
    {
        double *dbl_par = static_cast<double *>(params);
        double alpha = dbl_par[0],
               beta = dbl_par[1],
               kappa = dbl_par[2],
               c = dbl_par[3],
               t = dbl_par[4];
        return (
            alpha * a * a / 4.0
            -
            1.5 * beta * std::log(a)
            +
            kappa * t
            -
            c
        );
    }

    ///\brief The derivative of the equation that should be solved in order to
    ///get the semimamjor axis at a given time for the locked evolution when the
    ///wind is not saturated.
    ///
    ///The paramaters should be: alpha, beta.
    double locked_unsat_deriv(double a, void *params)
    {
        double *dbl_par = static_cast<double *>(params);
        double alpha = dbl_par[0],
               beta = dbl_par[1];
        return 0.5 * alpha * std::pow(a, 4) - 1.5 * beta * std::pow(a, 2);
    }

    ///\brief The derivative of the equation that should be solved in order to
    ///get the semimamjor axis at a gien time for the locked evolution when the
    ///wind is saturated.
    ///
    ///The paramaters should be: alpha, beta.
    double locked_sat_deriv(double a, void *params)
    {
        double *dbl_par = static_cast<double *>(params);
        double alpha = dbl_par[0],
               beta = dbl_par[1];
        return alpha * a / 2.0 - 1.5 * beta / a;
    }

    ///\brief The equation and its derivative that should be solved in order to
    ///get the semimamjor axis at a given time for the locked evolution when the
    ///wind is not saturated.
    ///
    ///The paramaters should be: alpha, beta, kappa, C, t
    void locked_unsat_eq_deriv(double a, void *params, double *f, double *df)
    {
        double *dbl_par = static_cast<double *>(params);
        double alpha = dbl_par[0],
               beta = dbl_par[1],
               kappa = dbl_par[2],
               c = dbl_par[3],
               t = dbl_par[4];
        *f = (
            0.1 * alpha * std::pow(a, 5)
            -
            0.5 * beta * std::pow(a, 3)
            +
            kappa * t
            -
            c
        );
        *df = (0.5 * alpha * std::pow(a, 4)
               -
               1.5 * beta * std::pow(a, 2));
    }

    ///\brief The equation and its derivative that should be solved in order to
    ///get the semimamjor axis at a given time for the locked evolution when the
    ///wind is saturated.
    ///
    ///The paramaters should be: alpha, beta, kappa, C, t
    void locked_sat_eq_deriv(double a, void *params, double *f, double *df)
    {
        double *dbl_par = static_cast<double *>(params);
        double alpha = dbl_par[0],
               beta = dbl_par[1],
               kappa = dbl_par[2],
               c = dbl_par[3],
               t = dbl_par[4];
        *f = (
            alpha * a * a / 4.0
            -
            1.5 * beta * std::log(a)
            +
            kappa * t
            -
            c
        );
        *df = alpha * a / 2.0 - 1.5 * beta / a;
    }

    void test_OrbitSolver::test_locked_evolution()
    {
        try {
            ExpectedEvolutionMode<Core::EvolModeType> expected_mode;
            expected_mode.add_break(TSTART, Core::BINARY);
            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);
            const double Ic = 0.001, Kwind = 1e-3, Kwind_s = 1;
            StellarEvolution::MockStellarEvolution *
                no_evol = StellarEvolution::make_no_evolution(1.0, Ic);

            double a1 = 3,
                   Lscale = (Core::AstroConst::jupiter_mass
                             /
                             std::pow(Core::AstroConst::solar_radius, 1.5)
                             *
                             std::sqrt(
                                 Core::AstroConst::G
                                 /
                                 (Core::AstroConst::jupiter_mass
                                  +
                                  Core::AstroConst::solar_mass)
                             )
                             *
                             Core::AstroConst::day),
                   beta = (
                       Ic * std::sqrt(
                           Core::AstroConst::G
                           *
                           (
                               Core::AstroConst::solar_mass
                               +
                               Core::AstroConst::jupiter_mass
                           )
                       )
                       *
                       Core::AstroConst::day
                       /
                       std::pow(Core::AstroConst::solar_radius, 1.5)
                   ),
                   wsat_s = 0.1, //Must be adjusted if a1 is adjusted
//                   wsat = Core::Inf,
                   kappa=(
                       Kwind
                       *
                       std::pow(
                           (
                               Core::AstroConst::G
                               *
                               (
                                   Core::AstroConst::solar_mass
                                   +
                                   Core::AstroConst::jupiter_mass
                               )
                               /
                               std::pow(Core::AstroConst::solar_radius, 3)
                           ),
                           1.5
                       )
                       *
                       std::pow(Core::AstroConst::day, 3)
                   ),
                   kappa_s = (
                       Kwind_s * wsat_s * wsat_s *
                       std::sqrt(
                           Core::AstroConst::G
                           *
                           (
                               Core::AstroConst::solar_mass
                               +
                               Core::AstroConst::jupiter_mass
                           )
                           /
                           std::pow(Core::AstroConst::solar_radius, 3)
                       )
                       *
                       Core::AstroConst::day
                   ),
                   int_const = (Lscale / 10.0 * std::pow(a1, 5)
                                -
                                beta / 2.0 * std::pow(a1, 3)
                                +
                                kappa),
                   int_const_s = (Lscale / 4.0 * std::pow(a1, 2)
                                  -
                                  3.0 * beta / 2.0 * std::log(a1)
                                  +
                                  kappa_s);
/*            double solver_params[] = {Lscale, beta, kappa, int_const, 0.0};
            double a0 = solve(a1,
                              0.0,
                              1e-9,
                              &locked_unsat_eq,
                              &locked_unsat_deriv,
                              &locked_unsat_eq_deriv,
                              static_cast<void*>(solver_params));
            solver_params[4] = TSTART;
            double astart = solve(a0,
                                  0.0,
                                  1e-9,
                                  &locked_unsat_eq,
                                  &locked_unsat_deriv,
                                  &locked_unsat_eq_deriv,
                                  static_cast<void*>(solver_params));
            solver_params[2] = kappa_s;
            solver_params[3] = int_const_s;
            solver_params[4] = 0;
            double a0_s = solve(a1,
                                0.0,
                                1e-9,
                                &locked_sat_eq,
                                &locked_sat_deriv,
                                &locked_sat_eq_deriv,
                                static_cast<void*>(solver_params));
            solver_params[4] = TSTART;
            double astart_s = solve(a1,
                                    0.0,
                                    1e-9,
                                    &locked_sat_eq,
                                    &locked_sat_deriv,
                                    &locked_sat_eq_deriv,
                                    static_cast<void*>(solver_params));
            solver_params[4] = 1.0;

    		double w0=std::sqrt(
                           AstroConst::G*
                           (AstroConst::solar_mass+AstroConst::jupiter_mass)/
                           std::pow(a0*AstroConst::solar_radius, 3))*
                       AstroConst::day,
                   w0_s=std::sqrt(
                           AstroConst::G*
                           (AstroConst::solar_mass+AstroConst::jupiter_mass)/
                           std::pow(a0_s*AstroConst::solar_radius, 3))*
                       AstroConst::day,
                   w1=std::sqrt(
                           AstroConst::G*
                           (AstroConst::solar_mass+AstroConst::jupiter_mass)/
                           std::pow(a1_check*AstroConst::solar_radius, 3))*
                       AstroConst::day,
                   w1_s=std::sqrt(
                           AstroConst::G*
                           (AstroConst::solar_mass+AstroConst::jupiter_mass)/
                           std::pow(a1_check_s*AstroConst::solar_radius, 3))*
                       AstroConst::day;*/
            std::valarray<double> a_transform_coef(0.0, 6),
                                  t_coef(2),
                                  t_coef_s(2);
            a_transform_coef[5] = Lscale / 10.0;
            a_transform_coef[3] = -beta / 2.0;
            t_coef[0] = int_const;
            t_coef[1] = -kappa;
            t_coef_s[0] = int_const_s;
            t_coef_s[1] = -kappa_s;
            std::valarray<double> Lconv_term1_coef(0.0, 2),
                                  Lconv_term2_coef(0.0, 3),
                                  identity_coef(0.0, 2),
                                  a_term1_coef_s(0.0, 3),
                                  Lconv_term1_coef_s(0.0, 2),
                                  Lc_beta_coef_s(0.0, 2);
            Lconv_term1_coef[1] = 1.0 / beta * std::pow(10.0 / Lscale,
                                                        3.0 / 10.0);
            Lconv_term2_coef[2] = -2.0 / std::pow(beta, 3);
            identity_coef[1] = 1;
            a_term1_coef_s[2] = Lscale / 4.0;
            Lconv_term1_coef_s[1] = 1.0 / beta * std::pow(4.0 / Lscale,
                                                          3.0 / 4.0);
            Lc_beta_coef_s[1] = 1.0 / beta;

            StellarEvolution::PolynomialEvolutionQuantity
                a_transform(a_transform_coef, 0.0, Core::Inf),
                a_transform1_s(a_term1_coef_s, 0.0, Core::Inf),
                transformed_a_evol(t_coef, TSTART, 1.0),
                transformed_a_evol_s(t_coef_s, TSTART, 1.0),
                Lconv_term1_poly(Lconv_term1_coef, 0.0, Core::Inf),
                Lconv_term2_poly(Lconv_term2_coef, 0.0, Core::Inf),
                Lconv_term1_poly_s(Lconv_term1_coef_s, 0.0, Core::Inf),
                Lc_beta_s(Lc_beta_coef_s, 0.0, Core::Inf),
                identity(identity_coef, 0.0, Core::Inf);
            FunctionToPower L_transform1(&Lconv_term1_poly, -10.0 / 3.0),
                            L_transform2(&Lconv_term2_poly, -1.0),
                            L_transform1_s(&Lconv_term1_poly_s, -4.0 / 3.0);
            LogFunction log_a(&identity),
                        log_Lc_beta_s(&Lc_beta_s);
            ScaledFunction a_transform2_s(&log_a, -1.5 * beta),
                           L_transform2_s(&log_Lc_beta_s, beta);
            FuncPlusFunc L_transform(&L_transform1, &L_transform2),
                         a_transform_s(&a_transform1_s, &a_transform2_s),
                         L_transform_s(&L_transform1_s, &L_transform2_s);
#if 0
            std::valarray<double> initial_orbit(2);

            initial_orbit[0]=astart;
            initial_orbit[1]=0.0;
            std::cout.precision(16);
            std::cout.setf(std::ios_base::scientific);
            Star star_not_saturated_wind_no_coupling(1.0, 0.0, Kwind, wsat, Inf,
                    0.0, 0.0, 0.0, no_evol);
            Planet planet1(&star_not_saturated_wind_no_coupling, 1.0, 1.0, 1.0);
            StellarSystem system1(&star_not_saturated_wind_no_coupling, &planet1);
            OrbitSolver solver(tstart, 1.0, 1e-9);
            solver(system1, Inf, 0.0, a0, tstart, LOCKED_TO_PLANET, initial_orbit, true);
            TransformedSolution to_check(a_transform, L_transform, identity, -Inf);
            to_check(solver);
            test_solution(to_check, transformed_a_evol, transformed_a_evol,
                    zero_func, tstart, 1.0, expected_mode);

            initial_orbit[0]=astart_s;
            initial_orbit[1]=0.0;
            Star star_saturated_wind_no_coupling(1.0, 0.0, Kwind_s, wsat_s, Inf,
                    0.0, 0.0, 0.0, no_evol);
            Planet planet2(&star_saturated_wind_no_coupling, 1.0, 1.0, 1.0);
            StellarSystem system2(&star_saturated_wind_no_coupling, &planet2);
            solver(system2, Inf, 0.0, a0_s, tstart, LOCKED_TO_PLANET,
                    initial_orbit, true);
            TransformedSolution to_check_s(a_transform_s, L_transform_s,
                    identity, -Inf);
            to_check_s(solver);
            test_solution(to_check_s, transformed_a_evol_s, transformed_a_evol_s,
                    zero_func, tstart, 1.0, expected_mode);
#endif
            delete no_evol;
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()).c_str());
        }
    }

    void test_OrbitSolver::test_disklocked_to_locked_to_noplanet()
    {
        try {
            const double Ic = 0.001,
                         Kwind = 2e-4,
                         tdisk = 1,
                         tfinal = 2;
            StellarEvolution::MockStellarEvolution *
                no_evol = StellarEvolution::make_no_evolution(1.0, Ic);

            double afinal = 1.75,
                   Lscale = (Core::AstroConst::jupiter_mass
                             /
                             std::pow(Core::AstroConst::solar_radius, 1.5)
                             *
                             std::sqrt(
                                 Core::AstroConst::G
                                 /
                                 (
                                     Core::AstroConst::jupiter_mass
                                     +
                                     Core::AstroConst::solar_mass
                                 )
                             )
                             *
                             Core::AstroConst::day),
                   beta = (
                       Ic
                       *
                       std::sqrt(
                           Core::AstroConst::G
                           *
                           (
                               Core::AstroConst::solar_mass
                               +
                               Core::AstroConst::jupiter_mass
                           )
                       )
                       *
                       Core::AstroConst::day
                       /
                       std::pow(Core::AstroConst::solar_radius, 1.5)
                   ),
                   wsat = 1e100,
                   kappa = (
                       Kwind
                       *
                       std::pow(
                           Core::AstroConst::G
                           *
                           (
                               Core::AstroConst::solar_mass
                               +
                               Core::AstroConst::jupiter_mass
                           )
                           /
                           std::pow(Core::AstroConst::solar_radius, 3),
                           1.5
                       )
                       *
                       std::pow(Core::AstroConst::day, 3)
                   ),
                   int_const = (Lscale / 10.0 * std::pow(afinal, 5)
                                -
                                beta / 2.0 * std::pow(afinal, 3)
                                +
                                kappa * tfinal);
            double solver_params[]={Lscale, beta, kappa, int_const, tdisk};
            double ainitial = solve(2.0 * afinal,
                                    0.0,
                                    1e-9,
                                    &locked_unsat_eq,
                                    &locked_unsat_deriv,
                                    &locked_unsat_eq_deriv,
                                    static_cast<void*>(solver_params)),
                   wdisk = (
                       std::sqrt(
                           Core::AstroConst::G
                           *
                           (
                               Core::AstroConst::solar_mass
                               +
                               Core::AstroConst::jupiter_mass
                           )
                           /
                           std::pow(ainitial * Core::AstroConst::solar_radius, 3)
                       )
                       *
                       Core::AstroConst::day
                   );

            __star = make_const_lag_star(
                *no_evol,
                Kwind,
                wsat,
                Core::Inf,//core-env coupling timescale
                lag_from_lgQ(0.1, (Core::AstroConst::jupiter_mass
                                   /
                                   Core::AstroConst::solar_mass))
            );
            Planet::Planet planet(Mjup_to_Msun, Rjup_to_Rsun);
            planet.configure(true, //init
                             tdisk, //age
                             __star->mass(), //mass
                             ainitial * (1.0 - 1e-14),//planet formation semimajor
                             0.0, //eccentricity
                             &zero, //spin angmom
                             NULL, //inclination
                             NULL, //periapsis
                             false, //locked surface
                             true, //zero outer inclination
                             true);//zero outer periapsis

            __system = new Evolve::DiskBinarySystem(
                *__star,
                planet,
                ainitial * (1.0 - 1e-14), //semimajor
                0.0, //eccentricity
                0, //inclination
                wdisk, //Wdisk
                tdisk, //disk dissipation age
                tdisk //planet formation age
            );

            __system->configure(true, //init
                                TSTART,
                                Core::NaN, //semimajor
                                Core::NaN, //eccentricity
                                &zero, //spin angmom
                                NULL, //inclination
                                NULL, //periapsis
                                Core::LOCKED_SURFACE_SPIN);


            double adestr = __system->minimum_semimajor(),
                   tdestr = ((beta / 2.0 * std::pow(adestr, 3)
                              +
                              int_const
                              -
                              Lscale / 10.0 * std::pow(adestr, 5)) / kappa),
                   Lc_at_destr = (beta / std::pow(adestr, 1.5)
                                  +
                                  Lscale * std::sqrt(adestr));


            std::valarray<double> a_transform_coef(0.0, 6), t_coef(2);
            a_transform_coef[5] = Lscale / 10.0;
            a_transform_coef[3] = -beta / 2.0;
            t_coef[0] = int_const;
            t_coef[1] = -kappa;
            std::valarray<double> Lconv_term1_coef(0.0, 2),
                                  Lconv_term2_coef(0.0, 3),
                                  identity_coef(0.0, 2),
                                  noplanet_Lconv_m2_coef(2);
            Lconv_term1_coef[1] = 1.0 / beta * std::pow(10.0 / Lscale, 0.3);
            Lconv_term2_coef[2] = -2.0 / std::pow(beta, 3);
            noplanet_Lconv_m2_coef[0] = (
                std::pow(Lc_at_destr, -2)
                -
                2.0 * Kwind * tdestr / std::pow(Ic, 3)
            );
            noplanet_Lconv_m2_coef[1] = 2.0 * Kwind / std::pow(Ic, 3);
            identity_coef[1] = 1;
            StellarEvolution::PolynomialEvolutionQuantity
                identity(identity_coef, -Core::Inf, Core::Inf),
                disk_nan_evol(std::valarray<double>(Core::NaN, 2),
                              TSTART,
                              tdisk),
                locked_a_transform(a_transform_coef, -Core::Inf, Core::Inf),
                locked_aLconv_evol(t_coef, tdisk, tdestr),
                locked_e_evol(std::valarray<double>(), tdisk, tdestr),
                noplanet_nan_evol(std::valarray<double>(Core::NaN, 2),
                                  tdestr,
                                  tfinal),
                disk_Lconv_evol(std::valarray<double>(wdisk * Ic, 1),
                                TSTART,
                                tdisk),
                noplanet_Lconv_m2_evol(noplanet_Lconv_m2_coef, tdestr, tfinal),
                Lconv_term1_poly(Lconv_term1_coef, -Core::Inf, Core::Inf),
                Lconv_term2_poly(Lconv_term2_coef, -Core::Inf, Core::Inf),
                Lrad_evol(std::valarray<double>(), TSTART, tfinal);
            FunctionToPower L_transform1(&Lconv_term1_poly, -10.0 / 3.0),
                            L_transform2(&Lconv_term2_poly, -1.0),
                            noplanet_Lconv_evol(&noplanet_Lconv_m2_evol, -0.5);
            LogFunction log_a(&identity);
            FuncPlusFunc locked_Lconv_transform(&L_transform1, &L_transform2);
            PiecewiseFunction a_evol, e_evol, Lconv_evol;
            a_evol.add_piece(&disk_nan_evol);
            a_evol.add_piece(&locked_aLconv_evol);
            a_evol.add_piece(&noplanet_nan_evol);
            e_evol.add_piece(&disk_nan_evol);
            e_evol.add_piece(&locked_e_evol);
            e_evol.add_piece(&noplanet_nan_evol);
            Lconv_evol.add_piece(&disk_Lconv_evol);
            Lconv_evol.add_piece(&locked_aLconv_evol);
            Lconv_evol.add_piece(&noplanet_Lconv_evol);

            std::vector< const Core::OneArgumentDiffFunction * >
                transformations(NUM_REAL_QUANTITIES - 1, &identity);
            TransformedSolution to_check(transformations, TSTART);
            transformations[SEMIMAJOR] = &locked_a_transform;
            transformations[CONV_ANGMOM] = &locked_Lconv_transform;
            to_check.add_transformation(transformations, tdisk);
            transformations[SEMIMAJOR] = &identity;
            transformations[CONV_ANGMOM] = &identity;
            to_check.add_transformation(transformations, tdestr);

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
            expected_evol_mode.add_break(tdisk, Core::BINARY);
            expected_evol_mode.add_break(tdestr, Core::SINGLE);

            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            __star->detect_saturation();
            __solver = new Evolve::OrbitSolver(tfinal, 1e-8);
            (*__solver)(*__system,
                        (tfinal - __system->age()) / 10000.0, //time step
                        std::list<double>()); //no required ages*/

            std::vector< const std::list<double> * >
                tabulated_evolution = get_evolution();
            const std::vector< const std::list<double> * > &
                transformed_evolution = to_check(tabulated_evolution);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);
            expected_real_quantities[SEMIMAJOR] = &a_evol;
            expected_real_quantities[ECCENTRICITY] = &e_evol;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = &Lconv_evol;
            expected_real_quantities[RAD_ANGMOM] = &Lrad_evol;

            test_solution(transformed_evolution,
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          tfinal);
            delete no_evol;

#if 0
            Star star_not_saturated_wind_no_coupling(
                1.0,//M*
                1.0e-10,//Q*
                Kwind,
                wsat,
                Inf,//core-env coupling timescale
                0.0,//Q transition width
                wdisk,
                tdisk,
                no_evol);
            Planet planet1(&star_not_saturated_wind_no_coupling, 1.0, 1.0, 1.0);
            StellarSystem system1(&star_not_saturated_wind_no_coupling, &planet1);
    /*		std::cout << std::endl << "alpha=" << Lscale
                << std::endl << "beta=" << beta
                << std::endl << "kappa=" << kappa
                << std::endl << "wdisk=" << wdisk << std::endl;
            std::cout << std::endl << "ainitial=" << ainitial << std::endl;
            std::cout << std::endl << "Destruction a=" << adestr
                << ", t=" << tdestr << ", " << ", Lc=" << Lc_at_destr
                << ", no planet time="
                << (Lscale*(std::pow(adestr, 5)-std::pow(afinal, 5))/10.0 -
                        beta*(std::pow(adestr, 3)-std::pow(afinal, 3))/2.0)/kappa
                << std::endl;*/

            OrbitSolver solver(tstart, tfinal, 1e-8);
            solver(system1,
                   Inf,//Max step
                   0.0,//Planet formation age
                   ainitial/AU_Rsun,//planet formation semimajor
                   tstart);//Start age

            test_solution(to_check,
                          a_evol,
                          Lconv_evol,
                          Lrad_evol,
                          tstart,
                          tfinal,
                          expected_mode);
#endif
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")
                                    +
                                    ex.what()).c_str());
        }

    }

    void test_OrbitSolver::test_disklocked_to_fast_to_noplanet()
    {
        try {
            const double TDISK = 1,
                         TDESTR = 2,
                         WDISK = 0.0,
                         I_CONV = 1,
                         L_SCALE = (
                             -Core::AstroConst::jupiter_mass
                             /
                             std::pow(Core::AstroConst::solar_radius, 1.5)
                             *
                             std::sqrt(
                                 Core::AstroConst::G
                                 /
                                 (
                                     Core::AstroConst::jupiter_mass
                                     +
                                     Core::AstroConst::solar_mass
                                 )
                             )
                             *
                             Core::AstroConst::day
                         );

            StellarEvolution::MockStellarEvolution *
                no_evol = StellarEvolution::make_no_evolution(1.0, I_CONV);

            double lgQ = 8,
                   alpha = (
                       -4.5
                       *
                       std::sqrt(
                           Core::AstroConst::G
                           /
                           (
                               Core::AstroConst::solar_radius
                               *
                               Core::AstroConst::solar_mass
                           )
                       )
                       *
                       Core::AstroConst::jupiter_mass / std::pow(10.0, lgQ)
                       *
                       Core::AstroConst::Gyr / Core::AstroConst::solar_radius
                   );

            __star = make_const_lag_star(
                *no_evol,
                0.0,
                1.0,
                Core::Inf,
                lag_from_lgQ(lgQ,
                             (Core::AstroConst::jupiter_mass
                              /
                              Core::AstroConst::solar_mass))
            );

            double adestr = (
                2.44
                *
                std::pow(
                    (
                        __star->mass() * Core::AstroConst::solar_mass
                        /
                        Core::AstroConst::jupiter_mass
                    ),
                    1.0 / 3.0
                )
                *
                Core::AstroConst::jupiter_radius
                /
                Core::AstroConst::solar_radius
            );
            double a6p5_offset = (std::pow(adestr, 6.5)
                                  -
                                  6.5 * alpha * TDESTR),
                   a_formation = std::pow(a6p5_offset + 6.5 * alpha * TDISK,
                                          1.0 / 6.5);

            evolve(WDISK,
                   TDISK,
                   a_formation,
                   &zero);//Initial L*

            ExpectedEvolutionMode<Core::EvolModeType> expected_mode;
            expected_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
            expected_mode.add_break(TDISK, Core::BINARY);
            expected_mode.add_break(TDESTR, Core::SINGLE);

            ExpectedEvolutionMode<bool> unsat_wind_mode;
            unsat_wind_mode.add_break(TSTART, false);

            std::valarray<double> a6p5_poly_coef(2);
            a6p5_poly_coef[0] = a6p5_offset;
            a6p5_poly_coef[1] = 6.5 * alpha;
            StellarEvolution::PolynomialEvolutionQuantity
                a6p5_evol(a6p5_poly_coef, TDISK, TDESTR),
                Lconv_disk(std::valarray<double>(I_CONV*WDISK, 1),
                           TSTART,
                           TDISK),
                Lconv_noplanet(
                    std::valarray<double>(
                        -L_SCALE * std::sqrt(a_formation) + I_CONV * WDISK,
                        1
                    ),
                    TDESTR,
                    MAX_AGE
                ),
                nan_disk(std::valarray<double>(Core::NaN, 2), TSTART, TDISK),
                nan_noplanet(std::valarray<double>(Core::NaN, 2),
                             TDESTR,
                             MAX_AGE),
                e_fast(std::valarray<double>(0.0, 1),
                       TDISK,
                       TDESTR),
                Lrad_evol(std::valarray<double>(), TSTART, MAX_AGE);
            FunctionToPower a_fast(&a6p5_evol, 1.0 / 6.5),
                            sqrta_evol(&a6p5_evol, 1.0 / 13.0);
            ExponentialPlusFunc Lconv_unscaled(
                &sqrta_evol,
                I_CONV * WDISK / L_SCALE - std::sqrt(a_formation),
                0
            );
            ScaledFunction Lconv_fast(&Lconv_unscaled, L_SCALE);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            PiecewiseFunction a_evol, e_evol, conv_angmom_evol;

            a_evol.add_piece(&nan_disk);
            a_evol.add_piece(&a_fast);
            a_evol.add_piece(&nan_noplanet);

            e_evol.add_piece(&nan_disk);
            e_evol.add_piece(&e_fast);
            e_evol.add_piece(&nan_noplanet);

            conv_angmom_evol.add_piece(&Lconv_disk);
            conv_angmom_evol.add_piece(&Lconv_fast);
            conv_angmom_evol.add_piece(
                &Lconv_noplanet
            );

            expected_real_quantities[SEMIMAJOR] = &a_evol;
            expected_real_quantities[ECCENTRICITY] = &e_evol;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_ANGMOM] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = &conv_angmom_evol;
            expected_real_quantities[RAD_ANGMOM] = &zero_func;

            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_mode,
                          unsat_wind_mode,
                          TSTART,
                          MAX_AGE);

            delete no_evol;
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
                    ex.what()+": "+ex.get_message()).c_str());
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
                    ex.what()).c_str());
        }
    }

    void test_OrbitSolver::test_disklocked_to_fast_to_locked()
    {
        try {
            const double lgQ = 8,
                         tdisk = 1,
                         async = 2.5,
                         tsync = 2.0,
                         tend = 3;

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities
                =
                calculate_expected_disklocked_to_fast_to_locked(
                    lgQ,
                    tdisk,
                    async,
                    tsync,
                    tend
                );

            double wsync = Core::orbital_angular_velocity(1.0,
                                                          Mjup_to_Msun,
                                                          async),
                   Lconv_sync = (*expected_real_quantities[CONV_ANGMOM])(
                       (tsync + tend) / 2.0
                   ),
                   Iconv = Lconv_sync / wsync,
                   Ldisk = (*expected_real_quantities[CONV_ANGMOM])(
                       (TSTART + tdisk) / 2.0
                   ),
                   wdisk = Ldisk / Iconv,
                   a_formation = (*expected_real_quantities[SEMIMAJOR])(tdisk),
                   phase_lag = lag_from_lgQ(lgQ,
                                            (Core::AstroConst::jupiter_mass
                                             /
                                             Core::AstroConst::solar_mass));

            StellarEvolution::MockStellarEvolution *
                no_evol = StellarEvolution::make_no_evolution(1.0, Iconv);

            __star = make_const_lag_star(
                *no_evol,//evolution
                0.0,//Kwind
                100.0,//wsat
                Core::Inf,//tcoup
                phase_lag
            );

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
            expected_evol_mode.add_break(tdisk, Core::BINARY);

            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            while(true) {
                evolve(wdisk,//wdisk,
                       tdisk,//tdisk
                       a_formation,//initial semimajor
                       (__star ? &zero : &Ldisk),//initial L*
                       0.0,//initial inclination
                       1.0,//planet_mass
                       Core::NaN,//time of secondary formation
                       tend,//max evolution age
                       0.01,//planet radius
                       1e-7);//precision

                test_solution(get_evolution(),
                              expected_real_quantities,
                              expected_evol_mode,
                              expected_wind_mode,
                              (__star ? TSTART : tdisk),
                              tend);
                if(__star) {
                    delete __star;
                    __star = NULL;

                    __primary_planet = new Planet::Planet(1.0, 1.0, Iconv);
                    __primary_planet->zone().setup(
                        std::vector<double>(),//Wtide breaks
                        std::vector<double>(),//W* breaks
                        std::vector<double>(1, 0.0),//Wtide pow.
                        std::vector<double>(1, 0.0),//W* pow.
                        phase_lag
                    );
                } else
                    break;
            }
            delete no_evol;

        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(
                false,
                (
                    std::string("Unexpected exception thrown: ")
                    +
                    ex.what() + ": " + ex.get_message()).c_str()
            );
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(
                false,
                (
                    std::string("Unexpected exception thrown: ")
                    +
                    ex.what()
                ).c_str()
            );
        }
    }

    void test_OrbitSolver::test_disklocked_to_locked_to_fast()
    {
        try {
            const double
                a0=3.2,
                abreak=3.0,
                adeath=2.0,
                tdisk=1,
                tbreak=2,
                tdeath=3,
                tend=4,
                Rp = (
                    adeath
                    *
                    Core::AstroConst::solar_radius
                    /
                    Core::AstroConst::jupiter_radius
                    /
                    (
                        2.44 * std::pow(Core::AstroConst::solar_mass
                                        /
                                        Core::AstroConst::jupiter_mass,
                                        1.0 / 3.0)
                    )
                ),
                beta = (
                    std::sqrt(
                        Core::AstroConst::G
                        *
                        (
                            Core::AstroConst::solar_mass
                            +
                            Core::AstroConst::jupiter_mass
                        )
                    )
                    *
                    Core::AstroConst::day
                    /
                    std::pow(Core::AstroConst::solar_radius, 1.5)
                ),
                wdisk = beta / std::pow(a0, 1.5),
                gamma = (
                    (std::pow(abreak, 6.5) - std::pow(adeath, 6.5))
                    /
                    (tdeath - tbreak)
                ),
                Q = (
                    9.0 * 13.0 / 4.0
                    *
                    std::sqrt(Core::AstroConst::G
                              /
                              Core::AstroConst::solar_mass)
                    *
                    Core::AstroConst::jupiter_mass
                    *
                    Core::AstroConst::Gyr
                    /
                    (gamma * std::pow(Core::AstroConst::solar_radius, 1.5))
                ),
                alphaL = (
                    Core::AstroConst::jupiter_mass
                    /
                    std::pow(Core::AstroConst::solar_radius, 1.5)
                    *
                    std::sqrt(Core::AstroConst::G
                              /
                              (
                                  Core::AstroConst::jupiter_mass
                                  +
                                  Core::AstroConst::solar_mass
                              )
                    )
                    *
                    Core::AstroConst::day
                ),
                Ic = (
                    (
                        (
                            (std::pow(a0, 5) - std::pow(abreak, 5))
                            *
                            std::pow(abreak, 3.5)
                            *
                            6.5
                        )
                        -
                        5.0 * gamma * abreak * abreak
                    )
                    /
                    (
                        (
                            (std::pow(a0, 3) - std::pow(abreak, 3))
                            *
                            std::pow(abreak, 3.5) * 6.5
                        )
                        -
                        3.0 * gamma
                    )
                    *
                    alphaL
                    /
                    (5.0 * beta)
                ),
                kappa = (
                    (std::pow(a0, 5) - std::pow(abreak, 5)) * alphaL / 10.0
                    -
                    (std::pow(a0, 3) - std::pow(abreak, 3)) * beta * Ic / 2.0
                ),
                Kwind = (
                    kappa
                    /
                    std::pow(
                        Core::AstroConst::G
                        *
                        (
                            Core::AstroConst::solar_mass
                            +
                            Core::AstroConst::jupiter_mass
                        )
                        /
                        std::pow(Core::AstroConst::solar_radius, 3),
                        1.5
                    )
                    /
                    std::pow(Core::AstroConst::day, 3)
                ),
                locked_a_int_const = (alphaL * std::pow(a0, 5) / 10.0
                                      -
                                      beta * Ic * std::pow(a0, 3) / 2.0
                                      +
                                      kappa * tdisk);

            StellarEvolution::MockStellarEvolution *
                no_evol = StellarEvolution::make_no_evolution(1.0, Ic);

            __star = make_const_lag_star(
                *no_evol,//evolution
                Kwind,//Kwind
                100.0,//wsat
                Core::Inf,//tcoup
                lag_from_lgQ(
                    std::log10(Q),
                    (Core::AstroConst::jupiter_mass
                     /
                     Core::AstroConst::solar_mass)
                )
            );
            evolve(wdisk,//wdisk,
                   tdisk,//tdisk
                   a0 * (1.0 + 1e-14),//initial semimajor
                   &zero,//initial L*
                   0.0,//initial inclination
                   1.0,//planet_mass
                   Core::NaN,//form the planet when disk dissipates
                   tend,//max evolution age
                   Rp);//planet radius

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
            expected_evol_mode.add_break(tdisk, Core::BINARY);
            expected_evol_mode.add_break(tdeath, Core::SINGLE);

            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            std::valarray<double> a_locked_transform_coef(0.0, 6),
                                  a_locked_evol_coef(2),
                                  identity_coef(0.0, 2),
                                  a6p5_fast_evol_coef(2),
                                  Lconv_locked_term1_coef(0.0, 2),
                                  Lconv_locked_term2_coef(0.0, 3);
            a_locked_transform_coef[5] = alphaL / 10.0;
            a_locked_transform_coef[3] = -beta * Ic / 2.0;
            a_locked_evol_coef[0] = locked_a_int_const;
            a_locked_evol_coef[1] = -kappa;
            identity_coef[1] = 1.0;
            a6p5_fast_evol_coef[0] = gamma * tbreak + std::pow(abreak, 6.5);
            a6p5_fast_evol_coef[1] = -gamma;
            Lconv_locked_term1_coef[1] = (1.0
                                          /
                                          (beta * Ic)
                                          *
                                          std::pow(10.0 / alphaL, 3.0 / 10.0));
            Lconv_locked_term2_coef[2] = -2.0 / std::pow(beta * Ic, 3);

            StellarEvolution::PolynomialEvolutionQuantity
                identity(identity_coef, -Core::Inf, Core::Inf),
                a_disk_transform = identity,
                nan_disk_evol(std::valarray<double>(Core::NaN, 2),
                            TSTART,
                            tdisk),
                a_locked_transform(a_locked_transform_coef,
                                   -Core::Inf,
                                   Core::Inf),
                a_locked_evol(a_locked_evol_coef, tdisk, tbreak),
                a_fast_transform = identity,
                a6p5_fast_evol(a6p5_fast_evol_coef, tbreak, tdeath),
                a_noplanet_transform = identity,
                nan_noplanet_evol(std::valarray<double>(Core::NaN, 2),
                                tdeath,
                                tend),

                Lconv_disk_transform = identity,
                Lconv_disk_evol(std::valarray<double>(wdisk * Ic, 1),
                                TSTART,
                                tdisk),
                Lconv_locked_term1_poly(Lconv_locked_term1_coef,
                                        -Core::Inf,
                                        Core::Inf),
                Lconv_locked_term2_poly(Lconv_locked_term2_coef,
                                        -Core::Inf,
                                        Core::Inf),
                Lconv_locked_evol = a_locked_evol,
                Lconv_fast_transform(std::valarray<double>(Core::NaN, 2),
                                     -Core::Inf,
                                     Core::Inf),
                Lconv_fast_evol(std::valarray<double>(Core::NaN, 2),
                                tbreak,
                                tdeath),
                Lconv_noplanet_transform(std::valarray<double>(Core::NaN, 2),
                                         -Core::Inf,
                                         Core::Inf),
                Lconv_noplanet_evol(std::valarray<double>(Core::NaN, 2),
                                    tdeath,
                                    tend),
                Lrad_transform = identity,
                Lrad_evol(std::valarray<double>(), TSTART, tend),
                zero_e(std::valarray<double>(), tdisk, tdeath);

            FunctionToPower
                a_fast_evol(&a6p5_fast_evol, 1.0 / 6.5),
                Lconv_locked_transform1(&Lconv_locked_term1_poly,
                                        -10.0 / 3.0),
                Lconv_locked_transform2(&Lconv_locked_term2_poly, -1.0);

            FuncPlusFunc Lconv_locked_transform(&Lconv_locked_transform1,
                                                &Lconv_locked_transform2);

            PiecewiseFunction a_evol, e_evol, Lconv_evol;

            a_evol.add_piece(&nan_disk_evol);
            a_evol.add_piece(&a_locked_evol);
            a_evol.add_piece(&a_fast_evol);
            a_evol.add_piece(&nan_noplanet_evol);

            e_evol.add_piece(&nan_disk_evol);
            e_evol.add_piece(&zero_e);
            e_evol.add_piece(&nan_noplanet_evol);

            Lconv_evol.add_piece(&Lconv_disk_evol);
            Lconv_evol.add_piece(&Lconv_locked_evol);
            Lconv_evol.add_piece(&Lconv_fast_evol);
            Lconv_evol.add_piece(&Lconv_noplanet_evol);

            std::vector< const Core::OneArgumentDiffFunction * >
                transformations(NUM_REAL_QUANTITIES - 1, &identity);
            TransformedSolution to_check(transformations, TSTART);

            transformations[SEMIMAJOR] = &a_disk_transform;
            transformations[CONV_ANGMOM] = &Lconv_disk_transform;
            to_check.add_transformation(transformations, TSTART);

            transformations[SEMIMAJOR] = &a_locked_transform;
            transformations[CONV_ANGMOM] = &Lconv_locked_transform;
            to_check.add_transformation(transformations, tdisk);

            transformations[SEMIMAJOR] = &a_fast_transform;
            transformations[CONV_ANGMOM] = &Lconv_fast_transform;
            to_check.add_transformation(transformations, tbreak);

            transformations[SEMIMAJOR] = &a_fast_transform;
            transformations[CONV_ANGMOM] = &Lconv_fast_transform;
            to_check.add_transformation(transformations, tbreak);

            transformations[SEMIMAJOR] = &a_noplanet_transform;
            transformations[CONV_ANGMOM] = &Lconv_noplanet_transform;
            to_check.add_transformation(transformations, tdeath);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            std::vector< const std::list<double> * >
                tabulated_evolution = get_evolution();

            const std::vector< const std::list<double> * > &
                transformed_evolution = to_check(tabulated_evolution);

            expected_real_quantities[SEMIMAJOR] = &a_evol;
            expected_real_quantities[ECCENTRICITY] = &e_evol;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = &Lconv_evol;
            expected_real_quantities[RAD_ANGMOM] = &zero_func;

            test_solution(transformed_evolution,
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          tend);
            delete no_evol;
        } catch (Core::Error::General &ex) {
            TEST_ASSERT_MSG(
                false,
                (
                    std::string("Unexpected exception thrown: ")
                    +
                    ex.what()+": "+ex.get_message()
                ).c_str()
            );
        } catch (std::exception &ex) {
            TEST_ASSERT_MSG(
                false,
                (
                    std::string("Unexpected exception thrown: ")
                    +
                    ex.what()
                ).c_str()
            );
        }

    }

    void test_OrbitSolver::test_polar_1_0_evolution()
    {
        StellarEvolution::MockStellarEvolution *
            no_evol = StellarEvolution::make_no_evolution();

        const double TDISK = 0.1,
                     WSTAR = 0.01,
                     WORB = 0.1;
        double initial_L = 1.0;
        std::vector<const Core::OneArgumentDiffFunction *>
            expected_real_quantities = calculate_expected_polar_1_0(
                TDISK,
                WSTAR,
                WORB
            );
        ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
        expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
        expected_evol_mode.add_break(TDISK, Core::BINARY);

        ExpectedEvolutionMode<bool> expected_wind_mode;
        expected_wind_mode.add_break(TSTART, false);

        double semimajor = (*expected_real_quantities[SEMIMAJOR])(
            (TDISK + MAX_AGE) / 2.0
        );

        make_single_component_star(*no_evol,
                                   0.0,//Kw
                                   1.0,//Wsat
                                   Core::Inf,
                                   0.9 * WSTAR,
                                   1.1 * WSTAR,
                                   0.1 * WSTAR,
                                   1.0);//tcoup

        while(true) {
            evolve(WSTAR,//wdisk
                   TDISK,//tdisk
                   semimajor,//initial semimajor
                   &initial_L,//initial L*
                   M_PI / 2.0);//initial inclination

            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          (__star ? TSTART : TDISK),
                          MAX_AGE);

            if(__star) {
                delete __star;
                delete __system;
                delete __solver;
                __star = NULL;
                __system = NULL;
                __solver = NULL;

                __primary_planet = new Planet::Planet(1.0, 1.0, 1.0);
                set_single_component_dissipation(
                    0.9 * WSTAR,
                    1.1 * WSTAR,
                    0.1 * WSTAR,
                    1.0
                );
                initial_L = WSTAR;
            } else
                break;
        }

        delete no_evol;
    }

    void test_OrbitSolver::test_polar_2_0_evolution()
    {
        StellarEvolution::MockStellarEvolution *
            no_evol = StellarEvolution::make_no_evolution();

        const double TDISK = 0.1,
                     WSTAR = 0.01,
                     WORB = 0.1,
                     PHASE_LAG = 1e-2;

        double lconv_decay_rate, semimajor, initial_L = 1.0;

        std::vector<const Core::OneArgumentDiffFunction *>
            expected_real_quantities = calculate_expected_polar_2_0(TDISK,
                                                                    WSTAR,
                                                                    WORB,
                                                                    PHASE_LAG,
                                                                    lconv_decay_rate,
                                                                    semimajor);


        std::valarray<double> identity_coef(0.0, 2);
        identity_coef[1] = 1.0;
        StellarEvolution::PolynomialEvolutionQuantity
            identity(identity_coef, -Core::Inf, Core::Inf),
            one_func(std::valarray<double>(1.0, 1), -Core::Inf, Core::Inf);

        CosFunction cos_transform(&identity);
        FuncPlusFunc one_plus_cos_transform(&one_func,
                                            &cos_transform);


        make_single_component_star(
            *no_evol,
            0.0,//Kw
            1.0,//Wsat
            Core::Inf,//tcoup
            2.0 * (WSTAR - lconv_decay_rate * MAX_AGE),
            2.0 * (WSTAR + lconv_decay_rate * MAX_AGE),
            0.1 * WSTAR,
            PHASE_LAG
        );

        while(true) {
            evolve(WSTAR,//wdisk
                   TDISK,//tdisk
                   semimajor,//initial semimajor
                   &initial_L,//initial L*
                   M_PI / 2.0);//initial inclination

            std::vector< const Core::OneArgumentDiffFunction * >
                transformations(NUM_REAL_QUANTITIES - 1, &identity);
            transformations[CONV_INCLINATION] = &one_plus_cos_transform;
            transformations[RAD_INCLINATION] = &one_plus_cos_transform;
            TransformedSolution to_check(transformations, TSTART);

            std::vector< const std::list<double> * >
                tabulated_evolution = get_evolution();
            const std::vector< const std::list<double> * > &
                transformed_evolution = to_check(tabulated_evolution);

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
            expected_evol_mode.add_break(TDISK, Core::BINARY);

            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            test_solution(transformed_evolution,
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          (__star ? TSTART : TDISK),
                          MAX_AGE);

            if(__star) {
                delete __star;
                delete __system;
                delete __solver;
                __star = NULL;
                __system = NULL;
                __solver = NULL;

                __primary_planet = new Planet::Planet(1.0, 1.0, 1.0);
                set_single_component_dissipation(
                    2.0 * (WSTAR - lconv_decay_rate * MAX_AGE),
                    2.0 * (WSTAR + lconv_decay_rate * MAX_AGE),
                    0.1 * WSTAR,
                    PHASE_LAG
                );
                initial_L = WSTAR;
            } else
                break;
        }

        delete no_evol;
    }
    void test_OrbitSolver::test_oblique_1_0_evolution()
    {
        StellarEvolution::MockStellarEvolution *
            no_evol = StellarEvolution::make_no_evolution();

        const double PHASE_LAG = 0.1,
                     TDISK = 0.1,
                     WORB = 0.1,
                     INITIAL_INC = 0.25 * M_PI,
                     INITIAL_WSTAR = 0.01;

        double initial_angmom = 1.0,
               min_wstar;

        std::vector<const Core::OneArgumentDiffFunction *>
            expected_real_quantities = calculate_expected_oblique_m_0(
                1,
                TDISK,
                WORB,
                INITIAL_INC,
                INITIAL_WSTAR,
                PHASE_LAG,
                min_wstar
            );

        double wstar_tolerance = 0.01 * (INITIAL_WSTAR - min_wstar),
               semimajor = (*expected_real_quantities[SEMIMAJOR])(
                   (TDISK + MAX_AGE) / 2.0
               );

        ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
        expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
        expected_evol_mode.add_break(TDISK, Core::BINARY);

        ExpectedEvolutionMode<bool> expected_wind_mode;
        expected_wind_mode.add_break(TSTART, false);

        make_single_component_star(
            *no_evol,
            0.0,//Kw
            1.0,//Wsat
            Core::Inf,//tcoup
            min_wstar - wstar_tolerance,
            INITIAL_WSTAR + wstar_tolerance,
            0.1 * INITIAL_WSTAR,
            PHASE_LAG
        );

        while(true) {
            evolve(INITIAL_WSTAR,//wdisk
                   TDISK,//tdisk
                   semimajor,//initial semimajor
                   &initial_angmom,//initial L*
                   INITIAL_INC);//initial inclination

            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          (__star ? TSTART : TDISK),
                          MAX_AGE);

            if(__star) {
                delete __star;
                delete __system;
                delete __solver;
                __star = NULL;
                __system = NULL;
                __solver = NULL;

                __primary_planet = new Planet::Planet(1.0, 1.0, 1.0);
                set_single_component_dissipation(
                    min_wstar - wstar_tolerance,
                    INITIAL_WSTAR + wstar_tolerance,
                    0.1 * INITIAL_WSTAR,
                    PHASE_LAG
                );

                initial_angmom = INITIAL_WSTAR;
            } else
                break;
        }
    }

    void test_OrbitSolver::test_oblique_2_0_evolution()
    {
        StellarEvolution::MockStellarEvolution *
            no_evol = StellarEvolution::make_no_evolution();

        const double PHASE_LAG = 0.1,
                     TDISK = 0.1,
                     WORB = 0.1,
                     INITIAL_INC = 0.25 * M_PI,
                     INITIAL_WSTAR = 0.01;

        double initial_angmom = 1.0,
               min_wstar;

        std::vector<const Core::OneArgumentDiffFunction *>
            expected_real_quantities = calculate_expected_oblique_m_0(
                2,
                TDISK,
                WORB,
                INITIAL_INC,
                INITIAL_WSTAR,
                PHASE_LAG,
                min_wstar
            );

        double wstar_tolerance = 0.01 * (INITIAL_WSTAR - min_wstar),
               semimajor = (*expected_real_quantities[SEMIMAJOR])(
                   (TDISK + MAX_AGE) / 2.0
               );

        ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
        expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
        expected_evol_mode.add_break(TDISK, Core::BINARY);

        ExpectedEvolutionMode<bool> expected_wind_mode;
        expected_wind_mode.add_break(TSTART, false);

        make_single_component_star(
            *no_evol,
            0.0,//Kw
            1.0,//Wsat
            Core::Inf,//tcoup
            2.0 * min_wstar - wstar_tolerance,
            2.0 * INITIAL_WSTAR + wstar_tolerance,
            0.1 * INITIAL_WSTAR,
            PHASE_LAG
        );

        while(true) {
            evolve(INITIAL_WSTAR,//wdisk
                   TDISK,//tdisk
                   semimajor,//initial semimajor
                   &initial_angmom,//initial L*
                   INITIAL_INC);//initial inclination

            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          (__star ? TSTART : TDISK),
                          MAX_AGE);
            if(__star) {
                delete __star;
                delete __system;
                delete __solver;
                __star = NULL;
                __system = NULL;
                __solver = NULL;

                __primary_planet = new Planet::Planet(1.0, 1.0, 1.0);
                set_single_component_dissipation(
                    2.0 * min_wstar - wstar_tolerance,
                    2.0 * INITIAL_WSTAR + wstar_tolerance,
                    0.1 * INITIAL_WSTAR,
                    PHASE_LAG
                );

                initial_angmom = INITIAL_WSTAR;
            } else
                break;
        }

        delete no_evol;

    }

    void test_OrbitSolver::setup()
    {
        assert(__star == NULL);
        assert(__primary_planet == NULL);
        assert(__system == NULL);
        assert(__solver == NULL);
        assert(__temp_functions.empty());
    }

    void test_OrbitSolver::tear_down()
    {
        if(__star) {
            std::cout << "Deleting star" << std::endl;
            delete __star;
            __star = NULL;
        }
        if(__primary_planet) {
            std::cout << "Deleting primary planet" << std::endl;
            delete __primary_planet;
            __primary_planet = NULL;
        }
        if(__system) {
            std::cout << "Deleting system" << std::endl;
            delete __system;
            __system = NULL;
        }
        if(__solver) {
            std::cout << "Deleting solver" << std::endl;
            delete __solver;
            __solver = NULL;
        }
        for(
            std::vector< const Core::OneArgumentDiffFunction* >::iterator
                temp_func_i = __temp_functions.begin();
            temp_func_i != __temp_functions.end();
            ++temp_func_i
        ) {
            std::cout << "Deleting function" << std::endl;
            delete *temp_func_i;
        }
        __temp_functions.clear();
    }

    test_OrbitSolver::test_OrbitSolver() :
        __solver(NULL),
        __system(NULL),
        __star(NULL),
        __primary_planet(NULL)
    {
        TEST_ADD(test_OrbitSolver::test_disk_locked_no_stellar_evolution);
        TEST_ADD(test_OrbitSolver::test_disk_locked_with_stellar_evolution);
        TEST_ADD(test_OrbitSolver::test_no_planet_evolution);
        TEST_ADD(test_OrbitSolver::test_unlocked_evolution);
//        TEST_ADD(test_OrbitSolver::test_locked_evolution);//NOT REVIVED!!!
        TEST_ADD(test_OrbitSolver::test_disklocked_to_locked_to_noplanet);
        TEST_ADD(test_OrbitSolver::test_disklocked_to_fast_to_noplanet);
        TEST_ADD(test_OrbitSolver::test_disklocked_to_fast_to_locked);
        TEST_ADD(test_OrbitSolver::test_disklocked_to_locked_to_fast);
        TEST_ADD(test_OrbitSolver::test_polar_1_0_evolution);
        TEST_ADD(test_OrbitSolver::test_polar_2_0_evolution);
        TEST_ADD(test_OrbitSolver::test_oblique_1_0_evolution);
        TEST_ADD(test_OrbitSolver::test_oblique_2_0_evolution);
    }

}//End Evolve namespace.
