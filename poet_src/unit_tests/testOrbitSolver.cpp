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

    void TransformedSolution::add_transformation(
        const std::vector<const Core::OneArgumentDiffFunction *> &transforms,
        double change_age
    )
    {
        for(unsigned q = 0; q < NUM_REAL_QUANTITIES - 1; ++q)
            __transforms[q].push_back(transforms[q]);

        __change_ages.push_back(change_age);
    }

    const std::vector< const std::list<double> * > &
        TransformedSolution::operator()(
            const std::vector< const std::list<double> * > &solution
        )
    {
        if(__transformed_orbit[AGE]) delete __transformed_orbit[AGE];
        __transformed_orbit[AGE] = solution[AGE];

        for(unsigned q = 0; q < NUM_REAL_QUANTITIES - 1; ++q) {
            std::list<double>::const_iterator
                change_i = __change_ages.begin();
            change_i++;
            std::list<const Core::OneArgumentDiffFunction *>::const_iterator
                transform_i = __transforms[q].begin();
            std::list<double> *transformed_quantity = new std::list<double>();
            for(
                std::list<double>::const_iterator
                    var_i = solution[q]->begin(),
//                    deriv_i = deriv->begin(),
                    t_i = solution[AGE]->begin();
                var_i != solution[q]->end();
                ++var_i,
//                ++deriv_i,
                ++t_i
            ) {
                if(change_i != __change_ages.end() && *t_i >= *change_i) {
                    ++transform_i;
                    ++change_i;
                }
                transformed_quantity->push_back((**transform_i)(*var_i));
/*                const FunctionDerivatives *dvar = (*transform_i)->deriv(*var_i);
                __transformed_orbit[var_type].push_back(dvar->order(0));
                __transformed_deriv[var_type].push_back(dvar->order(1)
                                                        *
                                                        (*deriv_i));*/
            }
            if(__transformed_orbit[q]) delete __transformed_orbit[q]; 
            __transformed_orbit[q] = transformed_quantity;
        }
        return __transformed_orbit;
    }

    TransformedSolution::~TransformedSolution()
    {
        for(unsigned q = 0; q < NUM_REAL_QUANTITIES - 1; ++q)
            if(__transformed_orbit[q]) delete __transformed_orbit[q];
    }

    void test_OrbitSolver::make_const_lag_star(
        const StellarEvolution::Interpolator &evolution,
        double wind_strength,
        double wind_sat_freq,
        double coupling_timescale,
        double phase_lag
    )
    {
        __star = new Star::InterpolatedEvolutionStar(1.0,//mass
                                                     0.0,//feh
                                                     wind_strength,
                                                     wind_sat_freq,
                                                     coupling_timescale,
                                                     evolution);
        Evolve::BrokenPowerlawPhaseLagZone *zone = &(__star->envelope());
        for(int i = 0; i < 2; ++i) {
            zone->setup(std::vector<double>(),//Wtide breaks
                        std::vector<double>(),//W* breaks
                        std::vector<double>(1, 0.0),//Wtide pow.
                        std::vector<double>(1, 0.0),//W* pow.
                        (i==0 ? phase_lag : 0.0));
            zone = &(__star->core());
        }
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
        //phase lag(min_frequency - decay_scale)
        //= 
        //phase lag(max_frequency + decay_scale)
        //=
        //suppression_factor * phase_lag
        const double suppression_factor = 0.01;

        __star = new Star::InterpolatedEvolutionStar(1.0,//mass
                                                     0.0,//feh
                                                     wind_strength,
                                                     wind_sat_freq,
                                                     coupling_timescale,
                                                     evolution);
        __star->core().setup(std::vector<double>(),//Wtide breaks
                             std::vector<double>(),//W* breaks
                             std::vector<double>(1, 0.0),//Wtide pow.
                             std::vector<double>(1, 0.0),//W* pow.
                             0.0);//Phase lag (no dissipation in the core).

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

        __star->envelope().setup(breaks,
                                 std::vector<double>(),
                                 powerlaw_indices,
                                 std::vector<double>(1, 0.0),
                                 phase_lag);
    }

    StellarEvolution::MockStellarEvolution *
        test_OrbitSolver::make_no_evolution(double Rstar, double Iconv)
    {
        return new StellarEvolution::MockStellarEvolution(
            0.0,
            std::valarray< std::valarray<double> >(//R
                std::valarray<double>(Rstar, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Iconv
                std::valarray<double>(Iconv, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Irad
                std::valarray<double>(1.0, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Rcore
                std::valarray<double>(1.0, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Mcore
                std::valarray<double>(1.0, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Lum
                std::valarray<double>(1.0, 1),
                1
            )
        );
    }

    StellarEvolution::MockStellarEvolution *
        test_OrbitSolver::make_linear_I_evolution()
    {
        return new StellarEvolution::MockStellarEvolution(
            0.0,
            std::valarray< std::valarray<double> >(//R
                std::valarray<double>(1.0, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Iconv
                std::valarray<double>(1.0, 1),
                2
            ),
            std::valarray< std::valarray<double> >(//Irad
                std::valarray<double>(1.0, 1),
                2
            ),
            std::valarray< std::valarray<double> >(//Rcore
                std::valarray<double>(1.0, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Mcore
                std::valarray<double>(1.0, 1),
                1
            ),
            std::valarray< std::valarray<double> >(//Lum
                std::valarray<double>(1.0, 1),
                1
            )
        );
    }

    void test_OrbitSolver::evolve(double wdisk,
                                  double tdisk,
                                  double initial_a,
                                  const double *initial_Lstar,
                                  double initial_incl,
                                  double planet_mass,
                                  double tplanet,
                                  double max_age,
                                  double planet_radius,
                                  double precision,
                                  double max_step_factor)
    {
        __star->select_interpolation_region(TSTART);
        

        if(std::isnan(tplanet)) tplanet = tdisk;

        Planet::LockedPlanet planet(planet_mass,
                                    (planet_mass ? planet_radius : 0.0));
        planet.configure(true, //init
                         tplanet, //age
                         __star->mass(), //mass
                         initial_a, //semimajor
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
            initial_a, //semimajor
            0.0, //eccentricity
            initial_incl, //inclination
            wdisk, //Wdisk
            tdisk, //disk dissipation age
            tplanet //planet formation age
        );
        if(tdisk <= TSTART) {
            double zeros[] = {0.0, 0.0};
            if(tplanet <= TSTART) {
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

            tabulated_real_quantities[CONV_INCLINATION] = &(
                __star->envelope().get_evolution_real(Evolve::INCLINATION)
            );

            tabulated_real_quantities[RAD_INCLINATION] = 
                &(__star->core().get_evolution_real(Evolve::INCLINATION));

            tabulated_real_quantities[CONV_PERIAPSIS] = 
                &(__star->envelope().get_evolution_real(Evolve::PERIAPSIS));

            tabulated_real_quantities[RAD_PERIAPSIS] = 
                &(__star->core().get_evolution_real(Evolve::PERIAPSIS));

            tabulated_real_quantities[CONV_ANGMOM] = &(
                __star->envelope().get_evolution_real(
                    Evolve::ANGULAR_MOMENTUM
                )
            );

            tabulated_real_quantities[RAD_ANGMOM] = &(
                __star->core().get_evolution_real(Evolve::ANGULAR_MOMENTUM)
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

        const std::list<bool> &tabulated_wind_sat =
            __star->wind_saturation_evolution();

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

        msg << tabulated_modes.size() << " tabulated modes, "
            << tabulated_wind_sat.size() << " tabulated wind saturations";

        all_same_size = (all_same_size
                         &&
                         tabulated_modes.size() == num_ages
                         &&
                         tabulated_wind_sat.size() == num_ages);


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
            real_tabulated_iter[q] = tabulated_real_quantities[q]->begin();
        std::list<Core::EvolModeType>::const_iterator
            tabulated_mode_iter = tabulated_modes.begin();
        std::list<bool>::const_iterator
            tabulated_wind_sat_iter = tabulated_wind_sat.begin();
        double last_checked_age = Core::NaN;

        for(
            std::list<double>::const_iterator
                age_i = tabulated_real_quantities[AGE]->begin();
            age_i != tabulated_real_quantities[AGE]->end();
            ++age_i
        ) {
            std::vector<double> expected_real_values(AGE);
            for(unsigned q = 0; q < AGE; ++q) {
                std::cerr << "Getting "
                          << static_cast<RealEvolutionQuantity>(q)
                          << std::endl
                          << "expected @ "
                          << expected_real_quantities[q]
                          << std::endl;
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
                          << ", mode = " << *tabulated_mode_iter
                          << ", wind is ";
            if(!(*tabulated_wind_sat_iter)) age_msg_start << " not ";
            age_msg_start << "saturated";
            for(unsigned q = 0; q < AGE; ++q)
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
            for(unsigned q = 0; q < AGE; ++q) ++(real_tabulated_iter[q]);
            ++tabulated_wind_sat_iter;
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

            make_const_lag_star(stellar_evol,
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
            delete __star;
            delete __system;
            delete __solver;

            if(initial_Lstar[0] == 0) return;

            for(double phase_lag = 0.0; phase_lag < 1.5; phase_lag += 1.0)
                for(
                    double mplanet = 0.0;
                    mplanet < 1.5 - phase_lag;
                    mplanet += 1.0
                ) {
                    make_const_lag_star(stellar_evol,
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

                    delete __star;
                    delete __system;
                    delete __solver;
                }
    }

    void test_OrbitSolver::test_disk_locked_no_stellar_evolution()
    {
        try {
            StellarEvolution::MockStellarEvolution *no_evol =
                make_no_evolution();
            make_const_lag_star(*no_evol,
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
            delete __solver;
            delete __system;
            delete __star;
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
                make_linear_I_evolution();

            make_const_lag_star(*evol1, 1.0, 1.0, 1.0);

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

            make_const_lag_star(evol2, 1.0, 1.0, 1.0);

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
            delete __solver;
            delete __system;
            delete __star;
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
                *stellar_evol = make_no_evolution();
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
            stellar_evol = make_linear_I_evolution();
            
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
            std::cerr << "Lc1(" << TSTART << ") = "
                      << Lc1(TSTART)
                      << std::endl;
            std::cerr << "Lc2(" << TSTART << ") = "
                      << Lc2(TSTART)
                      << std::endl;
            std::cerr << "Lconv(" << TSTART << ") = "
                      << (*(expected_real_quantities[CONV_ANGMOM]))(TSTART)
                      << std::endl;
            delete stellar_evol;
            stellar_evol = make_no_evolution();
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
                *no_evol = make_no_evolution(1.0);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            const double mplanet = 100;
            double lag = 1e-8 / mplanet,
                   mplanet_si = mplanet * Core::AstroConst::jupiter_mass,
                   mstar_si = Core::AstroConst::solar_mass,
                   alpha = (
                       -2.4 * M_PI
                       *
                       std::sqrt(
                           Core::AstroConst::G * (mplanet_si + mstar_si)
                           /
                           Core::AstroConst::solar_radius
                       )
                       *
                       mplanet_si / mstar_si
                       *
                       lag * Core::AstroConst::Gyr
                       /
                       Core::AstroConst::solar_radius
                   ),
                   a6p5_offset = std::pow(2.0, 6.5) - 6.5 * alpha,
                   a0 = std::pow(a6p5_offset, 1.0 / 6.5),
                   Lscale = (
                       -Core::AstroConst::jupiter_mass*mplanet
                       /
                       std::pow(Core::AstroConst::solar_radius, 1.5)
                       *
                       std::sqrt(
                           Core::AstroConst::G
                           /
                           (
                               mplanet * Core::AstroConst::jupiter_mass
                               +
                               Core::AstroConst::solar_mass
                           )
                       )
                       *
                       Core::AstroConst::day
                   );

            std::valarray<double> a6p5_poly_coef(2);
            a6p5_poly_coef[0] = a6p5_offset;
            a6p5_poly_coef[1] = 6.5 * alpha;
            StellarEvolution::PolynomialEvolutionQuantity a6p5_evol(
                a6p5_poly_coef,
                TSTART,
                1.0
            );
            FunctionToPower sqrta_evol(&a6p5_evol, 1.0/13.0);
            ExponentialPlusFunc Lconv_unscaled(&sqrta_evol,
                                               -std::sqrt(a0),
                                               0);

            expected_real_quantities[SEMIMAJOR] = new FunctionToPower(
                &a6p5_evol,
                1.0 / 6.5
            );
            expected_real_quantities[ECCENTRICITY] = &zero_func;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = new ScaledFunction(
                &Lconv_unscaled,
                Lscale
            );
            expected_real_quantities[RAD_ANGMOM] = &zero_func;

            double initial_a = (*expected_real_quantities[SEMIMAJOR])(TSTART);

            std::valarray<double> initial_L(3);
            initial_L[0] = (*expected_real_quantities[CONV_ANGMOM])(
                TSTART
            );
            initial_L[1] = 0.0;
            initial_L[2] = 0.0;

            make_const_lag_star(*no_evol,
                                0.0,//wind K
                                100.0,//wsat
                                Core::Inf,//tcoup
                                lag);

            evolve(0.0,//wdisk
                   0.0,//tdisk
                   (*expected_real_quantities[SEMIMAJOR])(TSTART),
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
                          1.0,
                          true);

            delete __system;
            delete __solver;
            delete expected_real_quantities[SEMIMAJOR];
            delete expected_real_quantities[CONV_ANGMOM];

            a0 = 2.6;
            a6p5_offset = std::pow(a0, 6.5);
            a6p5_poly_coef[0] = a6p5_offset;
            a6p5_poly_coef[1] = -6.5*alpha;
            StellarEvolution::PolynomialEvolutionQuantity a6p5_evol_slow(
                a6p5_poly_coef,
                TSTART,
                1.0
            );
            FunctionToPower a_evol_slow(&a6p5_evol_slow, 1.0 / 6.5),
                            sqrta_evol_slow(&a6p5_evol_slow, 1.0 / 13.0);
            ExponentialPlusFunc Lconv_unscaled_slow(&sqrta_evol_slow,
                                                    -1e5-std::sqrt(a0),
                                                    0);
            expected_real_quantities[SEMIMAJOR] = &a_evol_slow;
            expected_real_quantities[CONV_ANGMOM] = new ScaledFunction(
                &Lconv_unscaled_slow,
                Lscale
            );

            initial_a = (*expected_real_quantities[SEMIMAJOR])(TSTART);
            initial_L[0] = (*expected_real_quantities[CONV_ANGMOM])(TSTART);

            evolve(0.0,//wdisk
                   0.0,//tdisk
                   (*expected_real_quantities[SEMIMAJOR])(TSTART),
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

            delete __star;
            delete __system;
            delete __solver;
            delete expected_real_quantities[CONV_ANGMOM];
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
                no_evol = make_no_evolution(1.0, Ic);

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
            double solver_params[] = {Lscale, beta, kappa, int_const, 0.0};
/*            double a0 = solve(a1,
                              0.0,
                              1e-9,
                              &locked_unsat_eq,
                              &locked_unsat_deriv,
                              &locked_unsat_eq_deriv,
                              static_cast<void*>(solver_params));*/
            solver_params[4] = TSTART;
/*            double astart = solve(a0,
                                  0.0,
                                  1e-9,
                                  &locked_unsat_eq,
                                  &locked_unsat_deriv,
                                  &locked_unsat_eq_deriv,
                                  static_cast<void*>(solver_params));*/
            solver_params[2] = kappa_s;
            solver_params[3] = int_const_s;
            solver_params[4] = 0;
/*            double a0_s = solve(a1,
                                0.0,
                                1e-9,
                                &locked_sat_eq,
                                &locked_sat_deriv,
                                &locked_sat_eq_deriv,
                                static_cast<void*>(solver_params));*/
            solver_params[4] = TSTART;
/*            double astart_s = solve(a1,
                                    0.0,
                                    1e-9,
                                    &locked_sat_eq,
                                    &locked_sat_deriv,
                                    &locked_sat_eq_deriv,
                                    static_cast<void*>(solver_params));*/
            solver_params[4] = 1.0;

    /*		double w0=std::sqrt(
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
                no_evol = make_no_evolution(1.0, Ic);

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

            make_const_lag_star(*no_evol,
                                Kwind,
                                wsat,
                                Core::Inf,//core-env coupling timescale
                                lag_from_lgQ(1, (Core::AstroConst::jupiter_mass
                                                 /
                                                 Core::AstroConst::solar_mass)));
            Planet::LockedPlanet planet(1.0, 1.0);
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


            std::cerr << "C_L = " << Lscale << std::endl;
            std::cerr << "beta = " << beta << std::endl;
            std::cerr << "C_i = " << int_const << std::endl;
            std::cerr << "kappa = " << kappa << std::endl;
            std::cerr << "Initial semimajor axis: " << ainitial << std::endl;
            std::cerr << "Destruction semimajor: " << adestr << std::endl;
            std::cerr << "Destruction age: " << tdestr << std::endl;

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
                        (tfinal - __system->age()) / 1000.0, //time step
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
                          tfinal,
                          true);

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
            delete __star;
            delete __system;
            delete __solver;
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
                no_evol = make_no_evolution(1.0, I_CONV);

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

            make_const_lag_star(
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

            std::cerr << "a0 = " << a_formation << std::endl
                      << "a_death = " << adestr << std::endl
                      << "formation age = " << TDISK << std::endl
                      << "destruction age = " << TDESTR << std::endl; 

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
                          MAX_AGE,
                          true);

            delete __star;
            delete __system;
            delete __solver;
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
                         Q = std::pow(10.0, lgQ),
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
                         tdisk = 1,
                         async = 2.5,
                         tsync = 2.0,
                         tend = 3,
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
    		std::cout << std::endl << "arate=" << -alpha
                << std::endl << "alpha=" << Lscale
                << std::endl << "beta=" << beta
                << std::endl << "Ic=" << Ic
                << std::endl << "wdisk=" << wdisk
                << std::endl << "wlocked=" << wlocked
                << std::endl << "a0=" << a_formation
                << std::endl;

            StellarEvolution::MockStellarEvolution *
                no_evol = make_no_evolution(1.0, Ic);

            make_const_lag_star(
                *no_evol,//evolution
                0.0,//Kwind
                100.0,//wsat
                Core::Inf,//tcoup
                lag_from_lgQ(lgQ,
                             (Core::AstroConst::jupiter_mass
                              /
                              Core::AstroConst::solar_mass))
            );
            evolve(wdisk,//wdisk,
                   tdisk,//tdisk
                   a_formation,//initial semimajor
                   &zero,//initial L*
                   0.0,//initial inclination
                   1.0,//planet_mass
                   Core::NaN,//form the planet when disk dissipates
                   tend,//max evolution age
                   0.01);//planet radius

            ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
            expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
            expected_evol_mode.add_break(tdisk, Core::BINARY);

            std::valarray<double> a6p5_poly_coef(2);
            a6p5_poly_coef[0] = a6p5_offset;
            a6p5_poly_coef[1] = 6.5 * alpha;
            StellarEvolution::PolynomialEvolutionQuantity
                a6p5_evol(a6p5_poly_coef, tdisk, tsync),
                Lconv_disk(std::valarray<double>(Ic * wdisk, 1),
                           TSTART,
                           tdisk),
                Lconv_locked(std::valarray<double>(Ic * wlocked, 1),
                             tsync,
                             tend),
                nan_disk(std::valarray<double>(Core::NaN, 2),
                       TSTART,
                       tdisk),
                a_locked(std::valarray<double>(async, 1),
                         tsync,
                         tend),
                Lrad_evol(std::valarray<double>(),
                          TSTART,
                          tend),
                zero_e(std::valarray<double>(), tdisk, tend);
            FunctionToPower a_fast(&a6p5_evol, 1.0 / 6.5),
                            sqrta_evol(&a6p5_evol, 1.0 / 13.0);
            ExponentialPlusFunc Lconv_unscaled(
                &sqrta_evol,
                -Ic * wdisk / Lscale - std::sqrt(a_formation), 0
            );
            ScaledFunction Lconv_fast(&Lconv_unscaled, -Lscale);

            PiecewiseFunction a_evol, e_evol, Lconv_evol;

            a_evol.add_piece(&nan_disk);
            a_evol.add_piece(&a_fast);
            a_evol.add_piece(&a_locked);

            e_evol.add_piece(&nan_disk);
            e_evol.add_piece(&zero_e);

            Lconv_evol.add_piece(&Lconv_disk);
            Lconv_evol.add_piece(&Lconv_fast);
            Lconv_evol.add_piece(&Lconv_locked);

            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities(NUM_REAL_QUANTITIES - 1);

            expected_real_quantities[SEMIMAJOR] = &a_evol;
            expected_real_quantities[ECCENTRICITY] = &e_evol;
            expected_real_quantities[CONV_INCLINATION] = &zero_func;
            expected_real_quantities[RAD_INCLINATION] = &zero_func;
            expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
            expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
            expected_real_quantities[CONV_ANGMOM] = &Lconv_evol;
            expected_real_quantities[RAD_ANGMOM] = &zero_func;


            ExpectedEvolutionMode<bool> expected_wind_mode;
            expected_wind_mode.add_break(TSTART, false);

            test_solution(get_evolution(),
                          expected_real_quantities,
                          expected_evol_mode,
                          expected_wind_mode,
                          TSTART,
                          tend,
                          true);

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
    		std::cout << std::endl << "a0=" << a0
                << std::endl << "alphaL=" << alphaL
                << std::endl << "beta=" << beta
                << std::endl << "gamma=" << gamma
                << std::endl << "Ic=" << Ic
                << std::endl << "wdisk=" << wdisk
                << std::endl << "kappa=" << kappa
                << std::endl << "Kwind=" << Kwind
                << std::endl << "Q=" << Q
                << std::endl << "Rp=" << Rp
                << std::endl;

            StellarEvolution::MockStellarEvolution *
                no_evol = make_no_evolution(1.0, Ic);

            make_const_lag_star(
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
                          tend,
                          true);
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
            no_evol = make_no_evolution();

        const double TDISK = 0.1,
                     WSTAR = 0.01,
                     WORB = 0.1,
                     AORB = std::pow(
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
                             WORB / Core::AstroConst::day,
                             2
                         ),
                         1.0 / 3.0
                     ) / Core::AstroConst::solar_radius,
                     ONE = 1.0;

        std::cerr << "Polar orbit inertial mode only:" << std::endl
                  << "TDISK = " << TDISK << std::endl
                  << "W* = " << WSTAR << std::endl
                  << "WORB = " << WORB << std::endl
                  << "a = " << AORB << std::endl;

        make_single_component_star(*no_evol,
                                   0.0,//Kw
                                   1.0,//Wsat
                                   Core::Inf,
                                   0.9 * WSTAR,
                                   1.1 * WSTAR,
                                   0.1 * WSTAR);//tcoup

        evolve(WSTAR,//wdisk
               TDISK,//tdisk
               AORB,//initial semimajor
               &ONE,//initial L*
               M_PI / 2.0);//initial inclination

        StellarEvolution::PolynomialEvolutionQuantity
            disk_nan_evol(std::valarray<double>(Core::NaN, 1),
                          TSTART,
                          TDISK),
            fixed_a_evol(std::valarray<double>(AORB, 1),
                         TDISK,
                         MAX_AGE),
            fixed_e_evol(std::valarray<double>(),
                         TDISK,
                         MAX_AGE),
            disk_zero_evol(std::valarray<double>(),
                           TSTART,
                           TDISK),
            halfpi_evol(std::valarray<double>(M_PI / 2.0, 1),
                        TDISK,
                        MAX_AGE);

        PiecewiseFunction a_evol, e_evol, conv_incl_evol, rad_incl_evol;

        a_evol.add_piece(&disk_nan_evol);
        a_evol.add_piece(&fixed_a_evol);

        e_evol.add_piece(&disk_nan_evol);
        e_evol.add_piece(&fixed_e_evol);

        conv_incl_evol.add_piece(&disk_zero_evol);
        conv_incl_evol.add_piece(&halfpi_evol);

        rad_incl_evol.add_piece(&disk_zero_evol);
        rad_incl_evol.add_piece(&halfpi_evol);

        std::vector<const Core::OneArgumentDiffFunction *>
            expected_real_quantities(NUM_REAL_QUANTITIES - 1);

        expected_real_quantities[SEMIMAJOR] = &a_evol;
        expected_real_quantities[ECCENTRICITY] = &e_evol;
        expected_real_quantities[CONV_INCLINATION] = &conv_incl_evol;
        expected_real_quantities[RAD_INCLINATION] = &rad_incl_evol;
        expected_real_quantities[CONV_PERIAPSIS] = &zero_func;
        expected_real_quantities[RAD_PERIAPSIS] = &zero_func;
        expected_real_quantities[CONV_ANGMOM] =
            new StellarEvolution::PolynomialEvolutionQuantity(
                std::valarray<double>(WSTAR, 1),
                TSTART,
                MAX_AGE
            );
        expected_real_quantities[RAD_ANGMOM] = &one_func;

        ExpectedEvolutionMode<Core::EvolModeType> expected_evol_mode;
        expected_evol_mode.add_break(TSTART, Core::LOCKED_SURFACE_SPIN);
        expected_evol_mode.add_break(TDISK, Core::BINARY);

        ExpectedEvolutionMode<bool> expected_wind_mode;
        expected_wind_mode.add_break(TSTART, false);

        test_solution(get_evolution(),
                      expected_real_quantities,
                      expected_evol_mode,
                      expected_wind_mode,
                      TSTART,
                      MAX_AGE);

        delete __star;
        delete __system;
        delete __solver;
        delete no_evol;
        delete expected_real_quantities[CONV_ANGMOM];
    }

    test_OrbitSolver::test_OrbitSolver()
    {
/*        TEST_ADD(test_OrbitSolver::test_disk_locked_no_stellar_evolution);
        TEST_ADD(test_OrbitSolver::test_disk_locked_with_stellar_evolution);
        TEST_ADD(test_OrbitSolver::test_no_planet_evolution);
        TEST_ADD(test_OrbitSolver::test_unlocked_evolution);*/
//        TEST_ADD(test_OrbitSolver::test_locked_evolution);//NOT REVIVED!!!
/*        TEST_ADD(test_OrbitSolver::test_disklocked_to_locked_to_noplanet);
        TEST_ADD(test_OrbitSolver::test_disklocked_to_fast_to_noplanet);
        TEST_ADD(test_OrbitSolver::test_disklocked_to_fast_to_locked);
        TEST_ADD(test_OrbitSolver::test_disklocked_to_locked_to_fast);*/
        TEST_ADD(test_OrbitSolver::test_polar_1_0_evolution);
    }

}//End Evolve namespace.

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
