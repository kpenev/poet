/**\file
 *
 * \brief The definition of some of the methods of the StellarSystem class.
 *
 * \ingroup StellarSystem_group
 */

#define BUILDING_LIBRARY
#include "BinarySystem.h"

namespace Evolve {

    void BinarySystem::find_locked_zones()
    {
        __locked_zones.clear();
        for(short body_ind = 0; body_ind < 2; ++body_ind) {
            DissipatingBody &body = (body_ind == 0 ? __body1 : __body2);
            for(
                unsigned zone_ind = 0;
                zone_ind < body.number_zones();
                ++zone_ind
            )
                if(body.zone(zone_ind).locked()) {
                    body.zone(zone_ind).locked_zone_index() =
                        __locked_zones.size();
                    __locked_zones.push_back(zone_ind);
                }
        }
    }

    BinarySystem::lock_scenario_type BinarySystem::find_synchronized_zones(
        double precision
    )
    {
        lock_scenario_type result;
        double worb = orbital_frequency();
        for(unsigned short body_ind = 0; body_ind < 2; ++body_ind) {
            DissipatingBody &body = (body_ind == 0 ? __body1 : __body2);
            for(
                unsigned zone_ind = 0;
                zone_ind < body.number_zones();
                ++zone_ind
            ) {
                DissipatingZone &zone = body.zone(zone_ind);
                if(!zone.can_lock())
                    continue;
                double wspin = zone.spin_frequency();
                bool other_lock=false;
                unsigned scenario_zone_i = (zone_ind
                                            +
                                            body_ind * __body1.number_zones());
                do {
                    const SpinOrbitLockInfo &lock = zone.lock_monitored(
                        other_lock
                    );
                    if(
                        std::abs(
                            (
                                lock.orbital_frequency_multiplier() * worb
                                -
                                lock.spin_frequency_multiplier() * wspin
                            )
                            /
                            worb
                        ) < precision
                    ) {
                        if(other_lock)
                            zone.swap_monitored_locks();
                        result.push_back(
                            std::make_tuple(
                                scenario_zone_i,
                                lock.lock_direction(),
                                zone.angular_momentum()
                            )
                        );
                    } else
                        assert(!zone.locked()
                               ||
                               std::get<0>(result.back()) == scenario_zone_i);
                    other_lock = !other_lock;
                } while(other_lock);
            }
        }
        return result;
    }

    void BinarySystem::set_lock_scenario(
        const lock_scenario_type &lock_scenario
    )
    {
        assert(__evolution_mode == Core::BINARY);
#ifdef VERBOSE_DEBUG
        std::cerr << "Setting lock scenario: " << std::endl;
        describe_lock_scenario(std::cerr,
                               lock_scenario,
                               std::vector<bool>(lock_scenario.size(), false),
                               true);
#endif

        std::vector<double>
            angmom(number_zones()),
            inclinations(number_zones()),
            periapses(number_zones() - 1);

        std::vector<double>::iterator
            angmom_dest = angmom.begin(),
            inclination_dest = inclinations.begin(),
            periapsis_dest = periapses.begin();
        lock_scenario_type::const_iterator
            scenario_zone_iter = lock_scenario.begin();

        unsigned scenario_zone;
        short scenario_lock_dir;
        double scenario_angmom, zone_angmom;
        std::tie(scenario_zone,
                 scenario_lock_dir,
                 scenario_angmom) = *scenario_zone_iter;

        for(short body_i = 0; body_i < 2; ++body_i) {
            DissipatingBody &body = (body_i == 0 ? __body1 : __body2);
            for(unsigned zone_i = 0; zone_i < body.number_zones(); ++zone_i) {
                DissipatingZone &zone = body.zone(zone_i);
                const SpinOrbitLockInfo &current_lock = zone.lock_monitored();
                if(
                        zone_i + body_i * __body1.number_zones()
                        ==
                        scenario_zone
                ) {
                    zone_angmom = scenario_angmom;
                    if(
                            current_lock.lock_direction()
                            !=
                            scenario_lock_dir
                    ) {
                        if(zone.locked()) {
                            body.unlock_zone_spin(zone_i, scenario_lock_dir);
#ifdef VERBOSE_DEBUG
                            std::cerr << "Unlocking "
                                      << (body_i == 0 ? "primary" : "secondary")
                                      << " zone " << zone_i
                                      << "(overall zone " << scenario_zone
                                      << ")"
                                      << std::endl;
#endif
                        } else {
                            body.lock_zone_spin(
                                zone_i,
                                current_lock.orbital_frequency_multiplier(),
                                current_lock.spin_frequency_multiplier()
                            );
#ifdef VERBOSE_DEBUG
                            std::cerr << "Locking "
                                      << (body_i == 0 ? "primary" : "secondary")
                                      << " zone " << zone_i
                                      << "(overall zone " << scenario_zone
                                      << ")"
                                      << std::endl;
#endif

                            if(scenario_lock_dir != 0) {
#ifdef VERBOSE_DEBUG
                                std::cerr << "Unlocking "
                                          << (body_i == 0 ? "primary" : "secondary")
                                          << " zone " << zone_i
                                          << "(overall zone " << scenario_zone
                                          << ")"
                                          << std::endl;
#endif
                                body.unlock_zone_spin(zone_i,
                                                      scenario_lock_dir);
                            }
                        }
                    }
#ifdef VERBOSE_DEBUG
                    else {
                        std::cerr << (body_i == 0 ? "Primary" : "Secondary")
                                  << " zone " << zone_i
                                  << "(overall zone " << scenario_zone
                                  << ") already matches scenario."
                                  << std::endl;
                    }
#endif
                    ++scenario_zone_iter;
                    if(scenario_zone_iter != lock_scenario.end())
                        std::tie(scenario_zone,
                                 scenario_lock_dir,
                                 scenario_angmom) = *scenario_zone_iter;

                } else {
                    zone_angmom = zone.angular_momentum();
                }

                if(!zone.locked())
                    *angmom_dest++ = zone_angmom;
                *inclination_dest++ = zone.inclination();
                if(body_i || zone_i)
                    *periapsis_dest++ = zone.periapsis();
            }
        }
        configure(false,
                  __age,
                  __semimajor,
                  __eccentricity,
                  angmom.data(),
                  inclinations.data(),
                  periapses.data(),
                  Core::BINARY);
    }

    void BinarySystem::unlock_all_zones(
        const std::vector<short> &unlock_directions,
        const std::vector<double> &original_angmom
    )
    {
        std::vector<double>
            config_angmom(number_zones()),
            config_inclinations(number_zones()),
            config_periapses(number_zones() - 1);
        unsigned locked_i = 0;
        for(unsigned config_i = 0; config_i < number_zones(); ++config_i) {
            DissipatingBody &body = (config_i < __body1.number_zones()
                                     ? __body1
                                     : __body2);
            unsigned body_zone_i = (
                config_i < __body1.number_zones()
                ? config_i
                : config_i - __body1.number_zones()
            );
            DissipatingZone &zone = body.zone(body_zone_i);

            config_inclinations[config_i] = zone.inclination();
            if(config_i)
                config_periapses[config_i - 1] = zone.periapsis();

            if(zone.locked()) {
                config_angmom[config_i] = (original_angmom.size()
                                           ? original_angmom[locked_i]
                                           : zone.angular_momentum());
                body.unlock_zone_spin(body_zone_i,
                                      unlock_directions[locked_i]);
                ++locked_i;
            } else
                config_angmom[config_i] = zone.angular_momentum();
        }
        configure(false,
                  __age,
                  __semimajor,
                  __eccentricity,
                  config_angmom.data(),
                  config_inclinations.data(),
                  config_periapses.data(),
                  Core::BINARY);
    }

    bool BinarySystem::test_synchronized_unlocked_zone(unsigned test_zone_ind)
    {
        bool primary_zone = test_zone_ind < __body1.number_zones();
        unsigned body_zone_ind = (primary_zone
                                  ? test_zone_ind
                                  : test_zone_ind - __body1.number_zones());
        const SpinOrbitLockInfo &lock = (
            primary_zone
            ? __body1
            : __body2
        ).zone(
            body_zone_ind
        ).lock_monitored();
        std::valarray<double> orbit;
        fill_orbit(orbit);
        std::valarray<double> derivatives(orbit.size());
        differential_equations(__age,
                               &(orbit[0]),
                               Core::BINARY,
                               &(derivatives[0]));

        std::valarray<double> check_cond_deriv;
        SynchronizedCondition(
            lock,
            primary_zone,
            body_zone_ind,
            *this
        )(
            Core::BINARY,
            orbit,
            derivatives,
            check_cond_deriv
        );
        assert(check_cond_deriv.size() == 1);
#ifndef NDEBUG
        std::cerr << "Lock dir " << lock.lock_direction()
                  << "for zone " << test_zone_ind
                  << " (condition deriv: " << check_cond_deriv[0]
                  << (lock.lock_direction() * check_cond_deriv[0] > 0
                      ? ") is good"
                      : ") is bad")
                  << std::endl;
#endif
        return lock.lock_direction() * check_cond_deriv[0] < 0;
    }

#ifndef NDEBUG
    void BinarySystem::describe_lock_scenario(
        std::ostream &os,
        const lock_scenario_type &lock_scenario,
        const std::vector<bool> passed,
        bool add_header
    )
    {
        lock_scenario_type::const_iterator
            scenario_zone_i = lock_scenario.begin();
        std::vector<bool>::const_iterator passed_i = passed.begin();
        unsigned test_zone_ind;
        short lock_direction;
        std::tie(test_zone_ind, lock_direction, std::ignore) = *scenario_zone_i;

        if(add_header)
            os  << "||" << std::setw(9 * __body1.number_zones() - 1) << "Body1"
                << "||" << std::setw(9 * __body2.number_zones() - 1) << "Body2"
                << "||" << std::endl
                << std::string(
                    9 * (__body1.number_zones() + __body2.number_zones()) + 4,
                    '_'
                )
                << std::endl
                << "||";

        for(short body_i = 0; body_i < 2; ++body_i) {
            DissipatingBody &body = (body_i == 0 ? __body1 : __body2);
            for(unsigned zone_i = 0; zone_i < body.number_zones(); ++zone_i) {
                if(zone_i != 0)
                    os << "|";
                DissipatingZone &zone = body.zone(zone_i);
                os  << " ";
                if(
                        zone_i + body_i * __body1.number_zones()
                        ==
                        test_zone_ind
                ) {
                    if(lock_direction == 0)
                        os << "0";
                    else if(lock_direction > 0)
                        os << "+";
                    os  << lock_direction
                        << " (" << (*passed_i ? "V" : "X") << ")";
                    ++scenario_zone_i;
                    if(scenario_zone_i != lock_scenario.end())
                        std::tie(test_zone_ind,
                                 lock_direction,
                                 std::ignore) = *scenario_zone_i;
                } else {
                    if(zone.can_lock()) os << " NSZ  ";
                    else os << " NLZ  ";
                }
                os << " ";
            }
            os << "||";
        }
        os << std::endl;
    }
#endif

    bool BinarySystem::test_lock_scenario(
        const lock_scenario_type &lock_scenario
#ifndef NDEBUG
        , bool first_scenario
#endif
    )
    {
        set_lock_scenario(lock_scenario);
        double *above_lock_fraction_p = __above_lock_fractions[
            Dissipation::NO_DERIV
        ].data();
#ifndef DEBUG
        std::vector<bool> passed(lock_scenario.size());
        std::vector<bool>::iterator passed_i = passed.begin();
        bool result = true;
#endif
        for(
                lock_scenario_type::const_iterator
                    scenario_zone_iter = lock_scenario.begin();
                scenario_zone_iter != lock_scenario.end();
                ++scenario_zone_iter
        ) {
            unsigned test_zone_ind;
            short lock_direction;
            std::tie(test_zone_ind,
                     lock_direction,
                     std::ignore) = *scenario_zone_iter;
            if(lock_direction == 0) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Lock for zone " << test_zone_ind << " will ";
#endif
                if(
                        !(
                            *above_lock_fraction_p > 0
                            &&
                            *above_lock_fraction_p < 1
                        )
                ) {
#ifndef NDEBUG
                    *passed_i++ = false;
                    result = false;
#else
                   return false;
#endif
#ifdef VERBOSE_DEBUG
                    std::cerr << "not ";
#endif
                }
#ifndef NDEBUG
                else *passed_i++ = true;
#endif
#ifdef VERBOSE_DEBUG
                std::cerr << "hold (above lock frac: "
                          << *above_lock_fraction_p
                          << ")" << std::endl;
#endif
                ++above_lock_fraction_p;
            } else if(!test_synchronized_unlocked_zone(test_zone_ind)) {
#ifndef NDEBUG
                *passed_i++ = false;
                result = false;
#else
                return false;
#endif
            }
#ifndef NDEBUG
            else *passed_i++ = true;
#endif
        }
#ifndef NDEBUG
        describe_lock_scenario(std::cerr,
                               lock_scenario,
                               passed,
                               first_scenario || true);
        return result;
#else
        return true;
#endif
    }

    bool BinarySystem::explore_lock_scenarios(
        lock_scenario_type::const_iterator next_synchronized_zone,
        unsigned num_synchronized_zones,
        lock_scenario_type &lock_scenario
#ifndef NDEBUG
        , bool first_scenario
#endif

    )
    {
        if(lock_scenario.size() == num_synchronized_zones) {
#ifndef NDEBUG
            if(test_lock_scenario(lock_scenario, first_scenario)) {
                __selected_lock_scenario = lock_scenario;
                std::cerr << "Selected lock scenario: " << std::endl;
                describe_lock_scenario(
                    std::cerr,
                    lock_scenario,
                    std::vector<bool>(lock_scenario.size(), true),
                    true
                );
                describe_lock_scenario(
                    std::cerr,
                    __selected_lock_scenario,
                    std::vector<bool>(lock_scenario.size(), true),
                    true
                );

                return true;
            } else
                return false;
#else
            return test_lock_scenario(lock_scenario);
#endif
        }

        assert(lock_scenario.size() < num_synchronized_zones);

        lock_scenario.push_back(*next_synchronized_zone);
        short &scenario_dir = std::get<1>(lock_scenario.back());
        short orig_lock_direction = std::get<1>(*next_synchronized_zone);
#ifndef NDEBUG
        bool result = false;
#endif
        ++next_synchronized_zone;
        for(scenario_dir = -1; scenario_dir <= 1; ++scenario_dir) {
            if(
                    explore_lock_scenarios(
                        next_synchronized_zone,
                        num_synchronized_zones,
                        lock_scenario
#ifndef NDEBUG
                        , first_scenario
#endif
                    )
            ) {
#ifndef NDEBUG
                assert(!result);
                result = true;
#else
                return true;
#endif
            }
#ifndef NDEBUG
            first_scenario = false;
#endif
        }
        lock_scenario.pop_back();
#ifndef NDEBUG
        return result;
#else
        return false;
#endif
    }


    int BinarySystem::locked_surface_differential_equations(
        double *evolution_rates
    ) const
    {
        DissipatingZone &locked_zone = __body1.zone(0);
        locked_zone.set_evolution_rates(
            locked_zone.moment_of_inertia(1) * locked_zone.spin_frequency(),
            0.0,
            0.0
        );
        for(
            unsigned zone_index = 1;
            zone_index < __body1.number_zones();
            ++zone_index
        ) {
            evolution_rates[zone_index - 1] = __body1.nontidal_torque(
                zone_index
            )[2];

            __body1.zone(zone_index).set_evolution_rates(
                evolution_rates[zone_index - 1],
                0.0,
                0.0
            );

            assert(__body1.nontidal_torque(zone_index)[0] == 0);
            assert(__body1.nontidal_torque(zone_index)[1] == 0);

        }
        return 0;
    }

#ifdef ENABLE_DERIVATIVES
    void BinarySystem::locked_surface_jacobian(double *param_derivs,
                                               double *age_derivs) const
    {
        unsigned num_param = __body1.number_zones() - 1;
        double dR_dt = __body1.radius(1),
               dIabove_dt = __body1.zone(0).moment_of_inertia(1),
               dI_dt = __body1.zone(1).moment_of_inertia(1);
        for(unsigned row = 0; row < num_param; ++row) {
            double dIbelow_dt = 0;
            age_derivs[row] =
                __body1.nontidal_torque(row + 1, Dissipation::AGE)[2]
                +
                __body1.nontidal_torque(row + 1, Dissipation::RADIUS)[2] * dR_dt
                +
                __body1.nontidal_torque(row + 1,
                                        Dissipation::MOMENT_OF_INERTIA,
                                        0)[2] * dI_dt
                +
                __body1.nontidal_torque(row + 1,
                                        Dissipation::MOMENT_OF_INERTIA,
                                        -1)[2] * dIabove_dt;
            if(row < num_param - 1) {
                dIbelow_dt = __body1.zone(row + 2).moment_of_inertia(1);
                age_derivs[row] += __body1.nontidal_torque(
                    row + 1,
                    Dissipation::MOMENT_OF_INERTIA,
                    1)[2] * dIbelow_dt;
            }
            dIabove_dt = dI_dt;
            dI_dt = dIbelow_dt;

            assert(__body1.nontidal_torque(row + 1, Dissipation::AGE)[0] == 0);
            assert(__body1.nontidal_torque(row + 1, Dissipation::AGE)[1] == 0);
            assert(__body1.nontidal_torque(row + 1, Dissipation::RADIUS)[0]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1, Dissipation::RADIUS)[1]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1,
                                           Dissipation::MOMENT_OF_INERTIA,
                                           -1)[0]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1,
                                           Dissipation::MOMENT_OF_INERTIA,
                                           -1)[1]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1,
                                           Dissipation::MOMENT_OF_INERTIA,
                                           0)[0]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1,
                                           Dissipation::MOMENT_OF_INERTIA,
                                           0)[1]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1,
                                           Dissipation::MOMENT_OF_INERTIA,
                                           1)[0]
                   ==
                   0);
            assert(__body1.nontidal_torque(row + 1,
                                           Dissipation::MOMENT_OF_INERTIA,
                                           1)[1]
                   ==
                   0);

            for(unsigned col = 0; col < num_param; ++col) {
                double &dest = *(param_derivs + row * num_param + col);
                if(std::abs(static_cast<int>(row) - static_cast<int>(col)) < 1)
                    dest = 0;
                else {
                    dest = __body1.nontidal_torque(row + 1,
                                                   Dissipation::SPIN_ANGMOM,
                                                   col - row)[2];

                    assert(__body1.nontidal_torque(row + 1,
                                                   Dissipation::SPIN_ANGMOM,
                                                   col - row)[0]
                           ==
                           0);
                    assert(__body1.nontidal_torque(row + 1,
                                                   Dissipation::SPIN_ANGMOM,
                                                   col - row)[1]
                           ==
                           0);

                }
            }
        }
    }
#endif

    int BinarySystem::single_body_differential_equations(
        double *evolution_rates
    ) const
    {
        unsigned nzones = __body1.number_zones();
        double ref_angmom = __body1.zone(0).angular_momentum(),
               *inclination_evol = evolution_rates,
               *periapsis_evol = evolution_rates + (nzones - 1),
               *angmom_evol = periapsis_evol+(nzones - 1);

        Eigen::Vector3d reference_torque;
        for(unsigned zone_index = 0; zone_index < nzones; ++zone_index) {
            Eigen::Vector3d torque = __body1.nontidal_torque(zone_index);
            angmom_evol[zone_index] = torque[2];

            DissipatingZone &zone = __body1.zone(zone_index);
            if(zone_index) {
                zone.set_reference_zone_angmom(ref_angmom);
                inclination_evol[zone_index - 1] = zone.inclination_evolution(
                    zone_to_zone_transform(__body1.zone(0),
                                           zone,
                                           reference_torque),
                    torque
                );
                periapsis_evol[zone_index - 1] = zone.periapsis_evolution(
                    zone_to_zone_transform(__body1.zone(0),
                                           zone,
                                           reference_torque),
                    torque
                );
            } else reference_torque = torque;

            zone.set_evolution_rates(
                angmom_evol[zone_index],
                (zone_index ? inclination_evol[zone_index - 1] : 0.0),
                (zone_index ? periapsis_evol[zone_index -1] : 0.0)
            );
        }
#ifdef VERBOSE_DEBUG
        std::cerr << "rates: ";
        for(unsigned i = 0; i < 3 * nzones - 2; ++i) {
            if(i) std::cerr << ", ";
            std::cerr << evolution_rates[i];
        }
        std::cerr << std::endl;
#endif
        return 0;
    }

#ifdef ENABLE_DERIVATIVES
    void BinarySystem::fill_single_body_jacobian(
        double *inclination_param_derivs,
        double *periapsis_param_derivs,
        double *angmom_param_derivs,
        double *inclination_age_derivs,
        double *periapsis_age_derivs,
        double *angmom_age_derivs
    ) const
    {
        unsigned nzones=__body1.number_zones(), nparams = 3 * nzones - 2;
        Eigen::Vector3d ref_torque = __body1.nontidal_torque(0),
                        dref_torque_dincl = __body1.nontidal_torque(
                            0,
                            Dissipation::INCLINATION,
                            1
                        ),
                        dref_torque_dperi=__body1.nontidal_torque(
                            0,
                            Dissipation::PERIAPSIS,
                            1
                        ),
                        dref_torque_dangmom0=__body1.nontidal_torque(
                            0,
                            Dissipation::SPIN_ANGMOM,
                            0
                        ),
                        dref_torque_dangmom1=__body1.nontidal_torque(
                            0,
                            Dissipation::SPIN_ANGMOM,
                            1
                        ),
                        dref_torque_dage = __body1.nontidal_torque(
                            0,
                            Dissipation::AGE
                        ),
                        zero3d(0, 0, 0);
        for(
            unsigned deriv_zone_ind = 0;
            deriv_zone_ind < nzones;
            ++deriv_zone_ind
        ) {
            angmom_param_derivs[deriv_zone_ind - 1] =
                (deriv_zone_ind == 1 ? dref_torque_dincl[2] : 0);
            angmom_param_derivs[nzones+deriv_zone_ind - 2] =
                (deriv_zone_ind == 1 ? dref_torque_dperi[2] : 0);
            double &angmom_angmom_deriv =
                angmom_param_derivs[2 * nzones+deriv_zone_ind - 2];
            if(deriv_zone_ind==0) angmom_angmom_deriv = dref_torque_dangmom0[2];
            else if(deriv_zone_ind == 1)
                angmom_angmom_deriv = dref_torque_dangmom1[2];
            else angmom_angmom_deriv = 0;
        }
        angmom_age_derivs[0] = dref_torque_dage[2];
        DissipatingZone &ref_zone = __body1.zone(0);
        for(unsigned zone_ind = 1; zone_ind < nzones; ++zone_ind) {
            DissipatingZone &zone = __body1.zone(zone_ind);
            Eigen::Vector3d
                zone_ref_torque = zone_to_zone_transform(ref_zone,
                                                         zone,
                                                         ref_torque),
                zone_torque = __body1.nontidal_torque(zone_ind),
                ref_torque_age_deriv = zone_to_zone_transform(ref_zone,
                                                              zone,
                                                              dref_torque_dage),
                zone_torque_age_deriv=__body1.nontidal_torque(zone_ind,
                                                              Dissipation::AGE,
                                                              0);
            inclination_age_derivs[zone_ind - 1] =
                zone.inclination_evolution(zone_ref_torque,
                                           zone_torque,
                                           Dissipation::AGE,
                                           ref_torque_age_deriv,
                                           zone_torque_age_deriv);
            periapsis_age_derivs[zone_ind - 1] =
                zone.periapsis_evolution(zone_ref_torque,
                                         zone_torque,
                                         Dissipation::AGE,
                                         ref_torque_age_deriv,
                                         zone_torque_age_deriv);
            angmom_age_derivs[zone_ind] = zone_torque_age_deriv[2];
            for(
                unsigned quantity_ind = 0;
                quantity_ind < nparams;
                ++quantity_ind
            ) {
                unsigned offset = (zone_ind - 1) * nparams + quantity_ind;
                double &inclination_dest = inclination_param_derivs[offset],
                       &periapsis_dest = periapsis_param_derivs[offset],
                       &angmom_dest = angmom_param_derivs[offset+nparams];
                unsigned quantity_zone = quantity_ind;
                Dissipation::QuantityEntry with_respect_to;
                Eigen::Vector3d ref_torque_deriv(0, 0, 0);
                if(quantity_ind > 2 * nzones - 2) {
                    quantity_zone -= 2 * nzones - 2;
                    with_respect_to = Dissipation::SPIN_ANGMOM;
                    if(quantity_zone == 0)
                        ref_torque_deriv = dref_torque_dangmom0;
                    else if(quantity_zone == 1)
                        ref_torque_deriv = dref_torque_dangmom1;
                    if(quantity_zone <= 1)
                        ref_torque_deriv = zone_to_zone_transform(
                            ref_zone,
                            zone,
                            ref_torque_deriv
                        );
                } else {
                    if(quantity_ind >= nzones - 1) {
                        quantity_zone -= nzones - 2;
                        with_respect_to = Dissipation::PERIAPSIS;
                        if(quantity_zone == 1)
                            ref_torque_deriv = dref_torque_dperi;
                    } else {
                        quantity_zone += 1;
                        with_respect_to = Dissipation::INCLINATION;
                        if(quantity_zone == 1)
                            ref_torque_deriv = dref_torque_dincl;
                    }
                    if(quantity_zone == 1)
                        ref_torque_deriv = zone_to_zone_transform(
                            ref_zone,
                            zone,
                            ref_torque_deriv
                        );
                    ref_torque_deriv += zone_to_zone_transform(ref_zone,
                                                               zone,
                                                               ref_torque,
                                                               with_respect_to);
                }
                Eigen::Vector3d zone_torque_deriv;
                if(quantity_zone == zone_ind) {
                    zone_torque_deriv = __body1.nontidal_torque(zone_ind,
                                                                with_respect_to,
                                                                0);
                    inclination_dest=zone.inclination_evolution(
                        zone_ref_torque,
                        zone_torque,
                        with_respect_to,
                        ref_torque_deriv,
                        zone_torque_deriv
                    );
                    periapsis_dest=zone.periapsis_evolution(zone_ref_torque,
                                                            zone_torque,
                                                            with_respect_to,
                                                            ref_torque_deriv,
                                                            zone_torque_deriv);
                } else {
                    if(std::abs(static_cast<int>(quantity_zone-zone_ind)) == 1) {
                        zone_torque_deriv = __body1.nontidal_torque(
                            zone_ind,
                            with_respect_to,
                            quantity_zone - zone_ind
                        );
                    } else zone_torque_deriv = zero3d;
                    inclination_dest = zone.inclination_evolution(
                        ref_torque_deriv,
                        zone_torque_deriv
                    );
                    periapsis_dest = zone.periapsis_evolution(ref_torque_deriv,
                                                              zone_torque_deriv);
                }
                angmom_dest = zone_torque_deriv[2];
            }
        }
    }

    void BinarySystem::single_body_jacobian(double *param_derivs,
                                            double *age_derivs) const
    {
        unsigned nzones = __body1.number_zones(),
                 nparams=3*nzones-2;
        fill_single_body_jacobian(param_derivs,
                                  param_derivs + (nzones - 1) * nparams,
                                  param_derivs + 2 * (nzones - 1) * nparams,
                                  age_derivs,
                                  age_derivs + (nzones - 1),
                                  age_derivs + 2 * (nzones - 1));
    }
#endif

    double BinarySystem::semimajor_evolution(
        double orbit_power,
        double orbit_power_deriv
    ) const
    {
        if(std::isnan(orbit_power_deriv))
            return (orbit_power == 0
                    ? 0
                    : -__semimajor * orbit_power / __orbital_energy);

        else if(orbit_power == 0 && orbit_power_deriv == 0)
            return 0;

        return (
            -(2.0 * orbit_power
              +
              __semimajor * orbit_power_deriv
            )
            /
            __orbital_energy
        );
    }

    double BinarySystem::eccentricity_evolution(double orbit_power,
                                                double orbit_angmom_gain,
                                                double orbit_power_deriv,
                                                double orbit_angmom_gain_deriv,
                                                bool semimajor_deriv) const
    {
        if(__eccentricity <= 1e-8) return 0;
        double e2 = std::pow(__eccentricity, 2),
               factor = -(1.0 - e2) / (2.0 * __eccentricity);

        if(std::isnan(orbit_power_deriv))
            return factor * (orbit_power / __orbital_energy
                             +
                             2.0 * orbit_angmom_gain / __orbital_angmom);
        else if(semimajor_deriv)
            return (
                factor
                *
                (
                    (
                        orbit_power_deriv
                        +
                        orbit_power / __semimajor
                    )
                    /
                    __orbital_energy
                    +
                    (
                        2.0 * orbit_angmom_gain_deriv
                        -
                        orbit_angmom_gain / __semimajor
                    )
                    /
                    __orbital_angmom
                )
            );
        else return (
            factor
            *
            (
                orbit_power_deriv / __orbital_energy
                +
                2.0 * orbit_angmom_gain_deriv / __orbital_angmom
            )
            -
            2.0 * orbit_angmom_gain / __orbital_angmom
            -
            (1.0 + e2) / e2 * (
                orbit_power / __orbital_energy
                +
                2.0 * orbit_angmom_gain / __orbital_angmom
            )
        );
    }

    void BinarySystem::above_lock_problem_deriv_correction(
            Dissipation::QuantityEntry entry,
            bool body1_deriv,
            Eigen::MatrixXd &matrix,
            Eigen::VectorXd &rhs
    ) const
    {
        unsigned deriv_zone_index = (body1_deriv
                                     ? 0
                                     : __body1.number_locked_zones());
        DissipatingBody &deriv_body = (body1_deriv ? __body1 : __body2);
        DissipatingZone &deriv_zone = deriv_body.zone(0);

        if(
            entry == Dissipation::ORBITAL_FREQUENCY
            ||
            entry == Dissipation::SEMIMAJOR
        ) {
            double coef = (entry == Dissipation::SEMIMAJOR
                           ? 1
                           : Core::orbital_angular_velocity(__body1.mass(),
                                                            __body2.mass(),
                                                            __semimajor,
                                                            true));
            bool done = false;
            unsigned locked_ind = 0;
            for(DissipatingBody *body = &__body1; !done; body = &__body2) {
                for(
                    unsigned zone_ind = 0;
                    zone_ind < body->number_zones();
                    ++zone_ind
                ) {
                    rhs.array() += (1.5
                                    *
                                    (
                                        __body1.tidal_orbit_power()
                                        +
                                        __body2.tidal_orbit_power()
                                    )
                                    /
                                    __orbital_energy
                                    /
                                    __semimajor
                                    *
                                    coef);
                    if(body->zone(zone_ind).locked()) {
                        matrix.col(locked_ind).array() += (
                            1.5
                            *
                            (
                                body->tidal_power(zone_ind, true)
                                -
                                body->tidal_power(zone_ind, false)
                            )
                            /
                            __orbital_energy
                            /
                            __semimajor * coef
                        );

                        ++locked_ind;
                    }
                }
                done = (body == &__body2);
            }
        } else if(deriv_zone.locked()) {
            if(
                entry == Dissipation::SPIN_FREQUENCY
                ||
                entry == Dissipation::SPIN_ANGMOM
            ) {
                double coef = (
                    entry == Dissipation::SPIN_FREQUENCY
                    ? deriv_zone.moment_of_inertia()
                    : 1
                ) / std::pow(deriv_zone.angular_momentum(), 2);
                matrix(deriv_zone_index, deriv_zone_index) -=(
                    deriv_body.tidal_torque(0, true)[2]
                    -
                    deriv_body.tidal_torque(0, false)[2]
                ) * coef;
                rhs(deriv_zone_index) -= (
                    deriv_body.tidal_torque(0, false)[2]
                    *
                    coef
                );
            } else if(entry == Dissipation::MOMENT_OF_INERTIA) {
                rhs(deriv_zone_index) -= (
                    deriv_zone.moment_of_inertia(1)
                    /
                    std::pow(deriv_zone.moment_of_inertia(), 2)
                );
            }
        }
    }

    void BinarySystem::calculate_above_lock_fractions(
        Eigen::VectorXd &fractions,
        Dissipation::QuantityEntry entry,
        bool body1_deriv
    )
    {
        unsigned num_locked_zones = __locked_zones.size();

        assert(num_locked_zones == (__body1.number_locked_zones()
                                    +
                                    __body2.number_locked_zones()));
#ifndef NDEBUG
        std::cerr << "Setting "
                  << num_locked_zones
                  << " above lock fractions "
                  << entry
                  << std::endl;
#endif

        if(num_locked_zones == 0) {
            fractions.resize(0);
            return;
        }
        std::valarray<double> nontidal_torque(num_locked_zones),
                              tidal_torque_z_above(num_locked_zones),
                              tidal_torque_z_below(num_locked_zones),
                              tidal_power_difference(num_locked_zones);
        unsigned locked_zone_ind = 0;
        for(
            std::list<unsigned>::const_iterator zi = __locked_zones.begin();
            zi != __locked_zones.end();
            ++zi
        ) {
            DissipatingBody &body = (
                locked_zone_ind < __body1.number_locked_zones()
                ? __body1
                : __body2
            );
            nontidal_torque[locked_zone_ind] =
                body.nontidal_torque(*zi, entry)[2];
            tidal_torque_z_above[locked_zone_ind] =
                body.tidal_torque(*zi, true, entry)[2];
            tidal_torque_z_below[locked_zone_ind] =
                body.tidal_torque(*zi, false, entry)[2];
            if(!zone_specific(entry) || *zi == 0) {
                tidal_power_difference[locked_zone_ind] = (
                    body.tidal_power(*zi, true, entry)
                    -
                    body.tidal_power(*zi, false, entry)
                );
            }
            ++locked_zone_ind;
        }
        Eigen::MatrixXd matrix(num_locked_zones, num_locked_zones);
        Eigen::VectorXd rhs(num_locked_zones);
        rhs.setConstant(1.5
                        *
                        (
                            __body1.tidal_orbit_power(entry)
                            +
                            __body2.tidal_orbit_power(entry)
                        )
                        /
                        __orbital_energy);
        unsigned i = 0;
        for(
            std::list<unsigned>::const_iterator
            zi=__locked_zones.begin();
            zi!=__locked_zones.end();
            ++zi
        ) {
            if(!zone_specific(entry) || *zi == 0)
                matrix.col(i).setConstant(1.5
                                          *
                                          tidal_power_difference[i]
                                          /
                                          __orbital_energy);
            else matrix.col(i).setZero();
            DissipatingZone &zone = (i < __body1.number_locked_zones()
                                     ? __body1
                                     : __body2).zone(*zi);
            matrix(i, i) += (
                (tidal_torque_z_above[i] - tidal_torque_z_below[i])
                /
                zone.angular_momentum()
            );
            if(entry == Dissipation::NO_DERIV)
                rhs(i) += (zone.moment_of_inertia(1)
                           /
                           zone.moment_of_inertia());
            else if(entry == Dissipation::AGE)
                rhs(i) += (zone.moment_of_inertia(2)
                           /
                           zone.moment_of_inertia());
            rhs(i) -= ((nontidal_torque[i] + tidal_torque_z_below[i])
                       /
                       zone.angular_momentum());
            ++i;
        }
        above_lock_problem_deriv_correction(entry, body1_deriv, matrix, rhs);
        if(entry == Dissipation::NO_DERIV) {
            __above_lock_fractions_decomp.compute(matrix);
            fractions = __above_lock_fractions_decomp.solve(rhs);
        } else {
            fractions = __above_lock_fractions_decomp.solve(
                rhs - matrix * __above_lock_fractions[Dissipation::NO_DERIV]
            );
        }
    }

    Eigen::VectorXd BinarySystem::above_lock_fractions_deriv(
            Dissipation::QuantityEntry entry,
            DissipatingBody &body,
            unsigned zone_index)
    {
        assert(entry == Dissipation::INCLINATION
               ||
               entry == Dissipation::PERIAPSIS
               ||
               entry == Dissipation::MOMENT_OF_INERTIA
               ||
               entry == Dissipation::SPIN_ANGMOM);
        assert(zone_index > 0 || &body == &__body2);
        assert(body.number_zones() > zone_index);
        assert(number_locked_zones()
               ==
               __above_lock_fractions[Dissipation::NO_DERIV].size());

        unsigned num_locked_zones = number_locked_zones();
        if(num_locked_zones == 0) return Eigen::VectorXd();
        DissipatingZone &deriv_zone = body.zone(zone_index);
        unsigned locked_zone_index=(&body == &__body1
                                    ? 0
                                    : __body1.number_locked_zones());
        Eigen::VectorXd rhs;
        for(unsigned i = 0; i < zone_index; ++i)
            if(body.zone(i).locked()) ++locked_zone_index;
        if(deriv_zone.locked()) {
            double above_frac =
                __above_lock_fractions[Dissipation::NO_DERIV][locked_zone_index];
            rhs.setConstant(
                num_locked_zones,
                -1.5 * (above_frac*body.tidal_power(zone_index, true, entry)
                        +
                        (1.0 - above_frac) * body.tidal_power(zone_index,
                                                              false,
                                                              entry))
                -
                body.nontidal_torque(zone_index, entry)[2]
            );
            rhs(locked_zone_index) -=
                above_frac*body.tidal_torque(zone_index, true, entry)[2]
                +
                (1.0 - above_frac) * body.tidal_torque(zone_index,
                                                       false,
                                                       entry)[2];
            if(entry == Dissipation::MOMENT_OF_INERTIA)
                rhs(locked_zone_index) -= (
                    deriv_zone.moment_of_inertia(1)
                    /
                    std::pow(deriv_zone.moment_of_inertia(), 2)
                );
            else if(entry == Dissipation::SPIN_ANGMOM)
                rhs(locked_zone_index) += (
                    above_frac*body.tidal_torque(zone_index, true)[2]
                    +
                    (1.0 - above_frac)
                    *
                    body.tidal_torque(zone_index, false)[2]
                    +
                    body.nontidal_torque(zone_index)[2]
                ) / std::pow(deriv_zone.angular_momentum(), 2);
        } else {
            rhs.setConstant(num_locked_zones,
                            -1.5 * body.tidal_power(zone_index, false, entry));
        }
        if(zone_index > 0 && body.zone(zone_index - 1).locked())
            rhs(body.zone(zone_index - 1).locked_zone_index()) -=
                body.nontidal_torque(zone_index - 1, entry, 1)[2];
        if(
            zone_index + 1 < body.number_zones()
            &&
            body.zone(zone_index + 1).locked()
        )
            rhs(body.zone(zone_index + 1).locked_zone_index()) -=
                body.nontidal_torque(zone_index + 1, entry, -1)[2];
        return __above_lock_fractions_decomp.solve(rhs);
    }

    void BinarySystem::fill_above_lock_fractions_deriv()
    {
        unsigned num_zones = __body1.number_zones() + __body2.number_zones();
        DissipatingBody *body = &__body1;
        __above_lock_fractions_inclination_deriv.resize(num_zones);
        __above_lock_fractions_periapsis_deriv.resize(num_zones);
        __above_lock_fractions_inertia_deriv.resize(num_zones);
        __above_lock_fractions_angmom_deriv.resize(num_zones);
        __above_lock_fractions_inclination_deriv[0] =
            __above_lock_fractions[Dissipation::INCLINATION];
        __above_lock_fractions_periapsis_deriv[0] =
            __above_lock_fractions[Dissipation::PERIAPSIS];
        __above_lock_fractions_inertia_deriv[0] =
            __above_lock_fractions[Dissipation::MOMENT_OF_INERTIA];
        __above_lock_fractions_angmom_deriv[0] =
            __above_lock_fractions[Dissipation::SPIN_ANGMOM];
        unsigned body_zone_ind = 0;
        for(unsigned zone_ind = 1; zone_ind < num_zones; ++zone_ind) {
            ++body_zone_ind;
            if(body_zone_ind == body->number_zones()) {
                body_zone_ind = 0;
                body = &__body2;
            }
            __above_lock_fractions_inclination_deriv[zone_ind] =
                above_lock_fractions_deriv(Dissipation::INCLINATION,
                                           *body,
                                           body_zone_ind);
            __above_lock_fractions_periapsis_deriv[zone_ind] =
                above_lock_fractions_deriv(Dissipation::PERIAPSIS,
                                           *body,
                                           body_zone_ind);
            __above_lock_fractions_inertia_deriv[zone_ind] =
                above_lock_fractions_deriv(Dissipation::MOMENT_OF_INERTIA,
                                           *body,
                                           body_zone_ind);
            __above_lock_fractions_angmom_deriv[zone_ind] =
                above_lock_fractions_deriv(Dissipation::SPIN_ANGMOM,
                                           *body,
                                           body_zone_ind);
        }
    }

    void BinarySystem::update_above_lock_fractions()
    {
        std::valarray<Eigen::VectorXd>
            body2_above_lock_fractions(Dissipation::NUM_DERIVATIVES);
        for(
            int deriv = Dissipation::NO_DERIV;
            deriv < Dissipation::NUM_DERIVATIVES;
            ++deriv
        ) {
            calculate_above_lock_fractions(
                __above_lock_fractions[deriv],
                static_cast<Dissipation::QuantityEntry>(deriv)
            );
            if(!zone_specific(static_cast < Dissipation::QuantityEntry>(deriv)))
                body2_above_lock_fractions[deriv] = __above_lock_fractions[deriv];
            else calculate_above_lock_fractions(
                body2_above_lock_fractions[deriv],
                static_cast<Dissipation::QuantityEntry>(deriv),
                false
            );
        }
        __above_lock_fractions_body2_radius_deriv =
            body2_above_lock_fractions[Dissipation::RADIUS];
        __body1.set_above_lock_fractions(__above_lock_fractions);
        __body2.set_above_lock_fractions(body2_above_lock_fractions);
        fill_above_lock_fractions_deriv();
    }

    void BinarySystem::fill_orbit_torque_and_power()
    {
        __orbit_torque = (
            __body1.tidal_orbit_torque()
            +
            zone_to_zone_transform(__body2.zone(0),
                                   __body1.zone(0),
                                   __body2.tidal_orbit_torque())
        );

        __orbit_power = (__body1.tidal_orbit_power()
                         +
                         __body2.tidal_orbit_power());

        __orbit_angmom_gain = (__orbit_torque[0]
                               *
                               std::sin(__body1.zone(0).inclination())
                               +
                               __orbit_torque[2]
                               *
                               std::cos(__body1.zone(0).inclination()));
    }

    int BinarySystem::binary_differential_equations(
        double *differential_equations
    ) const
    {
        DissipatingZone &reference_zone = __body1.zone(0);
        unsigned num_body1_zones = __body1.number_zones(),
                 num_total_zones = num_body1_zones + __body2.number_zones();
        double *inclination_rates = differential_equations + 2,
               *periapsis_rates = inclination_rates + num_total_zones,
               *angmom_rates = periapsis_rates + num_total_zones - 1;
        unsigned angmom_skipped = 0;
        double reference_periapsis_rate = Core::NaN;
        for(unsigned zone_ind = 0; zone_ind < num_total_zones; ++zone_ind) {
            DissipatingBody &body = (zone_ind < num_body1_zones
                                     ? __body1
                                     : __body2);
            unsigned body_zone_ind = zone_ind;
            if(zone_ind >= num_body1_zones)
                body_zone_ind -= num_body1_zones;
            DissipatingZone &zone = body.zone(body_zone_ind);
            Eigen::Vector3d zone_orbit_torque,
                            total_zone_torque;
            total_zone_torque = body.nontidal_torque(body_zone_ind);

            zone_orbit_torque = (zone_ind
                                 ? zone_to_zone_transform(reference_zone,
                                                          zone,
                                                          __orbit_torque)
                                 : __orbit_torque);

            Dissipation::QuantityEntry entry = Dissipation::NO_DERIV;
            if(zone.locked()) {
                double above_frac = (__above_lock_fractions
                                     [Dissipation::NO_DERIV]
                                     [angmom_skipped]);
                total_zone_torque += (
                    above_frac
                    *
                    body.tidal_torque(body_zone_ind, true, entry)
                    +
                    (1.0 - above_frac)
                    *
                    body.tidal_torque(body_zone_ind, false, entry)
                );
                ++angmom_skipped;
            } else total_zone_torque += body.tidal_torque(body_zone_ind,
                                                          false,
                                                          entry);
            inclination_rates[zone_ind] = zone.inclination_evolution(
                zone_orbit_torque,
                total_zone_torque,
                entry
            );
            assert(!std::isnan(inclination_rates[zone_ind]));
            if(zone_ind) {
                periapsis_rates[zone_ind - 1] = zone.periapsis_evolution(
                    zone_orbit_torque,
                    total_zone_torque,
                    entry
                ) - reference_periapsis_rate;
                assert(!std::isnan(periapsis_rates[zone_ind - 1]));
            } else {
                reference_periapsis_rate = zone.periapsis_evolution(
                    zone_orbit_torque,
                    total_zone_torque,
                    entry
                );
                assert(!std::isnan(reference_periapsis_rate));

            }
            if(!zone.locked()) {
                angmom_rates[zone_ind - angmom_skipped] = total_zone_torque[2];
                assert(!std::isnan(angmom_rates[zone_ind - angmom_skipped]));
            }

            zone.set_evolution_rates(
                total_zone_torque[2],
                inclination_rates[zone_ind],
                (zone_ind ? periapsis_rates[zone_ind - 1] : 0.0)
            );

#ifdef VERBOSE_DEBUG
            std::cerr << "Zone " << zone_ind
                      << " torque: " << total_zone_torque
                      << " inclination rate";
            std::cerr << ": " << inclination_rates[zone_ind];
            if(zone_ind)
                std::cerr << " periapsis rate: "
                          << periapsis_rates[zone_ind - 1];

            std::cerr << std::endl;
#endif
            assert(!std::isnan(total_zone_torque.sum()));

        }

        differential_equations[0] = semimajor_evolution(__orbit_power);
        differential_equations[1] = eccentricity_evolution(__orbit_power,
                                                           __orbit_angmom_gain);
        __semimajor_rate = differential_equations[0];
        __eccentricity_rate = differential_equations[1];

        if(!angmom_skipped)
            differential_equations[0] *= 6.5 * std::pow(__semimajor, 5.5);

#ifdef VERBOSE_DEBUG
        std::cerr << "rates: ";
        for(
            unsigned i = 0;
            i < 3 * (__body1.number_zones() + __body2.number_zones()) - 1;
            ++i
        ) {
            if(i) std::cerr << ", ";
            std::cerr << differential_equations[i];
        }
        std::cerr << std::endl;
#endif

        return 0;
    }

    template<typename VALUE_TYPE>
    void BinarySystem::add_body_rate_deriv(
        const DissipatingBody &body,
        VALUE_TYPE (DissipatingBody::*func)(Dissipation::QuantityEntry,
                                            unsigned,
                                            const Eigen::VectorXd &) const,
        std::valarray<VALUE_TYPE> &orbit_rate_deriv,
        unsigned offset
    ) const
    {
        unsigned num_zones = __body1.number_zones() + __body2.number_zones(),
                 num_param = orbit_rate_deriv.size() - 1;
        orbit_rate_deriv[0] += (body.*func)(Dissipation::SEMIMAJOR,
                                            0,
                                            Eigen::VectorXd());
        orbit_rate_deriv[1] += (body.*func)(Dissipation::ECCENTRICITY,
                                            0,
                                            Eigen::VectorXd());
        orbit_rate_deriv[num_param] += (
            (body.*func)(Dissipation::AGE, 0, Eigen::VectorXd())
            +
            (body.*func)(Dissipation::RADIUS, 0, Eigen::VectorXd())
            *
            body.radius(1)
        );
        unsigned locked_zone_count = (offset ? __body1.number_locked_zones() : 0);
        for(unsigned zone_ind = 0; zone_ind < body.number_zones(); ++zone_ind) {
            orbit_rate_deriv[zone_ind+2+offset] =
                (body.*func)(
                    Dissipation::INCLINATION,
                    zone_ind,
                    __above_lock_fractions_inclination_deriv[zone_ind + offset]
                );
            if(zone_ind+offset > 0)
                orbit_rate_deriv[zone_ind + 1 + num_zones + offset] =
                    (body.*func)(
                        Dissipation::PERIAPSIS,
                        zone_ind,
                        __above_lock_fractions_periapsis_deriv[zone_ind+offset]
                    );
            const DissipatingZone &zone = body.zone(zone_ind);
            if(zone.locked()) ++locked_zone_count;
            else orbit_rate_deriv[
                zone_ind + 1 + 2 * num_zones - locked_zone_count
            ] = (body.*func)(
                Dissipation::SPIN_ANGMOM,
                zone_ind,
                __above_lock_fractions_angmom_deriv[zone_ind + offset]
            );
            orbit_rate_deriv[num_param] += (
                (body.*func)(
                    Dissipation::MOMENT_OF_INERTIA,
                    zone_ind,
                    __above_lock_fractions_inertia_deriv[zone_ind + offset]
                )
                *
                zone.moment_of_inertia(1)
            );
        }
    }

    void BinarySystem::fill_orbit_power_deriv(
            std::valarray<double> &orbit_power_deriv) const
    {
        orbit_power_deriv[0] = 0;
        orbit_power_deriv[1] = 0;
        orbit_power_deriv[orbit_power_deriv.size() - 1] = 0;
        add_body_rate_deriv(__body1, &DissipatingBody::tidal_orbit_power,
                            orbit_power_deriv,
                            0);
        add_body_rate_deriv(__body2,
                            &DissipatingBody::tidal_orbit_power,
                            orbit_power_deriv,
                            __body1.number_zones());
    }

    void BinarySystem::fill_orbit_angmom_gain_deriv(
            std::valarray<double> &orbit_angmom_gain_deriv
    ) const
    {
        std::valarray<Eigen::Vector3d>
            body1_orbit_torque_deriv(orbit_angmom_gain_deriv.size()),
            body2_orbit_torque_deriv(orbit_angmom_gain_deriv.size());
        body1_orbit_torque_deriv[0] =
            body1_orbit_torque_deriv[1] =
            body2_orbit_torque_deriv[0] =
            body2_orbit_torque_deriv[1] = Eigen::Vector3d(0, 0, 0);
        add_body_rate_deriv(__body1,
                            &DissipatingBody::tidal_orbit_torque,
                            body1_orbit_torque_deriv,
                            0);
        add_body_rate_deriv(__body2, &DissipatingBody::tidal_orbit_torque,
                            body2_orbit_torque_deriv,
                            __body1.number_zones());
        double body1_sin_inc = std::sin(__body1.zone(0).inclination()),
               body1_cos_inc = std::sin(__body1.zone(0).inclination()),
               body2_sin_inc = std::sin(__body2.zone(0).inclination()),
               body2_cos_inc = std::sin(__body2.zone(0).inclination());
        for(unsigned i = 0; i < orbit_angmom_gain_deriv.size(); ++i)
            orbit_angmom_gain_deriv[i] = (
                body1_sin_inc * body1_orbit_torque_deriv[i][0]
                +
                body1_cos_inc * body1_orbit_torque_deriv[i][2]
                +
                body2_sin_inc * body2_orbit_torque_deriv[i][0]
                +
                body2_cos_inc * body2_orbit_torque_deriv[i][2]
            );
        Eigen::Vector3d body1_orbit_torque = __body1.tidal_orbit_torque(),
                        body2_orbit_torque = __body2.tidal_orbit_torque();
        orbit_angmom_gain_deriv[2] +=(body1_cos_inc * body1_orbit_torque[0]
                                      -
                                      body1_sin_inc * body1_orbit_torque[2]
                                      +
                                      body2_cos_inc * body2_orbit_torque[0]
                                      -
                                      body2_sin_inc * body2_orbit_torque[2]);
    }

#ifdef ENABLE_DERIVATIVES
    void BinarySystem::semimajor_jacobian(
            const std::valarray<double> &orbit_power_deriv,
            bool a6p5,
            double *param_derivs,
            double &age_deriv
    ) const
    {
        param_derivs[0] = semimajor_evolution(__orbit_power,
                                              orbit_power_deriv[0]);
        if(a6p5) param_derivs[0] +=
            5.5 * semimajor_evolution(__orbit_power) / __semimajor;
        unsigned i = 1;
        for(; i < orbit_power_deriv.size() - 1; ++i)
            param_derivs[i] = semimajor_evolution(orbit_power_deriv[i]);
        age_deriv = semimajor_evolution(orbit_power_deriv[i]);
    }

    void BinarySystem::eccentricity_jacobian(
        const std::valarray<double> &orbit_power_deriv,
        const std::valarray<double> &orbit_angmom_gain_deriv,
        bool a6p5,
        double *param_derivs,
        double &age_deriv
    ) const
    {
        param_derivs[0] = eccentricity_evolution(__orbit_power,
                                                 __orbit_angmom_gain,
                                                 orbit_power_deriv[0],
                                                 orbit_angmom_gain_deriv[0],
                                                 true);
        if(a6p5) param_derivs[0] /= 6.5 * std::pow(__semimajor, 5.5);
        param_derivs[1] = eccentricity_evolution(__orbit_power,
                                                 __orbit_angmom_gain,
                                                 orbit_power_deriv[1],
                                                 orbit_angmom_gain_deriv[1],
                                                 false);
        unsigned i = 2;
        for(; i < orbit_power_deriv.size() - 1; ++i)
            param_derivs[i] = eccentricity_evolution(orbit_power_deriv[i],
                                                     orbit_angmom_gain_deriv[i]);
        age_deriv=eccentricity_evolution(orbit_power_deriv[i],
                                         orbit_angmom_gain_deriv[i]);
    }

    void BinarySystem::angle_evolution_age_deriv(DissipatingBody &body,
                                                 unsigned zone_ind,
                                                 double sin_inc,
                                                 double cos_inc,
                                                 unsigned locked_zone_ind,
                                                 double &inclination,
                                                 double &periapsis) const
    {
        double dR1_dt = __body1.radius(1),
               dR2_dt = __body2.radius(1);
        DissipatingZone &zone = body.zone(zone_ind);
        Eigen::Vector3d orbit_torque_age_deriv = (
            __body1.tidal_orbit_torque(zone, Dissipation::AGE)
            +
            __body2.tidal_orbit_torque(zone, Dissipation::AGE)
            +
            __body1.tidal_orbit_torque(zone, Dissipation::RADIUS) * dR1_dt
            +
            __body2.tidal_orbit_torque(zone, Dissipation::RADIUS) * dR2_dt
        );
        DissipatingBody *other_body = &__body1;
        unsigned body_zone_ind = 0,
                 num_zones = number_zones();
        for(
            unsigned other_zone_ind = 0;
            other_zone_ind<num_zones;
            ++other_zone_ind
        ) {
            orbit_torque_age_deriv += (
                other_body->tidal_orbit_torque(
                    zone,
                    Dissipation::MOMENT_OF_INERTIA,
                    body_zone_ind,
                    __above_lock_fractions_inertia_deriv[other_zone_ind]
                )
                *
                other_body->zone(body_zone_ind).moment_of_inertia(1)
            );
            ++body_zone_ind;
            if(body_zone_ind == other_body->number_zones()) {
                body_zone_ind = 0;
                other_body = &__body2;
            }
        }
        periapsis = (orbit_torque_age_deriv[1]
                     *
                     cos_inc
                     /
                     (__orbital_angmom * sin_inc));
        inclination = (orbit_torque_age_deriv[0] * cos_inc
                       -
                       orbit_torque_age_deriv[2] * sin_inc) / __orbital_angmom;
        Eigen::Vector3d to_add = -body.nontidal_torque(zone_ind,
                                                       Dissipation::AGE);
        if(zone.locked()) {
            double above_frac =
                __above_lock_fractions[Dissipation::NO_DERIV][locked_zone_ind];
            double above_frac_deriv =
                __above_lock_fractions[Dissipation::AGE][locked_zone_ind];
            to_add -= (
                body.tidal_torque(zone_ind, true, Dissipation::AGE)
                *
                above_frac
                +
                body.tidal_torque(zone_ind, false, Dissipation::AGE)
                *
                (1.0 - above_frac)
                +
                above_frac_deriv*(body.tidal_torque(zone_ind, true)
                                  -
                                  body.tidal_torque(zone_ind, false))
            );
        } else to_add += body.tidal_torque(zone_ind, false, Dissipation::AGE);
        to_add /= zone.angular_momentum();
        inclination += to_add[0];
        periapsis -= to_add[1] / sin_inc;
    }

    void BinarySystem::angle_evolution_orbit_deriv(Dissipation::QuantityEntry entry,
                                                   double angmom_deriv,
                                                   DissipatingBody &body,
                                                   unsigned zone_ind,
                                                   double sin_inc,
                                                   double cos_inc,
                                                   unsigned locked_zone_ind,
                                                   double &inclination,
                                                   double &periapsis) const
    {
        assert(entry == Dissipation::SEMIMAJOR
               ||
               entry == Dissipation::ECCENTRICITY);

        DissipatingZone &zone = body.zone(zone_ind);
        Eigen::Vector3d orbit_torque = (__body1.tidal_orbit_torque(zone)
                                        +
                                        __body2.tidal_orbit_torque(zone)),
                        orbit_torque_deriv = (
                            __body1.tidal_orbit_torque(zone, entry)
                            +
                            __body2.tidal_orbit_torque(zone, entry)
                        );
        inclination = (orbit_torque_deriv[0] * cos_inc
                       -
                       orbit_torque_deriv[2] * sin_inc
                       -
                       (orbit_torque[0] * cos_inc - orbit_torque[2] * sin_inc)
                       *
                       angmom_deriv
                       /
                       __orbital_angmom) / __orbital_angmom;
        periapsis = (
            (
                orbit_torque[1] * angmom_deriv / __orbital_angmom
                -
                orbit_torque_deriv[1]
            )
            *
            cos_inc
            /
            (__orbital_angmom * sin_inc)
        );
        Eigen::Vector3d to_add =- body.nontidal_torque(zone_ind, entry);
        if(zone.locked()) {
            double
                above_frac =
                __above_lock_fractions[Dissipation::NO_DERIV][locked_zone_ind],
                above_frac_deriv = __above_lock_fractions[entry][locked_zone_ind];
            to_add -= (
                body.tidal_torque(zone_ind, true, entry) * above_frac
                +
                body.tidal_torque(zone_ind, false, entry)
                *
                (1.0 - above_frac)
                +
                above_frac_deriv
                *
                (
                    body.tidal_torque(zone_ind, true)
                    -
                    body.tidal_torque(zone_ind, false)
                )
            );
        } else to_add -= body.tidal_torque(zone_ind, false, entry);
        to_add /= zone.angular_momentum();
        inclination += to_add[0];
        periapsis -= to_add[1] / sin_inc;
    }

    void BinarySystem::fill_orbit_torque_deriv(
        Dissipation::QuantityEntry entry,
        DissipatingBody &body,
        unsigned zone_ind,
        std::valarray<Eigen::Vector3d> &orbit_torque_deriv
    ) const
    {
        assert(entry == Dissipation::INCLINATION
               ||
               entry == Dissipation::PERIAPSIS
               ||
               entry == Dissipation::MOMENT_OF_INERTIA
               ||
               entry == Dissipation::SPIN_ANGMOM);
        assert(orbit_torque_deriv.size()
               ==
                __body1.number_zones() + __body2.number_zones());

        unsigned body_deriv_zone_ind = 0,
                 num_zones = orbit_torque_deriv.size();
        DissipatingBody *deriv_body = &__body1;
        const std::valarray<Eigen::VectorXd> *above_frac_deriv;
        if(entry == Dissipation::INCLINATION)
            above_frac_deriv = &__above_lock_fractions_inclination_deriv;
        else if(entry == Dissipation::PERIAPSIS)
            above_frac_deriv = &__above_lock_fractions_periapsis_deriv;
        else if(entry == Dissipation::MOMENT_OF_INERTIA)
            above_frac_deriv = &__above_lock_fractions_inertia_deriv;
        else
            above_frac_deriv = &__above_lock_fractions_angmom_deriv;

        DissipatingZone &zone = body.zone(zone_ind);
        for(
            unsigned deriv_zone_ind = 0;
            deriv_zone_ind < num_zones;
            ++deriv_zone_ind
        ) {
            orbit_torque_deriv[deriv_zone_ind] =
                deriv_body->tidal_orbit_torque(
                    zone,
                    entry,
                    body_deriv_zone_ind,
                    (*above_frac_deriv)[deriv_zone_ind]
                );
            if(
                (
                    entry == Dissipation::INCLINATION
                    ||
                    entry == Dissipation::PERIAPSIS
                )
                &&
                deriv_body == &body
                &&
                zone_ind == body_deriv_zone_ind
            ) {
                DissipatingBody &other_body = (&body == &__body1
                                               ? __body2
                                               : __body1);
                orbit_torque_deriv[deriv_zone_ind] +=
                    zone_to_zone_transform(other_body.zone(0),
                                           zone,
                                           other_body.tidal_orbit_torque(),
                                           entry,
                                           false);
            }
            ++body_deriv_zone_ind;
            if(body_deriv_zone_ind == deriv_body->number_zones()) {
                body_deriv_zone_ind = 0;
                deriv_body = &__body2;
            }
        }
    }

    void BinarySystem::fill_zone_torque_deriv(
        Dissipation::QuantityEntry entry,
        DissipatingBody &body,
        unsigned zone_ind,
        std::valarray<Eigen::Vector3d> &zone_torque_deriv
    ) const
    {
#ifndef NDEBUG
        if(body.zone(zone_ind).locked())
            assert(zone_torque_deriv.size() == 4);
        else
            assert(zone_torque_deriv.size() == 3);
#endif

        zone_torque_deriv[0] = (zone_ind == 0
                                ? Eigen::Vector3d::Zero()
                                : body.nontidal_torque(zone_ind, entry, -1));
        Eigen::Vector3d nontidal_torque = body.nontidal_torque(zone_ind,
                                                               entry,
                                                               0);
        zone_torque_deriv[1] = (nontidal_torque
                                +
                                body.tidal_torque(zone_ind, false, entry));
        zone_torque_deriv[2] = (zone_ind < body.number_zones() - 1
                                ? body.nontidal_torque(zone_ind, entry, 1)
                                : Eigen::Vector3d::Zero());
        if(body.zone(zone_ind).locked())
            zone_torque_deriv[3] = (nontidal_torque
                                    +
                                    body.tidal_torque(zone_ind, true, entry));
    }

    void BinarySystem::inclination_evolution_zone_derivs(
            Dissipation::QuantityEntry entry,
            DissipatingBody &body,
            unsigned zone_ind,
            double zone_x_torque_above,
            double zone_x_torque_below,
            const std::valarray<Eigen::Vector3d> &zone_torque_deriv,
            const Eigen::Vector3d &orbit_torque,
            const std::valarray<Eigen::Vector3d> &orbit_torque_deriv,
            const std::valarray<Eigen::VectorXd> &above_frac_deriv,
            double sin_inc,
            double cos_inc,
            unsigned locked_zone_ind,
            double *result
    ) const
    {
        assert(entry == Dissipation::INCLINATION
               ||
               entry == Dissipation::PERIAPSIS
               ||
               entry == Dissipation::SPIN_ANGMOM);
        assert(orbit_torque_deriv.size()
               ==
               __body1.number_zones() + __body2.number_zones());

        DissipatingZone &zone = body.zone(zone_ind);
        unsigned num_zones = orbit_torque_deriv.size();
        unsigned global_zone_ind = zone_ind;
        if(&body == &__body2) global_zone_ind += __body1.number_zones();

        double above_frac =
            __above_lock_fractions[Dissipation::NO_DERIV][locked_zone_ind];

        for(
            unsigned deriv_zone_ind = (entry == Dissipation::PERIAPSIS ? 1 : 0);
            deriv_zone_ind < num_zones;
            ++deriv_zone_ind
        ) {
            double &dest = (entry == Dissipation::PERIAPSIS
                            ? result[deriv_zone_ind - 1]
                            : result[deriv_zone_ind]);
            dest = (
                orbit_torque_deriv[deriv_zone_ind][0] * cos_inc
                -
                orbit_torque_deriv[deriv_zone_ind][2] * sin_inc
            ) / __orbital_angmom;
            if(zone.locked())
                dest -= (above_frac_deriv[deriv_zone_ind][locked_zone_ind]
                         *
                         (zone_x_torque_above - zone_x_torque_below)
                         /
                         zone.angular_momentum());
            if(std::abs(static_cast<int>(deriv_zone_ind - global_zone_ind)) <= 1) {
                double torque_deriv=
                    zone_torque_deriv[deriv_zone_ind + 1 - global_zone_ind][0];
                if(zone.locked() && deriv_zone_ind == global_zone_ind) {
                    dest-=(
                        above_frac * zone_torque_deriv[3][0]
                        +
                        (1.0 - above_frac) * torque_deriv
                    ) / zone.angular_momentum();
                } else dest -= torque_deriv / zone.angular_momentum();
            }
        }
        if(entry == Dissipation::INCLINATION)
            result[global_zone_ind] -= (
                orbit_torque[0] * sin_inc
                +
                orbit_torque[2] * cos_inc
            ) / __orbital_angmom;
        else if(entry == Dissipation::SPIN_ANGMOM)
            result[global_zone_ind] += (
                zone.locked()
                ? (above_frac * zone_x_torque_above
                   +
                   (1.0 - above_frac) * zone_x_torque_below)
                 : zone_x_torque_below
            ) / std::pow(zone.angular_momentum(), 2);
    }

    void BinarySystem::periapsis_evolution_zone_derivs(
            Dissipation::QuantityEntry entry,
            DissipatingBody &body,
            unsigned zone_ind,
            double zone_y_torque_above,
            double zone_y_torque_below,
            const std::valarray<Eigen::Vector3d> &zone_torque_deriv,
            double orbit_y_torque,
            const std::valarray<Eigen::Vector3d> &orbit_torque_deriv,
            const std::valarray<Eigen::VectorXd> &above_frac_deriv,
            double sin_inc,
            double cos_inc,
            unsigned locked_zone_ind,
            double *result
    ) const
    {
        assert(entry == Dissipation::INCLINATION
               ||
               entry == Dissipation::PERIAPSIS
               ||
               entry == Dissipation::MOMENT_OF_INERTIA
               ||
               entry == Dissipation::SPIN_ANGMOM);
        assert(orbit_torque_deriv.size()
               ==
                __body1.number_zones() + __body2.number_zones());

        DissipatingZone &zone=body.zone(zone_ind);
        unsigned num_zones = orbit_torque_deriv.size();
        unsigned global_zone_ind = zone_ind;
        if(&body == &__body2) global_zone_ind += __body1.number_zones();

        double above_frac =
            __above_lock_fractions[Dissipation::NO_DERIV][locked_zone_ind];

        for(
            unsigned deriv_zone_ind = 0;
            deriv_zone_ind < num_zones;
            ++deriv_zone_ind
        ) {
            if(entry == Dissipation::PERIAPSIS && deriv_zone_ind == 0) continue;
            double &dest = (entry == Dissipation::PERIAPSIS
                            ? result[deriv_zone_ind - 1]
                            : result[deriv_zone_ind]);
            dest = (-orbit_torque_deriv[deriv_zone_ind][1]
                    *
                    cos_inc
                    /
                    (__orbital_angmom * sin_inc));
            if(zone.locked())
                dest += (above_frac_deriv[deriv_zone_ind][locked_zone_ind]
                         *
                         (zone_y_torque_above - zone_y_torque_below)
                         /
                         (sin_inc * zone.angular_momentum()));
            if(
                std::abs(static_cast<int>(deriv_zone_ind - global_zone_ind))
                <=
                1
            ) {
                double torque_deriv =
                    zone_torque_deriv[deriv_zone_ind + 1 - global_zone_ind][1];
                if(zone.locked() && deriv_zone_ind == global_zone_ind) {
                    dest+=(
                        above_frac * zone_torque_deriv[3][1]
                        +
                        (1.0 - above_frac) * torque_deriv
                    ) / (sin_inc * zone.angular_momentum());
                } else dest += (torque_deriv
                                /
                                (sin_inc * zone.angular_momentum()));
            }
        }
        double zone_torque;
        if(zone.locked()) zone_torque = (above_frac*zone_y_torque_above
                                         +
                                         (1.0 - above_frac)
                                         *
                                         zone_y_torque_below);
        else zone_torque = zone_y_torque_below;
        if(entry == Dissipation::INCLINATION)
            result[global_zone_ind] += (
                orbit_y_torque / (__orbital_angmom * std::pow(sin_inc, 2))
                -
                zone_torque*cos_inc
                /
                (zone.angular_momentum() * std::pow(sin_inc, 2))
            );
        else if(entry == Dissipation::SPIN_ANGMOM)
            result[global_zone_ind] -=
                zone_torque / (std::pow(zone.angular_momentum(), 2) * sin_inc);
    }

    void BinarySystem::spin_angmom_evolution_zone_derivs(
            Dissipation::QuantityEntry entry,
            DissipatingBody &body,
            unsigned zone_ind,
            double zone_z_torque_above,
            double zone_z_torque_below,
            const std::valarray<Eigen::Vector3d> &zone_torque_deriv,
            const std::valarray<Eigen::VectorXd> &above_frac_deriv,
            unsigned locked_zone_ind,
            double *result
    ) const
    {
        unsigned global_zone_ind = zone_ind,
                 num_zones = number_zones();
        if(&body == &__body2) global_zone_ind += __body1.number_zones();
        double above_frac =
            __above_lock_fractions[Dissipation::NO_DERIV][locked_zone_ind];
        bool zone_is_locked = body.zone(zone_ind).locked();
        for(
            unsigned deriv_zone_ind = 0;
            deriv_zone_ind<num_zones;
            ++deriv_zone_ind
        ) {
            double &dest = (entry == Dissipation::PERIAPSIS
                            ? result[deriv_zone_ind - 1]
                            : result[deriv_zone_ind]);
            if(zone_is_locked)
                dest = (above_frac_deriv[deriv_zone_ind][locked_zone_ind]
                        *
                        (zone_z_torque_above - zone_z_torque_below));
            else dest = 0;
            if(std::abs(static_cast<int>(deriv_zone_ind
                                         -
                                         global_zone_ind)) <= 1) {
                double torque_deriv =
                    zone_torque_deriv[deriv_zone_ind + 1 - global_zone_ind][2];
                if(zone_is_locked) dest += (above_frac * zone_torque_deriv[3][2]
                                            +
                                            (1.0 - above_frac) * torque_deriv);
                else dest += torque_deriv;
            }
        }
    }

    void BinarySystem::binary_jacobian(double *param_derivs,
                                       double *age_derivs) const
    {
        unsigned body1_zones = __body1.number_zones(),
                 num_zones = body1_zones + __body2.number_zones(),
                 num_locked_zones = (__body1.number_locked_zones()
                                     +
                                     __body2.number_locked_zones()),
                 num_param = 1 + 6 * num_zones - num_locked_zones;
        double dangmom_da = __orbital_angmom / (2.0 * __semimajor),
               dangmom_de = (-__eccentricity
                             /
                             (1.0 - std::pow(__eccentricity, 2))
                             *
                             __orbital_angmom);
        std::valarray<double> orbit_power_deriv(num_param + 1),
                              orbit_angmom_gain_deriv(num_param + 1);
        fill_orbit_power_deriv(orbit_power_deriv);
        fill_orbit_angmom_gain_deriv(orbit_angmom_gain_deriv);
        semimajor_jacobian(orbit_power_deriv,
                           num_locked_zones == 0,
                           param_derivs,
                           age_derivs[0]);
        eccentricity_jacobian(orbit_power_deriv,
                              orbit_angmom_gain_deriv,
                              num_locked_zones == 0,
                              param_derivs + num_param,
                              age_derivs[1]);
        unsigned locked_zone_ind = 0;
        DissipatingBody *body = &__body1;
        unsigned body_zone_ind = 0;
        static Dissipation::QuantityEntry zone_deriv_list[]={
            Dissipation::INCLINATION,
            Dissipation::PERIAPSIS,
            Dissipation::SPIN_ANGMOM
        };
        std::valarray<double> reference_periapsis_param_deriv(num_param);
        double reference_periapsis_age_deriv;
        const std::valarray<Eigen::VectorXd> *above_frac_deriv[] = {
            &__above_lock_fractions_inclination_deriv,
            &__above_lock_fractions_periapsis_deriv,
            &__above_lock_fractions_angmom_deriv
        };
        for(unsigned zone_ind = 0; zone_ind < num_zones; ++zone_ind) {
            DissipatingZone &zone = body->zone(zone_ind);
            double sin_inc = std::sin(zone.inclination()),
                   cos_inc = std::cos(zone.inclination());
            double &periapsis_age_deriv = (
                zone_ind==0
                ? reference_periapsis_age_deriv
                : age_derivs[1 + num_zones + zone_ind]
            );
            angle_evolution_age_deriv(*body,
                                      zone_ind,
                                      sin_inc,
                                      cos_inc,
                                      locked_zone_ind,
                                      age_derivs[2 + zone_ind],
                                      periapsis_age_deriv);

            double *periapsis_row = (
                zone_ind==0
                ? &reference_periapsis_param_deriv[0]
                : param_derivs+(1 + zone_ind+num_zones) * num_param
            );

            angle_evolution_orbit_deriv(Dissipation::SEMIMAJOR,
                                        dangmom_da,
                                        *body,
                                        body_zone_ind,
                                        sin_inc,
                                        cos_inc,
                                        locked_zone_ind,
                                        param_derivs[(2 + zone_ind) * num_param],
                                        periapsis_row[0]);
            angle_evolution_orbit_deriv(
                Dissipation::ECCENTRICITY,
                dangmom_de,
                *body,
                body_zone_ind,
                sin_inc, cos_inc,
                locked_zone_ind,
                param_derivs[(2 + zone_ind) * num_param + 1],
                periapsis_row[1]
            );

            Eigen::Vector3d zone_torque_above = body->nontidal_torque(zone_ind),
                            zone_torque_below = zone_torque_above,
                            orbit_torque = body->tidal_orbit_torque(zone);
            zone_torque_above += body->tidal_torque(zone_ind, true);
            zone_torque_below += body->tidal_torque(zone_ind, false);
            std::valarray<Eigen::Vector3d>
                zone_torque_deriv(zone.locked() ? 4 : 3),
                orbit_torque_deriv(num_zones);
            unsigned param_offset = (2 + zone_ind) * num_param + 2;
            for(unsigned deriv_ind = 0; deriv_ind < 3; ++deriv_ind) {
                Dissipation::QuantityEntry deriv = zone_deriv_list[deriv_ind];
                fill_zone_torque_deriv(deriv,
                                       *body,
                                       zone_ind,
                                       zone_torque_deriv);
                fill_orbit_torque_deriv(deriv,
                                        *body,zone_ind,
                                        orbit_torque_deriv);
                inclination_evolution_zone_derivs(deriv,
                                                  *body, zone_ind,
                                                  zone_torque_above[0],
                                                  zone_torque_below[0],
                                                  zone_torque_deriv,
                                                  orbit_torque,
                                                  orbit_torque_deriv,
                                                  *(above_frac_deriv[deriv_ind]),
                                                  sin_inc,
                                                  cos_inc,
                                                  locked_zone_ind,
                                                  param_derivs + param_offset);
                periapsis_evolution_zone_derivs(deriv,
                                                *body,
                                                zone_ind,
                                                zone_torque_above[1],
                                                zone_torque_below[1],
                                                zone_torque_deriv,
                                                orbit_torque[1],
                                                orbit_torque_deriv,
                                                *(above_frac_deriv[deriv_ind]),
                                                sin_inc,
                                                cos_inc,
                                                locked_zone_ind,
                                                periapsis_row + 2);
                spin_angmom_evolution_zone_derivs(
                    deriv,
                    *body,
                    zone_ind,
                    zone_torque_above[2],
                    zone_torque_below[2],
                    zone_torque_deriv,
                    *(above_frac_deriv[deriv_ind]),
                    locked_zone_ind,
                    param_derivs + param_offset + 2 * num_zones - 1
                );
                param_offset += (deriv == Dissipation::PERIAPSIS
                                 ? num_zones - 1
                                 : num_zones);
            }

            if(zone.locked()) ++locked_zone_ind;
            ++body_zone_ind;
            if(body_zone_ind == body->number_zones()) {
                body_zone_ind = 0;
                body = &__body2;
            }
            if(zone_ind != 0) {
                periapsis_age_deriv -= reference_periapsis_age_deriv;
                for(unsigned i = 0; i < num_param; ++i)
                    periapsis_row[i] -= reference_periapsis_param_deriv[i];
            }
        }

    }
#endif

    void BinarySystem::fill_locked_surface_orbit(
        std::valarray<double> &orbit
    ) const
    {
        assert(__evolution_mode == Core::LOCKED_SURFACE_SPIN);
        assert(std::isnan(__semimajor));
        assert(std::isnan(__eccentricity));

        orbit.resize(__body1.number_zones() - 1);
        for(unsigned i = 1; i < __body1.number_zones(); ++i) {
            orbit[i - 1] = __body1.zone(i).angular_momentum();
            assert(__body1.zone(i).inclination() == 0);
            assert(__body1.zone(i).periapsis() == 0);
        }
    }

    void BinarySystem::fill_binary_orbit(std::valarray<double> &orbit) const
    {
        assert(__evolution_mode == Core::BINARY);
        assert(!std::isnan(__semimajor));
        assert(!std::isnan(__eccentricity));
        assert(__body1.zone(0).periapsis() == 0);

        orbit.resize(1 + 3 * number_zones() - number_locked_zones());
        if(number_locked_zones() == 0) orbit[0] = std::pow(__semimajor, 6.5);
        else orbit[0] = __semimajor;
        orbit[1] = __eccentricity;
        unsigned inclination_ind = 2,
                 periapsis_ind = 2 + number_zones(),
                 angmom_ind = periapsis_ind + number_zones() - 1;
        for(short body_ind = 0; body_ind < 2; ++body_ind) {
            DissipatingBody &body = (body_ind == 0 ? __body1 : __body2);
            for(
                unsigned zone_ind = 0;
                zone_ind < body.number_zones();
                ++zone_ind
            ) {
                DissipatingZone &zone = body.zone(zone_ind);
                orbit[inclination_ind++] = zone.inclination();
                if(body_ind || zone_ind)
                    orbit[periapsis_ind++] = zone.periapsis();
                if(!zone.locked())
                    orbit[angmom_ind++] = zone.angular_momentum();
            }
        }
    }

    void BinarySystem::fill_single_orbit(std::valarray<double> &orbit) const
    {
        assert(__evolution_mode == Core::SINGLE);
        assert(std::isnan(__semimajor));
        assert(std::isnan(__eccentricity));
        assert(__body1.zone(0).periapsis() == 0);
        assert(__body1.zone(0).inclination() == 0);
        assert(__body1.number_locked_zones() == 0);

        orbit.resize(3 * __body1.number_zones() - 2);
        unsigned inclination_ind = 0,
                 periapsis_ind = __body1.number_zones() - 1,
                 angmom_ind = periapsis_ind + __body1.number_zones() - 1;
        for(
            unsigned zone_ind = 0;
            zone_ind < __body1.number_zones();
            ++zone_ind
        ) {
            DissipatingZone &zone = __body1.zone(zone_ind);
            if(zone_ind) orbit[inclination_ind++] = zone.inclination();
            if(zone_ind) orbit[periapsis_ind++] = zone.periapsis();
            assert(!zone.locked());
            orbit[angmom_ind++] = zone.angular_momentum();
        }
    }

    int BinarySystem::configure(bool initialize,
                                double age,
                                double semimajor,
                                double eccentricity,
                                const double *spin_angmom,
                                const double *inclination,
                                const double *periapsis,
                                Core::EvolModeType evolution_mode)
    {
#ifndef NDEBUG
        if(initialize)
            std::cerr << "Initializing BinarySystem." << std::endl;
        if(evolution_mode != Core::BINARY) {
            assert(std::isnan(semimajor));
            assert(std::isnan(eccentricity));
        }
        if(evolution_mode == Core::LOCKED_SURFACE_SPIN) {
            assert(inclination == NULL);
            assert(periapsis == NULL);
        }

        std::cerr << "Configuring binary with a = " << semimajor
                  << ", e = " << eccentricity
                  << " in " << evolution_mode << " mode"
                  << std::endl;
#endif

        if(
            evolution_mode == Core::BINARY
            &&
            !(
                semimajor > 0.0
                &&
                eccentricity >= 0.0
                &&
                eccentricity < 1.0
            )
        )
            return GSL_EDOM;

        __evolution_mode = evolution_mode;
        double m1 = __body1.mass(),
               m2 = __body2.mass();
        __age = age;
        __semimajor = semimajor;
        __eccentricity = eccentricity;
        __orbital_energy = Core::orbital_energy(m1, m2, semimajor);
        __orbital_angmom = Core::orbital_angular_momentum(m1,
                                                          m2,
                                                          semimajor,
                                                          eccentricity);
#ifndef NDEBUG
        std::cerr << "Configuring primary." << std::endl;
#endif
        __body1.configure(initialize,
                          age,
                          m2,
                          semimajor,
                          eccentricity,
                          spin_angmom,
                          inclination,
                          periapsis,
                          evolution_mode == Core::LOCKED_SURFACE_SPIN,
                          evolution_mode != Core::BINARY,
                          true);


        if(evolution_mode == Core::BINARY) {
            unsigned offset = __body1.number_zones();
#ifndef NDEBUG
            std::cerr << "Configuring secondary." << std::endl;
#endif
            __body2.configure(initialize,
                              age,
                              m1,
                              semimajor,
                              eccentricity,
                              spin_angmom + offset - __body1.number_locked_zones(),
                              inclination + offset,
                              periapsis + offset - 1);
            find_locked_zones();
            update_above_lock_fractions();
        } else
            find_locked_zones();

        fill_orbit_torque_and_power();
#ifndef NDEBUG
    //	if(evolution_mode == BINARY)
    //		assert(__semimajor > minimum_semimajor());
#endif
        return 0;
    }

    int BinarySystem::configure(bool initialize,
                                double age,
                                const double *parameters,
                                Core::EvolModeType evolution_mode)
    {
        __evolution_mode = evolution_mode;
        double semimajor, eccentricity;
        const double *spin_angmom, *inclination, *periapsis;
        unsigned num_zones = number_zones();
        if(evolution_mode == Core::BINARY) {
            if(parameters[0] < 0) {
#ifndef NDEBUG
                std::cerr << "At t = " << age << " param: ";
                for(
                    unsigned i = 0;
                    i < 3 * num_zones + 1 - number_locked_zones();
                    ++i
                ) {
                    if(i) std::cerr << ", ";
                    std::cerr << parameters[i];
                }
                std::cerr << std::endl;
#endif
                return GSL_EDOM;
            }
            if(__body1.number_locked_zones() || __body2.number_locked_zones()) {
                semimajor = parameters[0];
            } else
                semimajor = std::pow(parameters[0], 1.0 / 6.5);
            eccentricity = parameters[1];
            inclination = parameters+2;
        } else {
            semimajor = eccentricity=Core::NaN;
            if(evolution_mode == Core::SINGLE) inclination = parameters;
            else inclination = NULL;
        }

        if(evolution_mode == Core::LOCKED_SURFACE_SPIN) {
            periapsis = NULL;
            spin_angmom = parameters;
        } else {
            assert(inclination != NULL);

            periapsis = inclination+num_zones;
            if(evolution_mode == Core::SINGLE) --periapsis;
            spin_angmom = periapsis + num_zones - 1;
        }
        return configure(initialize,
                         age,
                         semimajor,
                         eccentricity,
                         spin_angmom,
                         inclination,
                         periapsis,
                         evolution_mode);
    }

    Core::EvolModeType BinarySystem::fill_orbit(
        std::valarray<double> &orbit
    ) const
    {
        if(__evolution_mode == Core::LOCKED_SURFACE_SPIN)
            fill_locked_surface_orbit(orbit);
        else if(__evolution_mode == Core::BINARY) fill_binary_orbit(orbit);
        else fill_single_orbit(orbit);
        return __evolution_mode;
    }

    double BinarySystem::above_lock_fraction(unsigned locked_zone_index,
                                             Dissipation::QuantityEntry entry,
                                             unsigned deriv_zone_index,
                                             bool secondary_radius)
    {
        if(zone_specific(entry)) {
            assert(entry == Dissipation::INCLINATION
                   ||
                   entry == Dissipation::PERIAPSIS
                   ||
                   entry == Dissipation::SPIN_ANGMOM
                   ||
                   entry == Dissipation::MOMENT_OF_INERTIA);

            if(entry == Dissipation::INCLINATION)
                return __above_lock_fractions_inclination_deriv[deriv_zone_index]
                    [locked_zone_index];
            else if(entry == Dissipation::PERIAPSIS)
                return __above_lock_fractions_periapsis_deriv[deriv_zone_index]
                    [locked_zone_index];
            else if(entry == Dissipation::SPIN_ANGMOM)
                return __above_lock_fractions_angmom_deriv[deriv_zone_index]
                    [locked_zone_index];
            else return __above_lock_fractions_inertia_deriv[deriv_zone_index]
                    [locked_zone_index];
        } else {
            if(entry == Dissipation::RADIUS && secondary_radius)
                return
                    __above_lock_fractions_body2_radius_deriv[locked_zone_index];
            else return __above_lock_fractions[entry][locked_zone_index];
        }
    }

    int BinarySystem::differential_equations(double age,
                                             const double *parameters,
                                             Core::EvolModeType evolution_mode,
                                             double *differential_equations)
    {
#ifndef NDEBUG
        std::cerr << "Finding differential equations at t = " << age
                  << " in " << evolution_mode
                  << " mode, with orbit[0] = " << parameters[0]
                  << std::endl;
#endif
        int status = configure(false, age, parameters, evolution_mode);
        if(status != GSL_SUCCESS) return status;
        __semimajor_rate = __eccentricity_rate = Core::NaN;
        switch(evolution_mode) {
            case Core::LOCKED_SURFACE_SPIN :
                return locked_surface_differential_equations(
                    differential_equations
                );
            case Core::SINGLE :
                return single_body_differential_equations(
                    differential_equations
                );
            case Core::BINARY :
                return binary_differential_equations(differential_equations);
            default :
                throw Core::Error::BadFunctionArguments(
                    "Evolution mode other than LOCKED_SURFACE_SPIN, SINGLE or "
                    "BINARY encountered in "
                    "BinarySystem::differential_equations!"
                );
        }
    }

#ifdef ENABLE_DERIVATIVES
    int BinarySystem::jacobian(double age,
                               const double *parameters,
                               Core::EvolModeType evolution_mode,
                               double *param_derivs,
                               double *age_derivs)
    {
        configure(false, age, parameters, evolution_mode);
        switch(evolution_mode) {
            case Core::LOCKED_SURFACE_SPIN :
                locked_surface_jacobian(param_derivs, age_derivs);
                return 0;
            case Core::SINGLE : single_body_jacobian(param_derivs, age_derivs);
                                return 0;
            case Core::BINARY : binary_jacobian(param_derivs, age_derivs);
                                return 0;
            default : throw Core::Error::BadFunctionArguments(
                          "Evolution mode other than LOCKED_SURFACE_SPIN, "
                          "SINGLE or BINARY encountered in "
                          "BinarySystem::jacobian!"
                      );
        }
    }
#endif

    void BinarySystem::initialize_locks(double sync_precision)
    {
        lock_scenario_type lock_candidates = find_synchronized_zones(
            sync_precision
        );
        if(lock_candidates.empty())
            return;

        lock_scenario_type lock_scenario;
        if(!explore_lock_scenarios(lock_candidates.begin(),
                                   lock_candidates.size(),
                                   lock_scenario))
            throw Core::Error::Runtime("Failed to find viable lock scenario!");
#ifndef NDEBUG
        std::cerr << "Setting lock scenario: " << std::endl;
        describe_lock_scenario(
            std::cerr,
            __selected_lock_scenario,
            std::vector<bool>(__selected_lock_scenario.size(), true),
            true
        );
        set_lock_scenario(__selected_lock_scenario);
#endif
    }

    void BinarySystem::check_for_lock(int orbital_freq_mult,
                                      int spin_freq_mult,
                                      unsigned short body_index,
                                      unsigned zone_index,
                                      short direction)
    {
        DissipatingBody &body = (body_index ? __body2 : __body1);

        assert(body_index <= 1);
        assert(zone_index < body.number_zones());
        assert(spin_freq_mult);

        double original_angmom = body.zone(zone_index).angular_momentum();
        body.lock_zone_spin(zone_index, orbital_freq_mult, spin_freq_mult);
        unsigned num_zones = number_zones(),
                 num_locked_zones = number_locked_zones(),
                 locked_zone_ind = 0;
        std::vector<double> spin_angmom(num_zones - num_locked_zones),
                            inclinations(num_zones),
                            periapses(num_zones);
        for(unsigned zone_ind = 0; zone_ind < num_zones; ++zone_ind) {
            DissipatingZone &zone=(zone_ind < __body1.number_zones()
                                   ? __body1.zone(zone_ind)
                                   : __body2.zone(zone_ind
                                                  -
                                                  __body1.number_zones()));
            if(zone.locked()) ++locked_zone_ind;
            else spin_angmom[zone_ind-locked_zone_ind] = (
                zone.angular_momentum()
            );
            inclinations[zone_ind] = zone.inclination();
            if(zone_ind) periapses[zone_ind - 1] = zone.periapsis();
        }
        assert(locked_zone_ind == num_locked_zones);
        configure(false,
                  __age,
                  __semimajor,
                  __eccentricity,
                  &(spin_angmom[0]),
                  &(inclinations[0]),
                  &(periapses[0]),
                  Core::BINARY);
        DissipatingZone &locked_zone = body.zone(zone_index);
        double above_lock_fraction = __above_lock_fractions
                                     [Dissipation::NO_DERIV]
                                     [locked_zone.locked_zone_index()];
#ifndef NDEBUG
        std::cerr << "Holding lock requires above lock fraction of: "
                  << above_lock_fraction
                  << std::endl;
#endif
        if(above_lock_fraction > 0 && above_lock_fraction < 1) return;
        std::vector<double>::iterator check_zone_dest = (
            spin_angmom.begin()
            +
            (
                zone_index
                +
                (body_index ? __body1.number_zones() : 0)
                -
                locked_zone.locked_zone_index()
            )
        );
        spin_angmom.insert(check_zone_dest, original_angmom);
        if(direction == 0)
            throw Core::Error::Runtime(
                "Crossing spin-orbit synchronization with unspecified direction!"
            );
        body.unlock_zone_spin(zone_index, direction);
        configure(false,
                  __age,
                  __semimajor,
                  __eccentricity,
                  &(spin_angmom[0]),
                  &(inclinations[0]),
                  &(periapses[0]),
                  Core::BINARY);
    }

    double BinarySystem::minimum_separation(bool deriv) const
    {
        double rroche = (
            __body2.radius()
            ?(
                2.44
                *
                __body2.radius()
                *
                std::pow(__body1.mass() / __body2.mass(), 1.0 / 3.0)
            )
            :(
                0.0
            )
        );
        if(rroche > __body1.radius()) {
            if(deriv) return rroche * __body2.radius(1) / __body2.radius();
            else return rroche;
        } else {
            return __body1.radius(deriv ? 1 : 0);
        }
    }

    void BinarySystem::secondary_died()
    {
#ifndef NDEBUG
        std::cerr << "Handling secondary death!" << std::endl;
#endif
        unsigned num_zones = __body1.number_zones();
        std::valarray<double> spin_angmom(num_zones),
                              inclination(num_zones - 1),
                              periapsis(num_zones - 1);
        DissipatingZone &old_surface_zone = __body1.zone(0);
        double old_surface_inclination = old_surface_zone.inclination(),
               sin_inc = std::sin(old_surface_inclination),
               cos_inc = std::cos(old_surface_inclination),
               old_surface_angmom = old_surface_zone.angular_momentum(),
               angmom = std::sqrt(
                   std::pow(old_surface_angmom + __orbital_angmom * cos_inc, 2)
                   +
                   std::pow(__orbital_angmom * sin_inc, 2)
               ),
               new_surface_inclination=std::atan2(__orbital_angmom*sin_inc,
                                                  old_surface_angmom
                                                  +
                                                  __orbital_angmom*cos_inc);

        spin_angmom[0] = angmom;
    //	assert(num_zones == 2);
        if(old_surface_zone.locked()) __body1.unlock_zone_spin(0, 1);
        for(unsigned zone_ind = 1; zone_ind < num_zones; ++zone_ind) {
            DissipatingZone &zone = __body1.zone(zone_ind);
            assert(zone.periapsis() == 0);
            inclination[zone_ind - 1]=std::abs(zone.inclination()
                                               -
                                               old_surface_inclination
                                               +
                                               new_surface_inclination);
            periapsis[zone_ind - 1] = 0;
            spin_angmom[zone_ind] = zone.angular_momentum();
            if(zone.locked()) __body1.unlock_zone_spin(zone_ind, 1);
        }
        configure(false,
                  __age,
                  Core::NaN,
                  Core::NaN,
                  &spin_angmom[0],
                  &inclination[0],
                  &periapsis[0],
                  Core::SINGLE);
        __body1.spin_jumped();
    }

    void BinarySystem::release_lock(unsigned locked_zone_index, short direction)
    {
        DissipatingBody *body;
        if(locked_zone_index < __body1.number_locked_zones()) {
            body = &__body1;
            for(
                unsigned zone_ind = 0;
                zone_ind < __body2.number_locked_zones();
                ++zone_ind
            )
                if(__body2.zone(zone_ind).locked())
                    --__body2.zone(zone_ind).locked_zone_index();

        } else {
            body = &__body2;
            locked_zone_index -= __body1.number_locked_zones();
        }
        unsigned zone_ind = 0;
        while(true){
            if(body->zone(zone_ind).locked()) {
                if(locked_zone_index == 0) break;
                else --locked_zone_index;
            }
            assert(zone_ind < body->number_zones());
            ++zone_ind;
        }
        body->unlock_zone_spin(zone_ind, direction);
        for(; zone_ind < body->number_zones(); ++zone_ind)
            if(body->zone(zone_ind).locked())
                --body->zone(zone_ind).locked_zone_index();
    }

    void BinarySystem::add_to_evolution()
    {
        __semimajor_evolution.push_back(__semimajor);
        __eccentricity_evolution.push_back(__eccentricity);
        __semimajor_rate_evolution.push_back(__semimajor_rate);
        __eccentricity_rate_evolution.push_back(__eccentricity_rate);
        __body1.add_to_evolution();
        __body2.add_to_evolution();
    }

    void BinarySystem::reset_evolution()
    {
        __semimajor_evolution.clear();
        __eccentricity_evolution.clear();
        __semimajor_rate_evolution.clear();
        __eccentricity_rate_evolution.clear();
        __body1.reset_evolution();
        __body2.reset_evolution();
    }

    void BinarySystem::rewind_evolution(unsigned nsteps)
    {
        for(unsigned i = 0; i < nsteps; ++i) {
            __semimajor_evolution.pop_back();
            __eccentricity_evolution.pop_back();
            __semimajor_rate_evolution.pop_back();
            __eccentricity_rate_evolution.pop_back();
        }
        __body1.rewind_evolution(nsteps);
        __body2.rewind_evolution(nsteps);
    }

    CombinedStoppingCondition *BinarySystem::stopping_conditions()
    {
        CombinedStoppingCondition *result = new CombinedStoppingCondition();
        if(__evolution_mode == Core::BINARY)
            (*result) |= new SecondaryDeathCondition(*this);
        (*result) |= __body1.stopping_conditions(*this, true);
        if(__evolution_mode == Core::BINARY)
            (*result) |= __body2.stopping_conditions(*this, false);
        return result;
    }

    void BinarySystem::reached_critical_age(double age)
    {
        __body1.reached_critical_age(age);
        __body2.reached_critical_age(age);
    }

    double BinarySystem::next_stop_age() const
    {
        return std::min(__body1.next_stop_age(),
                        __body2.next_stop_age());
    }

    unsigned BinarySystem::expansion_order() const
    {
#ifndef NDEBUG
        int result = -1;
#endif
        DissipatingBody *body = &__body1;
        while(true) {
            for(
                unsigned zone_ind = 0;
                zone_ind < body->number_zones();
                ++zone_ind
            ) {
                DissipatingZone &zone = body->zone(zone_ind);
                if(zone.dissipative()) {
#ifndef NDEBUG
                    if(result < 1)
                        result = zone.expansion_order();
                    else
                        assert(static_cast<unsigned>(result)
                               ==
                               zone.expansion_order());
#else
                    return zone.expansion_order();
#endif
                }
            }
            if(body == &__body2) break;
            body = &__body2;
        }
#ifndef NDEBUG
        if(result >= 1)
            return result;
#endif
        throw Core::Error::BadFunctionArguments(
            "System contains no dissipative zones! "
            "Run as single objects instead!"
        );
    }

} //End Evolve namespace.
