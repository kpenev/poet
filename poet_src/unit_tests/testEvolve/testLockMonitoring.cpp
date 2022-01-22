/**\file
 *
 * \brief Implement the non-inline methods of test_LockMonitoring.
 *
 * \ingroup UnitTests_group
 */

#include "testLockMonitoring.h"

namespace Evolve {

    const double epsilon = std::numeric_limits<double>::epsilon();

    void test_LockMonitoring::set_expansion_order(int max_orbital_multiplier)
    {
        Planet::Planet body1(1.0, 1.0), body2(1.0, 1.0);
        BinarySystem system(body1, body2);
        __zone.change_expansion_order(max_orbital_multiplier, system, true, 0);
    }

    void test_LockMonitoring::set_orbital_spin_frequeqncy(
        double orbital_frequency,
        double spin_frequency
    )
    {
        __zone.configure(true,
                         0.0,
                         orbital_frequency,
                         0.0,
                         1.0,
                         spin_frequency,
                         0.0,
                         0.0,
                         true);
    }

    void test_LockMonitoring::add_locks_to_message(
        const SpinOrbitLockInfo &lock_below,
        const SpinOrbitLockInfo &lock_above,
        std::ostringstream &message
    )
    {
        message << lock_below.orbital_frequency_multiplier()
                << "/"
                << lock_below.spin_frequency_multiplier()
                << " * Worb "
                << " < W* < "
                << lock_above.orbital_frequency_multiplier()
                << "/"
                << lock_above.spin_frequency_multiplier()
                << " * Worb";
    }

    void test_LockMonitoring::add_zone_locks_to_message(
        std::ostringstream &message
    )
    {
        if(__zone.lock_monitored(false).lock_direction() > 0) {
            add_locks_to_message(__zone.lock_monitored(false),
                                 __zone.lock_monitored(true),
                                 message);

        } else {
            add_locks_to_message(__zone.lock_monitored(true),
                                 __zone.lock_monitored(false),
                                 message);
        }

    }

    void test_LockMonitoring::check_monitored_locks(
        double orbital_frequency,
        double spin_frequency,
        const SpinOrbitLockInfo &expected_lock_below,
        const SpinOrbitLockInfo &expected_lock_above
    )
    {
        set_orbital_spin_frequeqncy(orbital_frequency, spin_frequency);
        std::ostringstream message;
        message << "For Worb = " << orbital_frequency
                << ", W* = " << spin_frequency
                << ", should monitor ";
        add_locks_to_message(expected_lock_below,
                             expected_lock_above,
                             message);

        message << ". auto selected terms to monitor: ";
        add_zone_locks_to_message(message);
        message << ".";

        SpinOrbitLockInfo selected_lock_below, selected_lock_above;
        if(__zone.lock_monitored(false).lock_direction() > 0) {
            selected_lock_below = __zone.lock_monitored(false);
            selected_lock_above = __zone.lock_monitored(true);
        } else {
            selected_lock_below = __zone.lock_monitored(true);
            selected_lock_above = __zone.lock_monitored(false);
        }


        TEST_ASSERT_MSG(
            expected_lock_below.term(
                selected_lock_below.orbital_frequency_multiplier(),
                selected_lock_below.spin_frequency_multiplier()
            )
            &&
            (
                expected_lock_below.lock_direction()
                *
                selected_lock_below.lock_direction()
                >
                0
            ),
            message.str().c_str()
        );

        TEST_ASSERT_MSG(
            expected_lock_above.term(
                selected_lock_above.orbital_frequency_multiplier(),
                selected_lock_above.spin_frequency_multiplier()
            )
            &&
            (
                expected_lock_above.lock_direction()
                *
                selected_lock_above.lock_direction()
                >
                0
            ),
            message.str().c_str()
        );

    }

    double test_LockMonitoring::sample_range(std::pair<double, double> range,
                                             unsigned position)
    {
        std::cerr.precision(19);
        if(position == 0) {
            if(range.first < 0) {
                return range.first * (1.0 - epsilon);
            } else if(range.first == 0)
                return epsilon;
            else
                return range.first * (1.0 + epsilon);
        } else if(position == 1)
            return 0.5 * (range.first + range.second);
        else {
            if(range.second < 0)
                return range.second * (1.0 + epsilon);
            else if(range.second == 0)
                return -epsilon;
            else
                return range.second * (1.0 - epsilon);
        }
    }


    void test_LockMonitoring::check_tidal_frequency_range(
        double lock_set_orbital_frequency,
        double lock_set_spin_frequency,
        double test_orbital_frequency,
        std::pair<double, double> test_spin_frequency_range,
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        int expected_correction
    )
    {
        set_orbital_spin_frequeqncy(lock_set_orbital_frequency,
                                    lock_set_spin_frequency);

        for(unsigned spin_choice = 0; spin_choice <= 2; ++spin_choice) {
            double test_spin_frequency = sample_range(test_spin_frequency_range,
                                                      spin_choice);
            __zone.configure(false,
                             0.0,
                             std::numeric_limits<double>::quiet_NaN(),
                             0.0,
                             1.0,
                             test_spin_frequency,
                             0.0,
                             0.0,
                             true);
            double expected_tidal_frequency = (
                expected_correction == 0
                ? (
                    orbital_frequency_multiplier * test_orbital_frequency
                    -
                    spin_frequency_multiplier * test_spin_frequency
                )
                : expected_correction * epsilon
            );
            double calculated_tidal_frequency = __zone.forcing_frequency(
                orbital_frequency_multiplier,
                spin_frequency_multiplier,
                test_orbital_frequency
            );

            std::ostringstream message;
            message.precision(19);
            message << "With zone locks: ";
            add_zone_locks_to_message(message);
            message << ", Wtide(m=" << spin_frequency_multiplier
                    << ", s=" << orbital_frequency_multiplier
                    << ") for Worb = " << test_orbital_frequency
                    << ", W* = " << test_spin_frequency
                    << " is " << calculated_tidal_frequency
                    << " instead of " << expected_tidal_frequency;

            if(expected_correction == 0) {
                double difference = (expected_tidal_frequency
                                     -
                                     calculated_tidal_frequency);
                message << ", difference (expected - got): " << difference << ".";
                TEST_ASSERT_MSG(
                    (
                        std::abs(difference) < 1e-14
                        &&
                        (
                            expected_tidal_frequency
                            *
                            calculated_tidal_frequency
                            >
                            0
                            ||
                            (
                                expected_tidal_frequency == 0
                                &&
                                calculated_tidal_frequency == 0
                            )
                        )
                    ),
                    message.str().c_str()
                )
            } else {
                message << ".";
                TEST_ASSERT_MSG(
                    calculated_tidal_frequency == expected_tidal_frequency,
                    message.str().c_str()
                )
            }
        }
    }

    void test_LockMonitoring::check_unlocked_tidal_frequency(
        int expansion_order,
        int lower_lock_orbital_freuqency_multiplier,
        double orbital_frequency
    )
    {
        set_expansion_order(expansion_order);
        double lock_set_spin_frequency = (
            0.25
            +
            0.5 * lower_lock_orbital_freuqency_multiplier
        ) * orbital_frequency;
        for(
            int spin_range_s = -3*expansion_order;
            spin_range_s <= 3*expansion_order;
            ++spin_range_s
        ) {
            for(int test_m = -2; test_m <= 2; test_m += 1) {
                for(
                    int test_s = -expansion_order;
                    test_s <= expansion_order;
                    ++test_s
                ) {
                    std::pair<double, double> spin_frequency_range(
                        (0.5 * spin_range_s) * orbital_frequency,
                        0.5 * (spin_range_s + 1) * orbital_frequency
                    );
                    int expected_correction = 0;
                    if(test_m != 0) {
                        int scaled_test_s = 2 * test_s / test_m;
                        if(
                            (
                                scaled_test_s
                                >
                                lower_lock_orbital_freuqency_multiplier
                            )
                            &&
                            (
                                scaled_test_s <= spin_range_s
                            )
                        )
                            expected_correction = boost::math::sign(test_m);

                        if(
                            (
                                scaled_test_s
                                <=
                                lower_lock_orbital_freuqency_multiplier
                            )
                            &&
                            (
                                scaled_test_s
                                >
                                spin_range_s
                            )
                        )
                            expected_correction = -boost::math::sign(test_m);
                    }

                    check_tidal_frequency_range(
                        orbital_frequency,
                        lock_set_spin_frequency,
                        orbital_frequency,
                        spin_frequency_range,
                        test_s,
                        test_m,
                        expected_correction
                    );
                }
            }
        }
    }


    test_LockMonitoring::test_LockMonitoring()
        : __zone(1.0, 1.0, 1.0)
    {
        __zone.setup(
            std::vector<double>(),
            std::vector<double>(),
            std::vector<double>(1, 0.0),
            std::vector<double>(1, 0.0),
            1.0,
            1.0,
            1.0
        );

        TEST_ADD(test_LockMonitoring::test_lock_init);
        TEST_ADD(test_LockMonitoring::test_tidal_frequency_fix);
    }

    void test_LockMonitoring::test_lock_init()
    {
        set_expansion_order(3);
        check_monitored_locks(0.32,
                              0.10,
                              SpinOrbitLockInfo(0, 1, 1),
                              SpinOrbitLockInfo(1, 2, -1));

        check_monitored_locks(0.41,
                              0.30,
                              SpinOrbitLockInfo(1, 2, 1),
                              SpinOrbitLockInfo(1, 1, -1));

        check_monitored_locks(0.46,
                              0.57,
                              SpinOrbitLockInfo(1, 1, 1),
                              SpinOrbitLockInfo(3, 2, -1));


        check_monitored_locks(0.27,
                              0.67,
                              SpinOrbitLockInfo(2, 1, 1),
                              SpinOrbitLockInfo(3, 1, -1));

        check_monitored_locks(0.10,
                              0.69,
                              SpinOrbitLockInfo(3, 1, 1),
                              SpinOrbitLockInfo(1, 0, -1));

        check_monitored_locks(0.59,
                              -0.12,
                              SpinOrbitLockInfo(-1, 2, 1),
                              SpinOrbitLockInfo(0, 1, -1));

        check_monitored_locks(0.45,
                              -0.32,
                              SpinOrbitLockInfo(-1, 1, 1),
                              SpinOrbitLockInfo(-1, 2, -1));

        check_monitored_locks(0.37,
                              -0.45,
                              SpinOrbitLockInfo(-3, 2, 1),
                              SpinOrbitLockInfo(-1, 1, -1));

        check_monitored_locks(0.25,
                              -0.61,
                              SpinOrbitLockInfo(-3, 1, 1),
                              SpinOrbitLockInfo(-2, 1, -1));

        check_monitored_locks(0.12,
                              -0.67,
                              SpinOrbitLockInfo(-1, 0, 1),
                              SpinOrbitLockInfo(-3, 1, -1));
    }

    void test_LockMonitoring::test_tidal_frequency_fix()
    {
        const int expansion_order = 10;
        for(
            int orbital_frequency_multiplier = -3 * expansion_order;
            orbital_frequency_multiplier <= 3 * expansion_order;
            ++orbital_frequency_multiplier
        )
            check_unlocked_tidal_frequency(expansion_order,
                                           orbital_frequency_multiplier,
                                           0.32);
    }

} //End Evolve namespace.
