/**\file
 *
 * \brief Unit tests that check the machinery for monitoring potential tidal
 * locks.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_LOCK_MONITORING_H
#define __TEST_LOCK_MONITORING_H

#include <sstream>
#include <boost/math/special_functions/sign.hpp>

#include "../shared/Common.h"
#include "../../Planet/PlanetZone.h"
#include "../../Planet/Planet.h"

namespace Evolve {
    /**\brief The test suite that ensures the correct locks are selected for
     * monitoring and the fixes applied to the tidal frequency are correct.
     *
     * \ingroup UnitTests_group
     */
    class test_LockMonitoring : public Test::Suite {
    private:
        ///A dissipative zone to run the tests on.
        Planet::PlanetZone __zone;

        ///Change the order of the tidal potential expansion used by __zone.
        void set_expansion_order(int max_orbital_multiplier);

        ///Change the orbital and spin frequencies for __zone.
        void set_orbital_spin_frequeqncy(double orbital_frequency,
                                         double spin_frequency);

        ///\brief Add information about the given locks relative to the tidal
        ///frequency to the test message.
        void add_locks_to_message(
            const SpinOrbitLockInfo &lock_below,
            const SpinOrbitLockInfo &lock_above,
            std::ostringstream &message
        );

        ///Add the locks currently selected for monitoring in __zone to the
        ///message.
        void add_zone_locks_to_message(std::ostringstream &message);

        ///\brief Verify a single expectation of what locks would be monitored
        ///in a given situation.
        void check_monitored_locks(
            double orbital_frequency,
            double spin_frequency,
            const SpinOrbitLockInfo &expected_lock_below,
            const SpinOrbitLockInfo &expected_lock_above
        );

        ///Pick a just inside each boundary or in the middle of a range.
        double sample_range(
            ///The range to sample.
            std::pair<double, double> range,

            ///If 0 pick the smallest float above the range min
            ///If 1 pick the middle of the range
            ///If 2 pick the largest float below range max
            unsigned position
        );

        ///\brief Verify a single expectation for fixing of a tidal term in a
        ///given situation over a range of spin frequencies.
        void check_tidal_frequency_range(
            ///The orbital frequency to use when setting the locks of a zone.
            double lock_set_orbital_frequency,

            ///The spin frequency to use when setting locks of a zone.
            double lock_set_spin_frequency,

            ///The orbital frequency at which to calculate the possibly
            ///corrected tidal frequency.
            double test_orbital_frequency,

            ///The range of spin frequencies at which to calculate the possibly
            ///corrected tidal frequency. Comparison is done near each end of
            ///the range and the middle.
            std::pair<double, double> test_spin_frequency_range,

            ///The multiplier in front of the orbital frequency of the tidal
            //term being tested.
            int orbital_frequency_multiplier,

            ///The multiplier in front of the spin frequency of the tidal
            //term being calculated.
            int spin_frequency_multiplier,

            /// 1 if the tidal frequency should be corrected to positive
            /// 0 if the tidal frequency should not be corrected
            /// -1 if the tidal frequency should be corrected to negative
            int expected_correction
        );

        ///\brief Check that lock-corrected tidal frequency reproduces
        ///expectations for a single set of monitored locks.
        void check_unlocked_tidal_frequency(
            ///The largest orbital frequency multiplier to include in the
            ///expansion.
            int expansion_order,

            ///The orbital frequency multiplier of the tidal term defining the
            ///lower limit on the spin frequency to watch, assimung the spin
            ///frequency multiplier is 2.
            int lower_lock_orbital_freuqency_multiplier,

            ///All tests are performed at this fixed orbital frequency.
            double orbital_frequency
        );
    protected:
        ///No fixture at this point
        void setup() {};

        ///No fixture at this point
        void tear_down() {};

    public:
        ///Create the test suite.
        test_LockMonitoring();

        ///\brief Test whether the correct tidal terms are selected for lock
        ///monitoring when \f$\Omega_spin>0\f$ and
        void test_lock_init();

        ///\brief Test fix applied to the tidal frequency of tidal term(s)
        ///matching lower boundary lock.
        void test_tidal_frequency_fix();
    }; //End test_LockMonitoring class.

} //End Evolve namespace.

#endif
