#include "DissipatingBody.h"

namespace Evolve {

    std::ostream &operator<<(std::ostream &os, const SpinOrbitLockInfo &lock)
    {
        os << "Lock: " << lock.orbital_frequency_multiplier() << "*OrbFreq = "
            << lock.spin_frequency_multiplier() << "*SpinFreq ";
        if(!lock) os << "DISABLED";
        else os << (lock.lock_direction()>0 ? "from above" : "from below");
        return os;
    }

    void SpinOrbitLockInfo::set_lock(int orbital_freq_mult, int spin_freq_mult,
            short lock_direction)
    {
        __orbital_freq_mult=orbital_freq_mult;
        __spin_freq_mult=spin_freq_mult;
        __lock_direction=lock_direction;
    }

    bool SpinOrbitLockInfo::operator==(const SpinOrbitLockInfo &rhs) const
    {
        return rhs.__orbital_freq_mult==__orbital_freq_mult &&
            rhs.__spin_freq_mult==__spin_freq_mult &&
            bool(rhs.__lock_direction)==bool(__lock_direction);
    }

} //End Evolve namespace.
