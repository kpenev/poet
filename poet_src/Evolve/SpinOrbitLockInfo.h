/**\file
 * 
 * \brief Defines the SpinOrbitLockInfo class.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __SPIN_ORBIT_LOCK_INFO
#define __SPIN_ORBIT_LOCK_INFO

#include "../Core/SharedLibraryExportMacros.h"
#include <iostream>

namespace Evolve {

    ///\brief Defines a lock between the spin of a dissipating body and the
    ///orbit.
    ///
    ///With inclined and eccentric orbits, locks can occur at many different
    //frequencies, not only when the orbital and spin periods are the same. In
    ///general almost any rational ratio can result in a lock if the dissipation
    ///has the appropriate frequency dependence.
    ///
    ///\ingroup StellarSystem_group
    class LIB_LOCAL SpinOrbitLockInfo {
    private:
        ///The mutiplier in front of the orbital frequency in the lock.
        int __orbital_freq_mult,

            ///The multiplier in front of the spin frequency in the lock.
            __spin_freq_mult;

        ///\brief Should a lock be assumed, and if so from which direction is it
        ///approached?
        ///
        ///The values have the following meanings:
        /// - <0 	: the spin frequency of the body is slightly lagrer than
        ///			  necessary for a lock (or the orbital freqency is slightly
        ///		      smaller).
        ///
        /// - 0 	: The spin frequency and the orbital frequency have precisely
        ///			  the values to result in zero forcing frequency for this
        ///			  term.
        ///
        /// - >0 	: the spin frequency of the body is slightly smaller than
        ///			  necessary for a lock (or the orbital freqency is slightly
        ///		      larger).
        short __lock_direction;

    public:
        ///\brief Define which tidal dissipation term is in a lock.
        SpinOrbitLockInfo(
                ///The multiple of the orbital frequency at the lock.
                int orbital_freq_mult=0,

                ///The multiple of the spin frequency at the lock.
                int spin_freq_mult=0,

                ///The direction from which the spin frequency is approaching the
                ///lock. See #__lock_direction for the meaning of the values.
                short lock_direction=-1)
        {set_lock(orbital_freq_mult, spin_freq_mult, lock_direction);}

        ///\brief Copy the original to this.
        SpinOrbitLockInfo(const SpinOrbitLockInfo &orig)
        {set_lock(orig.orbital_frequency_multiplier(),
                orig.spin_frequency_multiplier(),
                orig.lock_direction());}

        ///\brief Make this the same lock as RHS.
        SpinOrbitLockInfo &operator=(const SpinOrbitLockInfo &rhs)
        {set_lock(rhs.orbital_frequency_multiplier(),
                rhs.spin_frequency_multiplier(),
                rhs.lock_direction()); return *this;}

        ///\brief Define which tidal dissipation term is in a lock.
        void set_lock(
                ///The multiple of the orbital frequency at the lock.
                int orbital_freq_mult,

                ///The multiple of the spin frequency at the lock.
                int spin_freq_mult,

                ///The sign of the forcing frequency for this term.
                ///See #__lock_direction for the meaning of the values.
                short lock_direction=0);

        ///\brief Spin frequency at exactly the lock that corresponds to the
        ///given orbital frequency.
        double spin(double orbital_frequency) const
        {return (orbital_frequency*__orbital_freq_mult)/__spin_freq_mult;}

        ///Is the given tidal dissipation term one of the locked terms?
        bool operator()(
                ///The multiple of the orbital frequency to consider.
                int orbital_freq_mult,

                ///The multiple of the spin frequency to consider.
                int spin_freq_mult) const
        {return (__lock_direction ? false :
                 term(orbital_freq_mult, spin_freq_mult));}

        ///Returns true if the lock is referring to the given term, regardless of
        ///whether it is locked or not.
        bool term(			
                ///The multiple of the orbital frequency to consider.
                int orbital_freq_mult,

                ///The multiple of the spin frequency to consider.
                int spin_freq_mult) const
        {return (orbital_freq_mult * __spin_freq_mult
                 ==
                 spin_freq_mult * __orbital_freq_mult);}

        ///Should this lock be assumed.
        operator bool() const {return __lock_direction==0;}

        ///\rbief The opposite of the sign of the forcing frequency associated
        ///with this component.
        short lock_direction() const {return __lock_direction;}

        ///Set the lock direction to the given value.
        void lock_direction(short value) {__lock_direction=value;}

        ///The multiplier in front of the orbital frequency in the lock.
        int orbital_frequency_multiplier() const {return __orbital_freq_mult;}

        ///The multiplier in front of the spin frequency in the lock.
        int spin_frequency_multiplier() const {return __spin_freq_mult;}

        ///Are the two locks for the same frequency ratio and in the same
        ///enabled/disabled state.
        bool operator==(const SpinOrbitLockInfo &rhs) const;
    };//End SpinOrbitLockInfo class.

    ///Civilized output for locks.
    LIB_LOCAL std::ostream &operator<<(std::ostream &os,
                                       const SpinOrbitLockInfo &lock);

} //End Evolve namespace.

#endif
