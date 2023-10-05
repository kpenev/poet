/**\file
 *
 * \brief Declares a class for zones with non-zere dissipation for only a single
 * tidal term (and it's negative coefficient opposite).
 *
 * \ingroup Evolve_group
 */

#ifndef __SINGLE_TERM_ZONE_H
#define __SINGLE_TERM_ZONE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "DissipatingZone.h"

namespace Evolve {
    ///brief A DissipatingZone where only a single term of the tidal potential
    ///contributes to the dissipation.
    class LIB_PUBLIC SingleTermZone : virtual public DissipatingZone {
    private:
        int
            ///The orbital frequency multiplier of the only dissiaptive term.
            __orbital_frequency_multiplier,

            ///The spin frequency multiplier of the only dissipative term
            __spin_frequency_multiplier;

        ///The phase lag of the only dissipative term.
        double __phase_lag;
    public:
        ///\brief create a zone specifying only a single dissipative term.
        SingleTermZone(
            ///The multiplier of the orbital frequency in the expression for
            ///the forcing frequency for the only dissipative term.
            int orbital_frequency_multiplier=0,

            ///The multiplier of the spin frequency in the expression for
            ///the forcing frequency for the only dissipative term.
            int spin_frequency_multiplier=0,

            ///The phase lag to assume for the only dissipative term.
            double phase_lag=0
        ) :
            __orbital_frequency_multiplier(orbital_frequency_multiplier),
            __spin_frequency_multiplier(spin_frequency_multiplier),
            __phase_lag(phase_lag)
        {}

        ///\brief Change which term is dissipative.
        void setup(
            ///See SingleTermzone::SingleTermZone()
            int orbital_frequency_multiplier,

            ///See SingleTermzone::SingleTermZone()
            int spin_frequency_multiplier,

            ///See SingleTermzone::SingleTermZone()
            double phase_lag
        );

        ///\brief Return the tidal phase lag times the love number for the given
        ///tidal term (or one of its derivatives).
        ///
        ///In case the forcing frequency is exactly zero, return the phase lag
        ///for the case of the spin frequency approaching the term from below.
        ///The lag for spin frequency approaching from above is be written to
        ///above_lock_value. If the forcing frequency is non-zero,
        ///above_lock_value is not modified.
        virtual double modified_phase_lag(
            ///The multiplier of the orbital frequency in the
            ///expression for the forcing frequency.
            int orbital_frequency_multiplier,

            ///The multiplier of the spin frequency in the
            ///expression for the forcing frequency.
            int spin_frequency_multiplier,

            ///The current forcing frequency in rad/day.
            double forcing_frequency,

            ///The return value should be either the phase lag itself
            ///(NO_DERIV) or its derivative w.r.t. the specified quantity.
            Dissipation::QuantityEntry entry,

            ///If the lag of a locked term is calculated this should be set
            ///to the lag assuming the spin frequency is just above the lock.
            ///Otherwise, leave untouched.
            double &above_lock_value
        ) const;

        ///\brief Return the corresponding component of the love coefficient
        ///(Lai 2012 Equation 24).
        virtual double love_coefficient(
            ///The multiplier of the orbital frequency in the
            ///expression for the forcing frequency.
            int,

            ///The multiplier of the spin frequency in the
            ///expression for the forcing frequency.
            int,

            ///The return value should be either the phase lag itself
            ///(NO_DERIV) or its derivative w.r.t. the specified quantity.
            Dissipation::QuantityEntry
        ) const
        {return 0;}

        ///\brief Return true iff disspation has been defined for the zone.
        virtual bool dissipative() const
        {
            return (
                (
                    __orbital_frequency_multiplier != 0
                    ||
                    __spin_frequency_multiplier != 0
                )
                &&
                __phase_lag != 0
            );
        }

        virtual bool can_lock() const
        {
            return (
                __orbital_frequency_multiplier != 0
                &&
                __spin_frequency_multiplier != 0
                &&
                __phase_lag != 0
            );
        }
    }; //End SingleTermZone class.

} //End Evove namespace

#endif

