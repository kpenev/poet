#ifndef __BROKEN_POWERLAW_PHASE_LAG_ZONE
#define __BROKEN_POWERLAW_PHASE_LAG_ZONE

/**\file
 *
 * \brief Declares the class that provides the phase lag function to
 * DissipatingZone objects.
 *
 * \ingroup StellarSystem_group
 */

#include "DissipatingZone.h"
#include "DissipatingBody.h"
#include "BinarySystem.h"
#include "CriticalForcingFrequencyCondition.h"
#include "CriticalSpinCondition.h"

namespace Evolve {

    class BrokenPowerlawPhaseLagZone : virtual public DissipatingZone {
    private:
        std::vector<double> 
            ///\brief The locations of the breaks in tidal frequency in rad/day.
            __tidal_frequency_breaks,

            ///The locations of the breaks in spin frequency in rad/day.
            __spin_frequency_breaks,
            
            ///\brief The powerlaw indices for the tidal frequency dependence.
            __tidal_frequency_powers,

            ///\brief The powerlaw indices for the spin frequency dependence.
            __spin_frequency_powers,
            
            ///\brief The phase lags at the tidal/spin frequency breaks.
            ///
            ///The tidal frequency break index changes faster and the spin
            ///frequency breaks index changes slower.
            __break_phase_lags;


        ///\brief The stopping condition monitoring for crossing of critical
        ///spin frequencies.
        CriticalSpinCondition *__spin_condition;

        ///Stoping conditions at the critical tidal frequencies for each
        ///order in eccentricity.
        std::list<CombinedStoppingCondition *> __tidal_frequency_conditions;

        ///\brief Make sure that the entries in __tidal_frequency_conditions
        ///are appropriate for the current eccentricity expansion order.
        void fill_tidal_frequency_conditions(
            ///The system being evolved.
            BinarySystem &system, 

            ///Is the body this zone is part of, the primary in the system.
            bool primary,

            ///The index of the zone in the body.
            unsigned zone_index
        );

        ///Print the configuration of the zone to stdlog.
        void print_configuration();

    public:
        ///\brief Create an unuseable zone. Must call setup() before use.
        BrokenPowerlawPhaseLagZone() :
            __spin_condition(NULL),
            __tidal_frequency_conditions()
        {}

        ///\brief Seup the zone with the given breaks/powers imposing
        ///continuity accress all breaks. Must only be called before use.
        void setup(
            ///The locations of the breaks in tidal frequency in rad/day.
            ///Entries should be sorted.
            std::vector<double> tidal_frequency_breaks,

            ///The locations of the breaks in spin frequency in rad/day.
            ///Entries should be sorted.
            std::vector<double> spin_frequency_breaks,
            
            ///The powerlaw indices for the tidal frequency dependence.
            ///Should be indexed in the same order as tidal_frequency_breaks,
            ///but must contain an additional starting entry for the powerlaw
            ///index before the first break.
            std::vector<double> tidal_frequency_powers,

            ///The powerlaw indices for the spin frequency dependence.
            ///Should be indexed in the same order as spin_frequency_breaks,
            ///but must contain an additional starting entry for the powerlaw
            ///index before the first break.
            std::vector<double> spin_frequency_powers,

            ///The phase lag at the first tidal and first spin frequency break.
            ///The rest are calculated by imposing continuity.
            double reference_phase_lag
        );
        
        ///\brief Should return the tidal phase lag times the love number for
        ///the given tidal term (or one of its derivatives).
        ///
        ///In case the forcing frequency is exactly zero, it should return
        ///the phase lag for the case of the spin frequency approaching the
        ///term from below. The lag for spin frequency approaching from above
        ///should be written to above_lock_value. If the forcing frequency is
        ///non-zero, leave above_lock_value untouched.
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
            Dissipation::Derivative deriv,

            ///If the lag of a locked term is calculated this should be set
            ///to the lag assuming the spin frequency is just above the lock.
            ///Otherwise, leave untouched.
            double &above_lock_value
        ) const;

        ///\brief Should return the corresponding component of the love
        ///coefficient (Lai 2012 Equation 24).
        virtual double love_coefficient(
            ///The multiplier of the orbital frequency in the
            ///expression for the forcing frequency.
            int,

            ///The multiplier of the spin frequency in the
            ///expression for the forcing frequency.
            int,

            ///The return value should be either the phase lag itself
            ///(NO_DERIV) or its derivative w.r.t. the specified quantity.
            Dissipation::Derivative
        ) const
        {return 0;}

        ///\brief Conditions detecting the next possible discontinuities in the
        ///evolution due to this zone.
        ///
        ///Must be deleted when no longer necessary.
        virtual CombinedStoppingCondition *stopping_conditions(
            ///The system being evolved.
            BinarySystem &system, 

            ///Is the body this zone is part of, the primary in the system.
            bool primary,

            ///The index of the zone in the body.
            unsigned zone_index
        );

        ///Changes the order of the eccentricity expansion performed.
        virtual void change_e_order(
            ///The new eccentricity expansion order.
            unsigned new_e_order,

            ///The system being evolved.
            BinarySystem &system, 

            ///Is the body this zone is part of, the primary in the system.
            bool primary,

            ///The index of the zone in the body.
            unsigned zone_index
        );

        ///Cleanup. 
        ~BrokenPowerlawPhaseLagZone();

    }; //End BrokenPowerlawPhaseLagZone class.

} //End BinarySystem namespace.

#endif
