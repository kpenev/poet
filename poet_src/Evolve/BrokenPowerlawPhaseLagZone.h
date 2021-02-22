/**\file
 *
 * \brief Declares the class that provides the phase lag function to
 * DissipatingZone objects.
 *
 * \ingroup StellarSystem_group
 */

#ifndef __BROKEN_POWERLAW_PHASE_LAG_ZONE
#define __BROKEN_POWERLAW_PHASE_LAG_ZONE

#include "../Core/SharedLibraryExportMacros.h"
#include "DissipatingZone.h"
#include "DissipatingBody.h"
#include "BinarySystem.h"
#include "LagForcingFrequencyBreakCondition.h"
#include "LagSpinBreakCondition.h"


namespace Evolve {

    ///\brief A DissipatingZone where the phase lag is described by a broken
    ///powerlaw.
    ///
    ///By default the zone is non-dissipative.
    class LIB_PUBLIC BrokenPowerlawPhaseLagZone : virtual public DissipatingZone {
        friend class LagForcingFrequencyBreakCondition;
        friend class LagSpinBreakCondition;
    private:
        ///Is the zone dissipative.
        bool __dissipative;

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

        double
            ///See setup()
            __inertial_mode_enhancement,

            ///See setup()
            __inertial_mode_sharpness;

        ///\brief The index within __spin_frequency_powers of the powerlaw
        ///now in effect.
        std::vector<double>::size_type __spin_index;

        ///\brief The indices within __tidal_frequency_powers of the
        ///powerlaw now in effect for each active tidal term.
        ///
        ///The index of the (m' * Worb - m * Wspin) term is given by:
        /// * (2 + m) + 5 * m' if m' > 0
        /// * (2 - m) - 5 * m' if m' < 0
        ///This way changing the eccentricity expansion coefficient only
        ///affects the tail of the vector.
        std::vector< std::vector<double>::size_type > __tidal_indices;

        ///The index of the given tidal term within __tidal_indices.
        inline std::vector< std::vector<double>::size_type >::size_type
            tidal_term_index(
                ///The orbital frequency multiplier of the forcing term being
                ///monitored.
                int orbital_frequency_multiplier,

                ///The spin frequency multiplier of the forcing term being
                ///monitored.
                int spin_frequency_multiplier
            ) const
            {
                return (
                    2
                    +
                    (orbital_frequency_multiplier > 0 ? 1 : -1)
                    *
                    (
                        spin_frequency_multiplier
                        +
                        5 * orbital_frequency_multiplier
                    )
                );
            }

        ///Return the current orbital frequency of the given system.
        double get_orbital_frequency(const BinarySystem &system) const;

        ///\brief Calculate the factor by which dissipation is enhanced due to
        ///inertial modes or one of its logarithmic derivatives.
        double get_inertial_mode_factor(
            ///Absolute value of the frequency of the tidal term for which
            //inertial mode enhancement is required.
            double abs_forcing_frequency,

            ///The multiplier of the orbital frequency in the expression for the
            ///forcing frequency. Only used if derivative is requested.
            int orbital_frequency_multiplier,

            ///The multiplier of the spin frequency in the expression for the
            ///forcing frequency. Only used if derivative is requested.
            int spin_frequency_multiplier,

            ///If Dissipation::NO_DERIV, return the enhancement factor. If
            //Dissipation::SPIN_FREQUENCY or Dissipation::ORBITAL_FREQUENCY,
            //return the derivative of the log of the enhancement factor w.r.t.
            ///the specified quantity. If anything else, derivative is zero.
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const;

        ///Cleanup either during destruction or in preparation for changing
        ///dissipation.
        void reset();

        ///\brief Properly set the value of __spin_index per the current spin
        ///of the zone.
        ///
        ///Should only be called once, when starting the evolution in a
        ///two-body configuration. After that, this values should be updated
        ///by stopping conditions.
        void set_spin_index();

        ///\brief The index within __tidal_frequency_powers to use for the
        ///given forcing frequency.
        ///
        ///Should only be called when starting the evolution in a two-body
        ///configuration. After that, this values should be updated by
        ///stopping conditions.
        std::vector<double>::size_type get_tidal_index(
            ///The absolute value of the forcing frequency.
            double abs_forcing_frequency
        ) const;

        ///\brief Make sure that the entries in __tidal_frequency_conditions
        ///are appropriate for the current eccentricity expansion order.
        void add_tidal_frequency_conditions(
            ///The system being evolved.
            BinarySystem &system,

            ///Is the body this zone is part of, the primary in the system.
            bool primary,

            ///The index of the zone in the body.
            unsigned zone_index,

            ///The final stopping condition returned by stopping_conditions()
            ///The newly constructed conditions get added to this argument.
            CombinedStoppingCondition &result
        );

        ///Print the configuration of the zone to stdlog.
        void print_configuration(std::ostream &out_stream = std::clog);

    public:
        ///\brief Create a non-dissipative zone. Must call setup() if
        ///the zone is dissipative.
        BrokenPowerlawPhaseLagZone() : __dissipative(false) {}

        ///\brief Seup the zone with the given breaks/powers and inertial mode
        ///enhancement. Continuous accress all breaks. Must only be called
        ///before use.
        void setup(
            ///The locations of the breaks in tidal frequency in rad/day.
            ///Entries should be sorted.
            const std::vector<double> &tidal_frequency_breaks,

            ///The locations of the breaks in spin frequency in rad/day.
            ///Entries should be sorted.
            const std::vector<double> &spin_frequency_breaks,

            ///The powerlaw indices for the tidal frequency dependence.
            ///Should be indexed in the same order as tidal_frequency_breaks,
            ///but must contain an additional starting entry for the powerlaw
            ///index before the first break.
            const std::vector<double> &tidal_frequency_powers,

            ///The powerlaw indices for the spin frequency dependence.
            ///Should be indexed in the same order as spin_frequency_breaks,
            ///but must contain an additional starting entry for the powerlaw
            ///index before the first break.
            const std::vector<double> &spin_frequency_powers,

            ///The phase lag at the first tidal and first spin frequency
            ///break. The rest are calculated by imposing continuity.
            double reference_phase_lag,

            ///The factor by which dissipation is enhanced in the inertial mode
            ///range.
            double inertial_mode_enhancement = 1.0,

            ///A parameter controlling how sharp the transition between
            ///non-enhanced and enhanced dissipation is near the inertial mode
            ///boundary.
            ///
            ///The expression for the enhancement factor is:
            /// \f$ \max\left[\left|\frac{2\Omega_\star}{\Omega_{m,s}}\right|^\beta, \gamma\right] \f$
            ///Where \f$\beta\f$ is inertial_mode_sharpness, and
            /// \f$\gamma\f$ is inertial_mode_enhancement.
            double inertial_mode_sharpness = 10.0
        );

        ///See DissipatingZone::configure().
        virtual void configure(
            bool initialize,
            double age,
            double orbital_frequency,
            double eccentricity,
            double orbital_angmom,
            double spin,
            double inclination,
            double periapsis,
            bool spin_is_frequency
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
            Dissipation::QuantityEntry entry,

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
            Dissipation::QuantityEntry
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

        ///\brief Return true iff disspation has been defined for the zone.
        virtual bool dissipative() const
        {return __dissipative;}

        ///Cleanup.
        ~BrokenPowerlawPhaseLagZone() {
#ifndef NDEBUG
            std::cerr << "Destroying powerlaw lag zone:"
                      << "Tidal breaks:";
            for(
                std::vector<double>::const_iterator
                    br_i = __tidal_frequency_breaks.begin();
                br_i != __tidal_frequency_breaks.end();
                ++br_i
            )
                std::cerr << " " << *br_i;
            std::cerr << "Spin breaks: ";
            for(
                std::vector<double>::const_iterator
                    br_i = __spin_frequency_breaks.begin();
                br_i != __spin_frequency_breaks.end();
                ++br_i
            )
                std::cerr << " " << *br_i;
            std::cerr << std::endl;
#endif
            reset();
        }

    }; //End BrokenPowerlawPhaseLagZone class.

} //End BinarySystem namespace.

#endif
