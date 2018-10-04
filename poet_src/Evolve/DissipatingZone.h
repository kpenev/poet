/**\file
 *
 * \brief Declares a class representing one zone of a body dissipative to
 * tidal distortions.
 *
 * \ingroup Evolve_group
 */

#ifndef __DISSIPATING_ZONE_H
#define __DISSIPATING_ZONE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "TidalPotentialTerms.h"
#include "ZoneOrientation.h"
#include "DissipationQuantities.h"
#include "SpinOrbitLockInfo.h"
#include "CombinedStoppingCondition.h"
#include "BreakLockCondition.h"
#include "SynchronizedCondition.h"
#include "../Core/Common.h"
#include <valarray>

namespace Evolve {

    ///IDs for quantities saved as part of the evolution.
    enum ZoneEvolutionQuantities {
        ///Angular momentum of the zone.
        ANGULAR_MOMENTUM,

        ///Inclination of the zone.
        INCLINATION,

        ///Periapsis of the zone.
        PERIAPSIS,

        ///Moment of inertia of the zone.
        MOMENT_OF_INERTIA,

        ///Age derivative of MOMENT_OF_INERTIA.
        MOMENT_OF_INERTIA_FIRST_DERIV,

        ///Age second deriv of MOMENT_OF_INERTIA
        MOMENT_OF_INERTIA_SECOND_DERIV,

        ///Outer radius boundary of the zone
        OUTER_RADIUS,

        ///First age deriv of OUTER_RADIUS
        OUTER_RADIUS_FIRST_DERIV,

        ///Second age deriv of OUTER_RADIUS
        OUTER_RADIUS_SECOND_DERIV,

        ///Outer mass boundary of the zone.
        OUTER_MASS,

        ///First age derivative of OUTER_MASS.
        OUTER_MASS_DERIV,

        ///Number of real values evolution quantities.
        NUM_REAL_EVOL_QUANTITIES,

        ///Eccentricity expansion order.
        E_ORDER=NUM_REAL_EVOL_QUANTITIES,

        ///For locked zones this is the orbital frequency multiple of the lock.
        ORBITAL_FREQ_MULTIPLIER,

        ///\brief For locked zones this is the spin frequency multiple of the
        ///lock.
        ///
        ///A value of zero indicates no lock was held.
        SPIN_FREQ_MULTIPLIER,

        ///The total number of quantities whose evolution is tracked.
        NUM_EVOL_QUANTITIES
    }; //End ZoneEvolutionQuantities enumeration.

    ///More civilized output for EvolVarType variables.
    std::ostream &operator<<(std::ostream &os,
                             const ZoneEvolutionQuantities &evol_var);

    ///Needed to break circular dependency.
    class BinarySystem;

    ///The collection of the (m-1, m'), (m, m') and (m+1, m') coefficients
    ///needed for updating the tidal torque and power.
    class TidalTermTriplet {
    public:
        double
            ///The (m-1, m') coefficient.
            m_minus_one,

            ///The (m, m') coefficient.
            m,

            ///The (m+1, m') coefficient.
            m_plus_one;

        TidalTermTriplet(double m_minus_one_value = 0.0,
                         double m_value = 0.0,
                         double m_plus_one_value = 0.0):
            m_minus_one(m_minus_one_value),
            m(m_value),
            m_plus_one(m_plus_one_value)
        {}

        TidalTermTriplet &operator=(const TidalTermTriplet &rhs)
        {
            m_minus_one = rhs.m_minus_one;
            m = rhs.m;
            m_plus_one = rhs.m_plus_one;
            return *this;
        }
    };



    ///\brief A layer of a system body for which the tidal bulge is not exactly
    ///in phase with the tidal potential.
    class LIB_PUBLIC DissipatingZone : public ZoneOrientation {
    private:
        ///The expansion order in eccentricity to use.
        unsigned __e_order;

        ///The expansion of the tidal potential in series.
        TidalPotentialTerms __potential_term;

        ///\brief \f$\kappa_{m,m'}^+/\kappa_{m,m'}\f$ as a function
        ///of \f$m=-2 \ldots 2\f$.
        static const double __torque_x_plus_coef[];

        ///\brief \f$\kappa_{m,m'}^-/\kappa_{m,m'}\f$ as a function
        ///of \f$m=-2 \ldots 2\f$.
        static const double __torque_x_minus_coef[];

        double
            ///The current angular momentum of the zone.
            __angular_momentum,

            ///The current spin frequency of the zone.
            __spin_frequency,

            ///The orbital frequency (rad/day).
            __orbital_frequency,

            ///The absolute value of the angular momentum of the orbit
            __orbital_angmom;

        
        ///\brief The dimensionless tidal power and its error and derivatives.
        ///
        ///Consists of pairs of numbers one for each derivative. The first number
        ///of each pair is always filled and if the zone is in a lock it is the
        ///tidal power calculated assuming the zone spin frequency approaches the
        ///lock from below. The second number is filled only if the zone is in a
        ///spin-orbit lock and is the tidal power assuming the zone spin
        ///frequency approaches the lock from above. After all derivatives the
        ///final pair of numbers give the error in the undifferentiated value.
        std::valarray<double> __power,

            ///\brief The dimensionless tidal torque in the x direction and its
            ///derivatives.
            ///
            ///See description of __power for details on the content.
            __torque_x,

            ///\brief The dimensionless tidal torque in the y direction and its
            ///derivatives.
            ///
            ///See description of __power for details on the content.
            __torque_y,

            ///\brief The dimensionless tidal torque in the z direction and its
            ///derivatives.
            ///
            ///See description of __power for details on the content.
            __torque_z;

        ///\brief The lock the zone is currently held at (disabled if not locked).
        ///
        ///If the zone is not locked, this is one of the two terms closest to the
        ///current spin-orbit ratio and its sign is correct.
        SpinOrbitLockInfo __lock,

                          ///\brief The term closest matched by the current 
                          ///spin-orbit ratio in the other direction from __lock.
                          __other_lock;

        ///The floating point quantities whose evolution is tracked.
        std::vector< std::list<double> > __evolution_real;

        ///The integer quantities whose evolution is tracked.
        std::vector< std::list<int> > __evolution_integer;

        ///\brief If this zone is locked, this is its index in the list of locked
        ///zones in the system.
        unsigned __locked_zone_index;

        ///Is the zone in the middle of initializing (disable lock check)?
        bool __initializing;

        ///\brief Ensure the forcing frequency has the correct sign per the given
        ///constraint.
        ///
        ///If the forcing frequency has the opposite sign of what is expected
        ///based on the given lock information it is set to the smallest possible
        ///value with the correct sign.
        void fix_forcing_frequency(
                ///A tidal term for which the sign is known.
                const SpinOrbitLockInfo &limit, 

                ///The multiplier of the orbital frequency in the
                ///expression for the forcing frequency.
                int orbital_frequency_multiplier,

                ///The multiplier of the spin frequency in the
                ///expression for the forcing frequency.
                int spin_frequency_multiplier,

                ///The current forcing frequency to be updated if it violates the
                ///limit.
                double &forcing_frequency
        ) const;

        ///\brief Updates a SpinOrbitLockInfo variable as appropriate when 
        ///decreasing the eccentricity expansion order.
        ///
        ///__e_order must already be updated to the new value.
        void update_lock_to_lower_e_order(SpinOrbitLockInfo &lock);

        ///Updates __lock and __other_lock to accomodate increasing __e_order.
        void update_locks_to_higher_e_order(unsigned new_e_order);

        ///\brief Initializes the locks the first time the zone is
        ///configure() -ed.
        ///
        ///The spin frequency and orbital frequency must already be set.
        void initialize_locks();

        ///Add a term to the tidal torque and power arrays.
        void add_tidal_term(
            ///The azimuthal number of the term to add.
            int m,

            ///The time frequnecy number of the term to add.
            int mp,

            ///The tidal frequency of the term to add.
            double tidal_frequency,

            ///The \f$\mathcal{U}_{m-1, m'}\f$,
            /// \f$\mathcal{U}_{m, m'}\f$, and
            /// \f$\mathcal{U}_{m+1, m'}\f$ terms
            const TidalTermTriplet &U_value,

            ///The derivative with respect to inclination of U_value.
            const TidalTermTriplet &U_i_deriv,

            ///The derivative with respect to eccentricity of U_value.
            const TidalTermTriplet &U_e_deriv,

            ///Estimate of the error in U_value due to truncating the
            ///eccentricity expansion.
            const TidalTermTriplet &U_error
        );

#ifndef NDEBUG
        ///\brief Runs a bunch of asserts to check the consistency of __lock and
        ///__other_lock.
        virtual void check_locks_consistency() const;
#endif

    protected:
        ///Notify the zone that it is in the process of initializing or not.
        void initializing(bool flag) {__initializing = flag;}

        ///Configures only the spin of the zone.
        void configure_spin(
            ///See same name argument to configure().
            double spin,

            ///See same name argument to configure().
            bool spin_is_frequency
        );

    public:
        DissipatingZone();

        ///\brief Defines the current orbit, triggering re-calculation of all
        ///quantities.
        virtual void configure(
            ///Is this the first time the zone is configure() -ed?
            bool initialize,

            ///The age to set the zone to.
            double age,

            ///The angular velocity of the orbit in rad/day.
            double orbital_frequency,

            ///The eccentricity of the orbit
            double eccentricity,

            ///The absolute value of the angular momentum of the orbit.
            double orbital_angmom,

            ///The angular momentum or spin frequency of the zone
            ///if the zone is not in a spin--orbit lock (ignored it if is).
            double spin,

            ///The inclination of the zone relative to the orbit.
            double inclination,

            ///The argument of periapsis of the orbit in the equatorial
            ///planet of the zone.
            double periapsis,

            ///Should the spin argument be interpreted as an angular momentum
            ///or a spin frequency?
            bool spin_is_frequency
        );

        ///\brief The tidal forcing frequency for the given term and orbital
        ///frequency.
        ///
        ///Makes sure the forcing frequency has the correct sign regardless of
        ///numerical round-off.
        double forcing_frequency(
            ///The multiplier of the orbital frequency in the
            ///expression for the forcing frequency.
            int orbital_frequency_multiplier,

            ///The multiplier of the spin frequency in the
            ///expression for the forcing frequency.
            int spin_frequency_multiplier,

            ///The orbital frequency.
            double orbital_frequency
        ) const;

        ///\brief Defines the angular momentum of the reference zone for single
        ///body evolution.
        void set_reference_zone_angmom(double reference_angmom)
        {__orbital_angmom=reference_angmom;}

        ///\brief The rate at which the periapsis of the orbit/reference zone in
        ///this zone's equatorial plane is changing.
        ///
        ///Either configure() or set_reference_zone_angmom() must already have
        ///been called, and inclination() and spin_frequency() must be current.
        double periapsis_evolution(
            ///The torque on the orbit due to all other zones in this zone's
            ///coordinate system.
            const Eigen::Vector3d &orbit_torque,

            ///All torques acting on this zone (i.e. tidale, angular
            ///momentum loss due to wind for the surface zone and coupling to
            ///neightboring zones due to differential rotation).
            const Eigen::Vector3d &zone_torque,

            ///If not Dissipation::NO_DERIV, the corresponding entry of the rate
            ///is returned. For zone-specific quantities, derivative with
            ///respect to this zone's quantity is computed. If derivatives with
            ///respect to other zone's quantities are required, those only come
            ///in through the orbit torque and external torque, so pass the
            ///corresponding derivative instead of the actual torques, and
            ///ignore this and subsequent arguments.
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV,

            ///This argument is required if dervi is neither NO_DERIV nor
            ///PERIAPSIS, and should contain the derivative of the orbital
            ///torque relative to the quantity identified by deriv.
            const Eigen::Vector3d &orbit_torque_deriv=
            Eigen::Vector3d(),

            ///This argument is required if deriv is neither NO_DERIV nor
            ///PERIAPSIS, and shoul contain the derivative of the zone
            ///torque relative to the quantity identified by deriv.
            const Eigen::Vector3d &zone_torque_deriv=Eigen::Vector3d()
        );

        ///\brief The rate at which the inclination between this zone and the 
        ///orbit is changing.
        ///
        ///configure() must already have been called, and inclination() and
        ///spin_frequency() must be current.
        double inclination_evolution(
                ///The torque on the orbit or reference zone due all
                ///zones (including this one) in this zone's coordinate system.
                const Eigen::Vector3d &orbit_torque,

                ///All torques acting on this zone (i.e. tidal, angular
                ///momentum loss due to wind for the surface zone and coupling to
                ///neightboring zones due to differential rotation).
                const Eigen::Vector3d &zone_torque,
                
                ///If not Dissipation::NO_DERIV, the derivative of the rate with
                ///respect to the given quantity is returned. For zone-specific
                ///quantities, derivative with respect to this zone's quantity is
                ///computed. If derivatives with respect to other zone's
                ///quantities are required, those only come in through the orbit
                ///torque and zone torque, so pass the corresponding derivative
                ///instead of the actual torques, and ignore this and subsequent
                ///arguments.
                Dissipation::QuantityEntry entry=Dissipation::NO_DERIV,
                
                ///This argument is required if deriv is neither NO_DERIV nor
                ///PERIAPSIS, and should contain the derivative of the orbital
                ///torque relative to the quantity identified by deriv.
                const Eigen::Vector3d &orbit_torque_deriv=
                    Eigen::Vector3d(),

                ///This argument is required if dervi is neither NO_DERIV nor
                ///PERIAPSIS, and shoul contain the derivative of the zone
                ///torque relative to the quantity identified by deriv.
                const Eigen::Vector3d &zone_torque_deriv=
                    Eigen::Vector3d());


        ///\brief Should return true iff the given term is presently locked.
        virtual bool locked(int orbital_frequency_multiplier,
                            int spin_frequency_multiplier) const
        {return __lock(orbital_frequency_multiplier, 
                       spin_frequency_multiplier);}

        ///Should return true iff any tidal term is locked.
        virtual bool locked() const
        {return __lock;}

        ///The currntly held lock.
        const SpinOrbitLockInfo &lock_held() const {return __lock;}

        ///\brief Update the zone as necessary when the held lock disappears from
        ///the expansion.
        void release_lock();

        ///Update the zone as necessary when the held lock is broken.
        void release_lock(
            ///The direction that the spin will evolve toward in the future.
            short direction
        );

        ///Locks the zone spin to the orbit in the given ratio.
        void set_lock(int orbital_frequency_multiplier,
                      int spin_frequency_multiplier)
        {
            assert(!__lock);
            __lock.set_lock(orbital_frequency_multiplier,
                            spin_frequency_multiplier);
        }

        ///\brief Should return the tidal phase lag time the love number for the
        ///given tidal term (or one of its derivatives).
        ///
        ///In case the forcing frequency is exactly zero, it should return the
        ///phase lag for the case of the spin frequency approaching the term from
        ///below. The lag for spin frequency approaching from above should be
        ///written to above_lock_value. If the forcing frequency is non-zero, 
        ///leave above_lock_value untouched.
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
                ///Usually there will be no error due to truncating the
                ///eccentricity expansion coefficient.
                Dissipation::QuantityEntry entry,

                ///If the lag of a locked term is calculated this should be set
                ///to the lag assuming the spin frequency is just above the lock.
                ///Otherwise, leave untouched.
                double &above_lock_value) const =0;

        ///\brief Should return the corresponding component of the love
        ///coefficient (Lai 2012 Equation 24).
        virtual double love_coefficient(
            ///The multiplier of the orbital frequency in the
            ///expression for the forcing frequency.
            int orbital_frequency_multiplier,

            ///The multiplier of the spin frequency in the
            ///expression for the forcing frequency.
            int spin_frequency_multiplier,

            ///The return value should be either the phase lag itself
            ///(NO_DERIV) or its derivative w.r.t. the specified quantity.
            ///Usually there will be no error due to truncating the
            ///eccentricity expansion coefficient.
            Dissipation::QuantityEntry entry
        ) const =0;

        ///\brief Moment of inertia of the zone or its age derivative at the age
        ///of last configure() call.
        virtual double moment_of_inertia(
            ///What to return:
            /// - 0 The moment of inertia in \f$M_\odot R_\odot^2\f$
            /// - 1 The rate of change of the moment of inertia in 
            ///     \f$M_\odot R_\odot^2/Gyr\f$
            /// - 2 The second derivative in \f$M_\odot R_\odot^2/Gyr^2\f$
            int deriv_order=0
        ) const =0;

        ///\brief The moment of inertia of the zone or its age derivative at a
        ///specified age (no configure necessary).
        virtual double moment_of_inertia(
            ///The age at which to evaluate the moment of inertia.
            double age,

            ///What to return:
            /// - 0 The moment of inertia in \f$M_\odot R_\odot^2\f$
            /// - 1 The rate of change of the moment of inertia in 
            ///     \f$M_\odot R_\odot^2/Gyr\f$
            /// - 2 The second derivative in \f$M_\odot R_\odot^2/Gyr^2\f$
            int deriv_order=0
        ) const =0;

        ///\brief The spin frequency of the given zone.
        double spin_frequency() const {return __spin_frequency;}

        ///The angular momentum of the given zone in \f$M_\odot R_\odot^2\f$.
        double angular_momentum() const {return __angular_momentum;}

        ///\brief The dimensionless tidal power or one of its derivatives.
        double tidal_power(
            ///If a spin-orbit lock is in effect and the time-lag is
            ///discontinuous near zero forcing frequency, two possible values
            ///can be calculated, assuming that the zone spin frequency
            ///approaches the lock from below (false) or from above (true).
            bool above,

            ///What to return
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(entry <= Dissipation::END_DIMENSIONLESS_DERIV);
            assert(2 * entry + 1 < static_cast<int>(__power.size() - 2));

            return __power[2 * entry + (above? 1 : 0)];
        }

        ///\brief Same as tidal_power(bool, Dissipation::QuantityEntry), but
        ///using the predefined mix of below/above contributions.
        double tidal_power(
            ///The fraction of the timestep to assume to have spin above the
            ///lock.
            double above_fraction,

            ///What to return
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {

            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(!locked() || (above_fraction >= 0 && above_fraction <= 1));
            assert(entry <= Dissipation::END_DIMENSIONLESS_DERIV);
            assert(2 * entry + 1 < static_cast<int>(__power.size() - 2));

            return (
                above_fraction * __power[2 * entry + 1]
                +
                (1.0 - above_fraction) * __power[2 * entry]
            );
        }

        ///\brief The dimensionless tidal torque along x.
        ///
        ///See tidal_power() for a description of the arguments.
        double tidal_torque_x(
            bool above,
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(entry <= Dissipation::END_DIMENSIONLESS_DERIV);
            assert(2 * entry + 1 < static_cast<int>(__torque_x.size()));

            return __torque_x[2 * entry + (above ? 1 : 0)];
        }

        ///\brief Same as tidal_torque_x(bool, Dissipation::QuantityEntry) but
        ///below and above contributions mixed.
        double tidal_torque_x(
            ///The fraction of the timestep to assume to have spin above the
            ///lock.
            double above_fraction,

            ///The entry required (ignores the derivative of
            ///above_fraction if derivative is required).
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(
                !locked() || (above_fraction >= 0 && above_fraction <= 1)
            );
            assert(entry <= Dissipation::END_DIMENSIONLESS_DERIV);
            assert(
                2 * entry + 1
                <
                static_cast<int>(__torque_x.size())
            );

            return (above_fraction * __torque_x[2 * entry + 1]
                    +
                    (1.0 - above_fraction) * __torque_x[2 * entry]);
        }

        ///\brief The dimensionless torque along y.
        ///
        ///See tidal_power() for a description of the arguments.
        double tidal_torque_y(
            bool above,
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(entry <= Dissipation::END_DIMENSIONLESS_DERIV);
            assert(
                2 * entry + 1
                <
                static_cast<int>(__torque_y.size())
            );

            return __torque_y[2 * entry + (above? 1 : 0)];
        }

        ///\brief Same as tidal_torque_y(bool, Dissipation::QuantityEntry) but
        ///below and above contributions mixed.
        double tidal_torque_y(
            ///The fraction of the timestep to assume to have spin above the
            ///lock.
            double above_fraction,

            ///The torque entry required (ignores the derivative of
            ///above_fraction for derivatives).
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(!locked() || (above_fraction >= 0 && above_fraction <= 1));
            assert(entry <= Dissipation::END_DIMENSIONLESS_DERIV);
            assert(2 * entry + 1 < static_cast<int>(__torque_y.size()));

            return (above_fraction * __torque_y[2 * entry + 1]
                    +
                    (1.0 - above_fraction) * __torque_y[2 * entry]);
        }

        ///\brief The dimensionless tidal torque along z.
        ///
        ///See tidal_power() for a description of the arguments.
        double tidal_torque_z(
            bool above,
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(entry < Dissipation::END_DIMENSIONLESS_DERIV);
            assert(2 * entry + 1 < static_cast<int>(__torque_z.size()));

            return __torque_z[2 * entry + (above? 1 : 0)];
        }

        ///\brief Same as tidal_torque_z(bool, Dissipation::QuantityEntry) but
        ///below and above contributions mixed.
        double tidal_torque_z(
            ///The fraction of the timestep to assume to have spin above the
            ///lock.
            double above_fraction,

            ///The entry required (ignores the derivative of
            ///above_fraction for derivatives).
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            if(entry == Dissipation::EXPANSION_ERROR)
                entry = Dissipation::END_DIMENSIONLESS_DERIV;

            assert(!locked() || (above_fraction >= 0 && above_fraction <= 1));
            assert(entry < Dissipation::END_DIMENSIONLESS_DERIV);
            assert(2 * entry + 1 < static_cast<int>(__torque_z.size()));

            return (above_fraction * __torque_z[2 * entry + 1]
                    +
                    (1.0 - above_fraction) * __torque_z[2 * entry]);
        }

        ///\brief Outer radius of the zone or its derivative (per last
        //configure()).
        ///
        ///The outermost zone's outer radius is considered to be the radius of
        ///the body.
        virtual double outer_radius(
            ///What to return:
            /// - 0 The boundary in \f$R_\odot\f$
            /// - 1 The rate of change of the boundary in \f$R_\odot/Gyr\f$
            /// - 2 The second derivative in \f$R_\odot/Gyr^2\f$
            int deriv_order = 0
        ) const =0;

        ///\brief Same as outer_radius(int) but may be evaluated at a different
        ///age than for last confgure(). 
        virtual double outer_radius(double age, int deriv_order=0) const =0;

        ///\brief Mass coordinate of the zone's outer ouboundary or its
        ///derivative.
        ///
        ///The outermost zone's boundary is considered to be the mass of the
        ///body and should be constant.
        virtual double outer_mass(
            ///What to return:
            /// - 0 The boundary in \f$M_\odot\f$
            /// - 1 The rate of change of the boundary in \f$M_\odot/Gyr\f$
            /// - 2 The second derivative in \f$M_\odot/Gyr^2\f$
            int deriv_order=0
        ) const =0;

        ///\brief Same as outer_mass(int), but may be evaluated at a different
        ///age than last configure().
        virtual double outer_mass(double age, int deriv_order=0) const =0;

        ///To what order should eccentricity expansion be performed for the given
        ///value of the eccentricity.
        virtual unsigned eccentricity_order() const {return __e_order;}

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

        ///Appends the state defined by last configure(), to the evolution.
        virtual void add_to_evolution();

        ///Discards the last steps from the evolution.
        virtual void rewind_evolution(
            ///How many steps of evolution to discard.
            unsigned nsteps
        );

        ///Discards all evolution.
        virtual void reset_evolution();

        ///The tabulated evolution of a real valued quantity so far.
        const std::list<double> &get_evolution_real(
            ZoneEvolutionQuantities quantity
        ) const
        {
            assert(quantity < NUM_REAL_EVOL_QUANTITIES);

            return __evolution_real[quantity];
        }

        ///The tabulated evolution of an integer quantity so far.
        const std::list<int> &get_evolution_integer(
            ZoneEvolutionQuantities quantity
        ) const
        {
            assert(quantity >= NUM_REAL_EVOL_QUANTITIES);

            return __evolution_integer[quantity - NUM_REAL_EVOL_QUANTITIES];
        }


        ///\brief The index of this zone in the list of locked zones (valid only
        ///if locked).
        unsigned locked_zone_index() const
        {
            assert(__lock);

            return __locked_zone_index;
        }

        ///\brief Reference to the locked_zone_index() of this zone.
        unsigned &locked_zone_index()
        {
            assert(__lock);

            return __locked_zone_index;
        }

        ///\brief Should return true iff the zone has some non-zero dissipation.
        virtual bool dissipative() const =0;

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

        ///Notifies the zone that its spin just jumped discontinously.
        virtual void spin_jumped()
        {initialize_locks();}

        ///\brief Change the body as necessary at the given age.
        ///
        ///Handles things like interpolation discontinuities. 
        virtual void reached_critical_age(double)
        {assert(false);}

        ///\brief The next age when the evolution needs to be stopped for a
        ///change in one of the bodies.
        virtual double next_stop_age() const
        {return Core::Inf;}
    }; //End DissipatingZone class.

} //End Evolve namespace.

#endif
