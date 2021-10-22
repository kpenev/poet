/**\file
 *
 * \brief Declares the DissipatingBody class.
 *
 * \ingroup Evolve_group
 */


#ifndef __DISSIPATING_BODY_H
#define __DISSIPATING_BODY_H

#include "../Core/SharedLibraryExportMacros.h"

#include "DissipatingZone.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/AstronomicalConstants.h"
#include "../Core/Common.h"
#include <valarray>
#include <cassert>

namespace Evolve {

    class BinarySystem;

    ///\brief A base class for any body contributing to tidal dissipation.
    ///
    ///All torques are in units of
    /// \f$M_\odot R_\odot^2 \mathrm{rad}\mathrm{day}^{-1}\mathrm{Gyr}^{-1}\f$
    ///and tidal power is in units of
    /// \f$M_\odot R_\odot^2 \mathrm{rad}^2\mathrm{day}^{-2}\mathrm{Gyr}^{-1}\f$
    ///
    ///Age derivatives are per Gyr
    ///
    ///Masses and radii are in solar units.
    ///
    ///Angular velocities are in rad/day.
    ///
    ///\ingroup StellarSystem_group
    class LIB_PUBLIC DissipatingBody {
    private:
        double
            ///The coefficient used to normalize tidal power.
            __power_norm,

            ///The mean angular velocity with which the orbit is traversed.
            __orbital_frequency,

            ///The derivative of the orbital frequency w.r.t. semimajor axis.
            __dorbital_frequency_da,

            ///The frequency at which the surface is locked (if any).
            __surface_lock_frequency;

        std::valarray< std::valarray<Eigen::Vector3d> >
            ///\brief The rate of angular momentum transfer between two
            ///neighboring zones due to zone boundaries moving.
            ///
            ///The outer index is the index of the outer zone and each entry
            ///consists of two sub-entries giving the rate at which each of
            ///the two zones (outer first, inner second) involved gains
            ///angular momentum due to the moving boundary in each zone's
            ///coordinate system.
            __angular_momentum_transfer,

            ///\brief Tidal torques and their derivatives on each zone for
            ///spin frequency approaching potential lock from above.
            ///
            ///See tidal_torque() for more details.
            __tidal_torques_above,

            ///\brief Tidal torques and their derivatives on each zone for
            ///spin frequency approaching potential lock from below.
            ///
            ///See tidal_torque() for more details.
            __tidal_torques_below;

        ///\brief The quantities w.r.t. which derivatives of the orbit energy
        ///gain and torque are pre-calculated.
        std::vector<Dissipation::QuantityEntry> __orbit_entries;

        std::valarray<double>
            ///\brief Total tidal power (and derivatives and error) gained by
            ///the orbit from this body.
            ///
            ///The derivatives are only w.r.t. the non-zone specific
            ///quantities listed in __orbit_entries.
            __orbit_power,

            ///\brief Corrections to __orbit_power_below
            ///(undifferentiated) if single zones switch to above.
            __orbit_power_correction;

        ///\brief Torque on the orbit (and its derivatives and error) due to
        ///tides on this body.
        ///
        ///In the coordinate system of the topmost zone.
        ///
        ///The derivatives are only w.r.t. the non-zone specific quantities
        ///listed in __orbit_entries.
        std::vector<Eigen::Vector3d> __orbit_torque,

            ///\brief Corrections to __orbit_torque_below (undifferentiated)
            ///if single zones switch to above.
            ///
            ///It has already been transformed to the coordinate system of the
            ///surface zone.
            __orbit_torque_correction;

        ///The number of zones currently in a spin-orbit lock.
        unsigned __num_locked_zones;

        ///\brief The fractional contribution of the above the lock rates for
        ///locked zones and their derivatives.
        std::valarray<Eigen::VectorXd> __above_lock_fractions;

        ///\brief Scales the dimensionless torques as appropriate and corrects
        ///the relevant derivatives returning the normalization used.
        double normalize_torques(
            ///The mass of the other object in the system.
            double companion_mass,

            ///The orbital semimajor axis.
            double semimajor,

            ///The orbital frequency (passed to avoid duplicate calculation).
            double orbital_frequency
        );

        ///\brief Calcuates the total energy and angular momentum loss rates for
        ///the orbit due to this body's tidal dissipation.
        void collect_orbit_rates(
            ///The orbital frequency.
            double orbital_frequency,

            ///The normalization constant for the torques, as returned by
            ///normalize_torques().
            double torque_norm
        );

        ///\brief Calculates the correction the orbital energy gain and torque
        ///due to switching locked zones from fully below to fully above.
        void calculate_orbit_rate_corrections();

        ///\brief Angular momentum transfer between two neighboring zones due to
        ///moving zone boundaries.
        ///
        ///The argument of periapsis and the inclination must be defined for both
        ///zones involved.
        void angular_momentum_transfer(
            ///The outer zone.
            const DissipatingZone &outer_zone,

            ///The inner zone.
            const DissipatingZone &inner_zone,

            ///The angular momentum change for the outer zone in the outer
            ///zone's coordinate system.
            Eigen::Vector3d &outer_angmom_gain,

            ///The angular momentum change for the inner zone in the inner
            ///zone's coordinate system.
            Eigen::Vector3d &inner_angmom_gain,

            ///Derivatives with respect to inclination and periapsis can be
            ///computed, in addition to the actual transfer. It is an error
            ///to request another derivative.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///If deriv is not NO_DERIV, derivatives can be computed with
            ///respect to quantities of the outer zone (if this argument is
            ///true) or the inner zone (if false).
            bool with_respect_to_outer=false
        ) const;

        ///\brief Rate of angular momentum transfer (or its derivatives) to a
        ///zone due to its top boundary moving.
        ///
        ///The result is in the coordinate system of the zone.
        ///
        ///The set_orbit() method must already have been called and all zones
        ///must have their inclinations and arguments of periapsis set.
        Eigen::Vector3d angular_momentum_transfer_from_top(
            ///The index of the zone to calculate angular momentum gain for.
            ///Must not be zero.
            unsigned zone_index,

            ///Whether to return the quantity or one of its derivatives.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///If deriv is a zone specific quantity this argument determines
            ///if derivative with respect to the quantity of the zone above
            ///(true) or this zone (false) should be calculated.
            bool with_respect_to_outer=false
        ) const;

        ///\brief Rate of angular momentum transfer (or its derivatives) to a
        ///zone due to its bottom boundary moving.
        ///
        ///The result is in the coordinate system of the zone.
        ///
        ///The configure() method must already have been called.
        Eigen::Vector3d angular_momentum_transfer_from_bottom(
            ///The index of the zone to calculate angular momentum gain for.
            ///Must not be less than number_zones()-1.
            unsigned zone_index,

            ///Whether to return the quantity or one of its derivatives.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///If deriv is a zone specific quantity this argument determines
            ///if derivative with respect to the quantity of the zone below
            ///(true) or this zone (false) should be calculated.
            bool with_respect_to_inner=false
        ) const;

        ///\brief Rate of angular momentum transfer (or its derivatives) to a
        ///zone due to moving zone boundaries.
        ///
        ///The result is in the coordinate system of the zone.
        ///
        ///The configure() method must already have been called.
        Eigen::Vector3d angular_momentum_transfer_to_zone(
            ///The index of the zone to calculate angular momentum gain for.
            unsigned zone_index,

            ///Whether to return the quantity or one of its derivatives.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///See matching argument of external_torque() for description.
            int deriv_zone=0
        ) const;

        ///\brief Calculates the non-tidal torques on all zones.
        void calculate_nontidal_torques();

        ///\brief Corrects __orbit_power for zone locks.
        void correct_orbit_power(
            ///The derivative w.r.t. age of the above lock fractions.
            Eigen::VectorXd &above_lock_fractions_age_deriv,

            ///The derivative w.r.t. semimajor axis of the above lock
            ///fractions.
            Eigen::VectorXd &above_lock_fractions_semimajor_deriv,

            ///The derivative w.r.t. eccentricity of the above lock
            ///fractions.
            Eigen::VectorXd &above_lock_fractions_eccentricity_deriv,

            ///The derivative w.r.t. the body radius of the above lock
            ///fractions.
            Eigen::VectorXd &above_lock_fractions_radius_deriv
        );

        ///\brief Corrects __orbit_torque for zone locks.
        void correct_orbit_torque(
            ///The same as the argument of set_above_lock_fractions.
            std::valarray<Eigen::VectorXd> &above_lock_fractions
        );
    public:
        ///Some initializations for new objects.
        DissipatingBody();

        ///\brief Defines the orbit this body is in.
        ///
        ///The inclinations and arguments of periapsis must be already set for
        ///all zones.
        virtual void configure(
            ///Is this the first time the body is configure() -ed?
            bool initialize,

            ///The age to set the body to.
            double age,

            ///The mass of the second body in the system.
            double companion_mass,

            ///The semimajor axis of the orbit in \f$R_\odot\f$.
            double semimajor,

            ///The eccentricity of the orbit
            double eccentricity,

            ///The spin angular momenta of the non-locked zones of the body
            ///(outermost zone to innermost).
            const double *spin_angmom,

            ///The inclinations of the zones of the body (same order as
            ///spin_angmom). If NULL, all inclinations are assumed zero.
            const double *inclination = NULL,

            ///The arguments of periapsis of the zones of the bodies (same
            ///order as spin_angmom). If NULL, all periapses are assumed
            ///zero.
            const double *periapsis = NULL,

            ///If true, the outermost zone's spin is assumed locked to a
            ///disk and spin_angmom is assumed to start from the next zone.
            bool locked_surface = false,

            ///If true, the outermost zone's inclination is assumed to be
            ///zero and the inclination argument is assumed to start from the
            ///next zone.
            bool zero_outer_inclination = false,

            ///If true, the outermost zone's periapsis is assumed to be
            ///zero and the inclination argument is assumed to start from the
            ///next zone.
            bool zero_outer_periapsis = false
        );

        ///Locks the given zone's spin to the orbit with the given frequency
        ///ratio.
        void lock_zone_spin(unsigned zone_index,
                            int orbital_frequency_multiplier,
                            int spin_frequency_multiplier)
        {
            zone(zone_index).set_lock(orbital_frequency_multiplier,
                                      spin_frequency_multiplier);
            ++__num_locked_zones;
        }

        ///Releases the given zone from a spin-orbit lock.
        void unlock_zone_spin(
            ///The zone whose spin orbit lock to release.
            unsigned zone_index,

            ///The direction in which the spin will evolve in the future
            ///relative to the lock.
            short direction
        )
        {
            if(direction == 0) zone(zone_index).release_lock();
            else zone(zone_index).release_lock(direction);
            --__num_locked_zones;
        }

        ///The number of zones currently in a spin-orbit lock.
        unsigned number_locked_zones() const {return __num_locked_zones;}

        ///\brief External torque acting on a single zone (last
        ///calculate_torques_power()).
        ///
        ///This includes torques coupling it to other zones (due to differential
        ///rotation and/or due to mass transfer) and wind loss for the
        ///surface zone.
        ///
        ///For each zone the torques are in a coordinate system with
        /// \f$\hat{x}\f$ along the ascending node of the orbit in the zone's
        ///equatorial plane, and \f$\hat{z}\f$ along the zone's angular momentum.
        Eigen::Vector3d nontidal_torque(
            ///The index of the zone whose torque is needed.
            unsigned zone_index,

            ///Which entry to return for the torque.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///Since external torques depend on neighboring zones, this
            ///parameter is used to distinguish those (it is ignored if :
            ///deriv is Dissipation::NO_DERIV or a quantity which is not
            ///zone specific):
            /// - -1 Return the derivative with respect to the quantity for
            ///   the zone above.
            /// - 0 Return the derivative with respect to the quantity for
            ///   this
            ///   zone.
            /// - 1 Return the derivative with respect to the quantity for
            ///   the zone below.
            int deriv_zone=0
        ) const;

        ///\brief Tidal torque acting on the given zone (last
        ///calculate_torques_power()).
        ///
        ///The coordinate system is the same as for nontidal_torque().
        const Eigen::Vector3d &tidal_torque(
            ///The index of the zone whose tidal torque is required.
            unsigned zone_index,

            ///If a zone is currently in a lock, this argument decides
            ///whetherthe torque for spin frequency above (true) or below
            ///(false) of the lock should be returned.
            bool above,

            ///Which derivative of the tidal torque is required.
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const
        {
            assert(zone_index<number_zones());

            return (above
                    ? __tidal_torques_above
                    : __tidal_torques_below)[zone_index][entry];
        }

        ///\brief Tidal power dissipated in the given zone.
        double tidal_power(
            ///The index of the zone whose tidal torque is required.
            unsigned zone_index,

            ///If a zone is currently in a lock, this argument decides
            ///whetherthe torque for spin frequency above (true) or below
            ///(false) of the lock should be returned.
            bool above,

            ///Which derivative of the tidal power is required.
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV
        ) const;

        ///\brief Corrects the tidal orbit energy gain and angular momentum gain
        ///for locked zones.
        void set_above_lock_fractions(
            ///The fractional contributions of the above the lock rates for
            ///each zone (the vector) and their derivatives (the outer
            ///index).
            std::valarray<Eigen::VectorXd> &above_lock_fractions
        );

        ///\brief Rate of increase of the orbital energy due to tides in this
        ///body (last calculate_torques_power()).
        ///
        ///If set_above_lock_fractions() has not been called since cofigure(),
        ///returns the rate assuming all locked zones are fully below the lock.
        ///If it has, returns the actual rate.
        double tidal_orbit_power(
            ///W.r.t. that quantity is the required derivative.
            Dissipation::QuantityEntry entry=Dissipation::NO_DERIV,

            ///If deriv is a zone-specific quantity, this argument
            ///determines which zone to differentiate w.r.t.
            unsigned deriv_zone_index=0,

            ///The derivatives of all lock fractions w.r.t. to the quantity
            ///identified by the above arguments.
            const Eigen::VectorXd &
            above_lock_fraction_deriv=Eigen::VectorXd()
        ) const;

        ///\brief The torque on the orbit due to tidal dissipation in the body.
        ///
        ///If set_above_lock_fractions() has not been called since cofigure(),
        ///returns the torue assuming all locked zones are fully below the lock.
        ///If it has, returns the actual rate.
        Eigen::Vector3d tidal_orbit_torque(
            ///The derivative of the angular momentum rate required.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///If deriv is a zone-specific quantity, this argument
            ///determines which zone to differentiate w.r.t.
            unsigned deriv_zone_index=0,

            ///The derivatives of all lock fractions w.r.t. to the quantity
            ///identified by the above arguments.
            const Eigen::VectorXd &
            above_lock_fraction_deriv=Eigen::VectorXd()
        ) const;

        ///\brief Same as tidal_orbit_torque(Dissipation::QuantityEntry,
        ///unsigned, const Eigen::VectorXd &) but allow specifying the zone
        ///whose coordinate system to use.
        Eigen::Vector3d tidal_orbit_torque(
            ///The zone whose coordinate system to express the result.
            const DissipatingZone &reference_zone,

            ///The derivative of the angular momentum rate required.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///If deriv is a zone-specific quantity, this argument
            ///determines which zone to differentiate w.r.t.
            unsigned deriv_zone_index=0,

            ///The derivatives of all lock fractions w.r.t. to the quantity
            ///identified by the above arguments.
            const Eigen::VectorXd &above_lock_fraction_deriv=
            Eigen::VectorXd()
        ) const;

        ///The number of zones the body consists of.
        virtual unsigned number_zones() const =0;

        ///A modifiable reference to one of the body's zones.
        virtual const DissipatingZone &zone(
            ///The index of the zone within the body. Sequential zones are
            ///ordered from outside to inside (0 is the surface zone,
            ///number_zones()-1 is the core).
            unsigned zone_index
        ) const=0;

        ///A modifiable reference to one of the body's zones.
        virtual DissipatingZone &zone(
            ///The index of the zone within the body. Sequential zones are
            ///ordered from outside to inside (0 is the surface zone,
            ///number_zones()-1 is the core).
            unsigned zone_index
        )=0;

        ///\brief Coupling torque for two neighboring zones in the coordinate
        ///system of the top zone.
        ///
        ///Units:
        /// \f$M_\odot R_\odot^2\mathrm{rad}\mathrm{day}^{-1}mathrm{Gyr}^{-1}\f$.
        virtual Eigen::Vector3d angular_momentum_coupling(
            ///The index of the outer of the two zones whose coupling is
            ///needed.
            unsigned top_zone_index,

            ///The derivative of the coupling torque required.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///For derivatives with respect to zone specific quantities,
            ///this determines which zone's quantity to differentiate
            ///with respect to (top zone if true, bottom zone if false).
            bool with_respect_to_top=false
        ) const =0;

        ///\brief Rate of angular momentum loss by the top zone of the body
        ///and its derivatives.
        ///
        ///The spin frequency of the topmost zone of the body must already be
        ///set.
        virtual double angular_momentum_loss(
            ///The derivative of the angular momentum loss required.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV
        ) const =0;

        ///\brief The current radius or its derivative with age of the body.
        double radius(
            ///The order of the derivative to return.
            int deriv_order=0
        ) const
        {return zone(0).outer_radius(deriv_order);}

        ///The mass of the body (constant with age).
        double mass() const
        {return zone(0).outer_mass(Dissipation::NO_DERIV);}

        ///The surface spin freuqency of the body.
        double spin_frequency() const {return zone(0).spin_frequency();}

        ///\brief Angular velocity of the surface zone when locked (assumed
        ///constant).
        ///
        ///For example this could be the frequency of the stellar convective
        ///zone when locked to the disk.
        double surface_lock_frequency() const
        {return __surface_lock_frequency;}

        ///\brief Sets the frequency at which the surface is locked (if any).
        void set_surface_lock_frequency(double frequency)
        {__surface_lock_frequency=frequency;}

        ///Appends the state defined by last configure(), to the evolution.
        virtual void add_to_evolution();

        ///Discards the last steps from the evolution.
        virtual void rewind_evolution(
            ///How many steps of evolution to discard.
            unsigned nsteps
        );

        ///Discards all evolution.
        virtual void reset_evolution();

        ///\brief Conditions detecting the next possible discontinuities in
        ///the evolution due to this body.
        ///
        ///Must be deleted when no longer necessary.
        virtual CombinedStoppingCondition *stopping_conditions(
            ///The system being evolved.
            BinarySystem &system,

            ///Is the body the primary in the system.
            bool primary
        );

        ///Notifies the body that its spin just discontinously jumped.
        virtual void spin_jumped() {zone(0).spin_jumped();}

        ///\brief Change the body as necessary at the given age.
        ///
        ///Handles things like interpolation discontinuities.
        virtual void reached_critical_age(double) {assert(false);}

        ///\brief The next age when the evolution needs to be stopped for a
        ///change in one of the bodies.
        virtual double next_stop_age() const {return Core::Inf;}

        ///Change the eccentricity expansion order for all dissipative zones.
        virtual void change_expansion_order(
            ///The new eccentricity expansion order.
            unsigned new_expansion_order,

            ///The system being evolved.
            BinarySystem &system,

            ///Is the body the primary in the system.
            bool primary
        );

        ///Virtual destructor.
        virtual ~DissipatingBody() {}
    }; //End DissipatingBody class.

} //End Evolve namespace.
#endif
