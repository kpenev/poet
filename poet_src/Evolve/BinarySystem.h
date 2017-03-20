/**\file
 *
 * \brief Defines the BinarySystem class.
 * 
 * \ingroup Evolve_group
 */

#ifndef __BINARY_SYSTEM_H
#define __BINARY_SYSTEM_H

#include "DissipatingBody.h"
#include "CombinedStoppingCondition.h"
#include "SecondaryDeathCondition.h"
#include "../Core/AstronomicalConstants.h"
#include "../Core/Common.h"
#include "../Core/OrbitalExpressions.h"
#include "../Core/Error.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_siman.h>
#include <string>
#include <limits>
#include <iostream>

namespace Evolve {

    ///\brief Describes a system of two bodies orbiting each other.
    ///
    ///The following variable are referred to throughout:
    /// - age: Gyr
    /// - a: semimajor axis in \f$R_\odot\f$ (average distance between the
    ///      bodies)
    /// - e: the eccentricity of the orbit 
    /// - Eorb: the energy of the orbit (potential + kinetic) treating the two
    ///         bodies as point masses in \f$GM_\odot^2/R_\odot\f$.
    /// - L: angular momentum of the orbit treating the two bodies as point
    ///      masses in
    ///      \f$M_\odot\cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$
    /// - \f$\theta^b_i\f$: inclination of the i-th zone of the b-th body
    ///					    relative to the orbit in radians.
    /// - \f$\omega^b_i\f$: argument of periapsis of the orbit for the i-th zone
    ///                     of the b-th body with a plane of reference
    ///                     perpendicular to that zone's spin.
    /// - \f$S^b_i\f$: spin angular momentum of the i-th zone of the b-th body in 
    ///             \f$M_\odot\cdot R_\odot^2 \cdot \mathrm{rad}/\mathrm{day}\f$.
    ///
    ///The index of the outermost zone of a body is zero and increases inward.
    ///
    ///All rates of change are per Gyr
    ///
    ///Positive tidal power and tidal torque mean that energy is being removed
    ///from the orbit and deposited into the bodies.
    ///
    ///\ingroup BinarySystem_group
    class BinarySystem {
    private:
        ///The name of the binary system (e.g. "HAT-P-20")
        std::string __name;

        ///\brief The evolution of the semimajor axis recorded by
        ///add_to_evolution() so far.
        std::list<double> __semimajor_evolution,

            ///\brief The evolution of the eccentricity recorded by
            ///add_to_evolution() so far
            __eccentricity_evolution;
        
        ///The present age of the stellar system in Gyrs.
        double __age,

               ///The current semimajor axis.
               __semimajor,

               ///The current eccentricity.
               __eccentricity,
               
               ///The current orbital energy.
               __orbital_energy,
               
               ///The current orbital angular momentum.
               __orbital_angmom,
               
               ///The rate at which the orbit gains energy due to tides.
               __orbit_energy_gain,
               
               ///\brief The rate at which the orbit gains angular momentum due
               ///to tides.
               __orbit_angmom_gain;

        ///The evolution mode from the last call to configure();
        Core::EvolModeType __evolution_mode;

        ///A list of indices of locked zones.
        std::list<unsigned> __locked_zones;

        ///\brief The above lock fractinos for the locked zones and their
        ///derivatives.
        ///
        ///For zone-dependent quantities, the derivatives are w.r.t. to the
        ///surface zone of __body1.
        std::valarray<Eigen::VectorXd> __above_lock_fractions,

            ///\brief The derivatives of the above lock fractions w.r.t. the
            ///inclinations of the zones.
            __above_lock_fractions_inclination_deriv,

            ///\brief The derivatives of the above lock fractions w.r.t. the
            ///periapses of the zones.
            __above_lock_fractions_periapsis_deriv,

            ///\brief The derivatives of the above lock fractions w.r.t. the
            ///moments of inertia of the zones.
            __above_lock_fractions_inertia_deriv,

            ///\brief The derivatives of the above lock fractions w.r.t. the
            ///angular momenta of the zones.
            __above_lock_fractions_angmom_deriv;

        ///\brief The derivative of the above lock fractions w.r.t. to the radius
        ///of the secondardy.
        ///
        ///The derivative w.r.t. to the radius of the primary is already
        ///available as __above_lock_fractions[Dissipation::RADIUS].
        Eigen::VectorXd __above_lock_fractions_body2_radius_deriv;

        ///\brief The matrix that defines the problem to solve for
        ///__above_lock_fractions.
        ///
        ///Useful when calculating derivatives.
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>
            __above_lock_fractions_decomp;

        DissipatingBody 
            ///\brief The first body in the system.
            ///
            ///This is the body that is allowed to exist both before and after
            ///the second body does. It is also the body that should support
            ///surface locks.
            &__body1, 
            
            ///\brief The second body in the system.
            ///
            ///This body is never allowed to exist alone in the system or to have
            ///its surface rotation locked.
            &__body2;

        ///Fills the __locked_zones list.
        void find_locked_zones();

        ///\brief Differential equations for the rotation of the zones of body 0
        ///with the topmost zone rotating with a fixed frequency.
        ///
        ///The age of the system must already be set appropriately by
        ///configure().
        int locked_surface_differential_equations(
                ///On output is set to the rates of change of \f$S^0_i\f$.
                double *evolution_rates) const;

        ///\brief Jacobian for the evolution of the rotation of the zones of 
        ///body 0 with the topmost zone rotating with a fixed frequency.
        ///
        ///The age of the system must already be set appropriately by
        ///configure().
        void locked_surface_jacobian(			
                ///On output is set to the Jacobian.
                double *param_derivs,
                
                ///On output is set to the partial age derivatives of the
                ///evolution equations.
                double *age_derivs) const;

        ///\brief Differential equations for the rotation of the zones of body 0
        ///if no other body is present.
        ///
        ///configure() must already have been called.
        int single_body_differential_equations(			
                ///On outputs is set to the rate of change of the orbital
                ///parameters.
                double *evolution_rates) const;

        ///Fills the jacobian for a system consisting of one isolated body.
        void fill_single_body_jacobian(
                ///The rows of the jacobian corresponding to the zone 
                ///inclinations.
                double *inclination_param_derivs,  

                ///The rows of the jacobian corresponding to the zone periapses.
                double *periapsis_param_derivs,

                ///The rows of the jacobian corresponding to the zone angular
                ///momenta.
                double *angmom_param_derivs, 

                ///The part of the age derivatives array corresponding to the
                ///zone inclinations.
                double *inclination_age_derivs, 

                ///The part of the age derivatives array corresponding to the
                ///zone periapses.
                double *periapsis_age_derivs, 

                ///The part of the age derivatives array corresponding to the
                ///zone angular momenta.
                double *angmom_age_derivs) const;

        ///\brief Jacobian for the evolution of the rotation of the zones of 
        ///body 0 if no other body is present.
        ///
        ///The age of the system must already be set appropriately by
        ///configure().
        void single_body_jacobian(
                ///On output is set to the Jacobian.
                double *param_derivs,
                
                ///On output is set to the partial age derivatives of the
                ///evolution equations.
                double *age_derivs) const;

        ///\brief Returns the rate of evolution of the semimajor axis or one of
        ///its derivatives.
        ///
        ///The age and the orbit of the system must already be set by calling
        ///configure().
        ///
        ///Supports differentiation with repesct to the semimajor axis. If
        ///derivatives with respect to other quantities are required, simply
        ///provide the derivative of the tidal power with respect to the quantity
        ///as the first argument.
        double semimajor_evolution(
                ///The rate at which the orbit gains energy (total for all zones 
                ///of all bodies) in
                /// \f$M_\odot R_\odot^2 \mathrm{day}^{-2}\mathrm{Gyr}^{-1}\f$
                double orbit_energy_gain,
                
                ///If not NaN, the derivative with respect to the semimajoir axis
                ///is returned, assuming that this is the derivative of
                ///orbit_energy_gain with respect to the semimajor axis.
                double orbit_energy_gain_deriv=Core::NaN) const;

        ///\brief Returns the rate of evolution of the eccentricity or one of its
        ///derivatives.
        ///
        ///The orbit and age of the system must already be set by calling
        ///configure().
        ///
        ///Supports differentiation with repesct to the semimajor axis and
        ///eccentricity. If derivatives with respect to other quantities are 
        ///required, simply provide the derivative of the tidal power and torque
        ///with respect to the quantity as the first and second arguments.
        double eccentricity_evolution(
                ///See semimajor_evolution()
                double orbit_energy_gain,

                ///The rate at which the orbit gains angular momentum (total for
                ///all zones of all bodies) in
                /// \f$M_\odot R_\odot^2 \mathrm{day}^{-2}\mathrm{Gyr}^{-1}\f$
                double orbit_angmom_gain,
                
                ///If this is not NaN, the derivative of the eccentricity
                ///evolution rate w.r.t. either the semimajor axis or the
                //eccentricity is returned instead of the rate itself. In this
                ///case, this value must be the derivative of orbit_power w.r.t.
                ///the same this as the desired derivative.
                double orbit_energy_gain_deriv=Core::NaN,
                
                ///If orbit_energy_gain_deriv is not NaN, this must be set to the
                ///derivative of orbit_torque with respect to the same variable
                ///as orbit_energy_gain_deriv.
                double orbit_angmom_gain_deriv=Core::NaN,
                
                ///If true the derivative calculated is assumed to be w.r.t. the
                ///semimajor axis.
                bool semimajor_deriv=true) const;

        ///\brief Makes corrections to the matrix and RHS to accomodate the given
        ///derivative for the linear problem that defines the above fractions.
        void above_lock_problem_deriv_correction(
                ///The derivative being calculated
                Dissipation::Derivative deriv,

                ///For zone-specific quantities, are we differentiating w.r.t. a
                ///zone part of body1?
                bool body1_deriv,
                
                ///The matrix to update.
                Eigen::MatrixXd &matrix, 

                ///The RHS vector to update.
                Eigen::VectorXd &rhs) const;

        ///\brief Calculates the fraction of a timestep above spin-orbit lock for
        ///all locked zones.
        ///
        ///The configure() method must already have been called.
        ///
        ///Only locked zones get an entry. The locked zones of body1 are first
        ///from outside to inside, followed by the locked zones for body2.
        void calculate_above_lock_fractions(
                ///The vector to fill with the calculated fractions.
                Eigen::VectorXd &fractions,
                
                ///The derivative of the above lock fractions to calculate. For
                ///zone-dependent quantities the derivatives are w.r.t. the
                ///surface zone quantity.
                Dissipation::Derivative deriv=Dissipation::NO_DERIV,

                ///Only used if deriv indicates a zone-dependent quantity. In
                ///that case, it specifies the body whose surface zone we are
                ///differentiating against.
                bool body1_deriv=true);


        ///\brief Calculates derivatives of the above lock fractions w.r.t.
        ///quantities of non-surface zones.
        Eigen::VectorXd above_lock_fractions_deriv(
                ///The derivative to calculate. Only INCLINATION, PERIAPSIS,
                ///MOMENT_OF_INERTIA and SPIN_ANGMOM derivatives are supported.
                Dissipation::Derivative deriv,

                ///The body whose zone's inclination we want the derivative 
                ///w.r.t.
                DissipatingBody &body,

                ///The zone whose inclination we are finding the derivative
                ///w.r.t.
                unsigned zone_index);

        ///Fills the __above_lock_fractinos_*_deriv members.
        void fill_above_lock_fractions_deriv();

        ///\brief Solves for and sets the above lock fractions and their
        ///derivatives.
        ///
        ///
        ///The system must already be configure() -ed in in BINARY evolution
        ///mode.
        void update_above_lock_fractions();

        ///\brief Fills an array with the evolution rates for the quantities
        ///describing a binary.
        ///
        ///The configure() method must already have been called.
        void fill_binary_evolution_rates(
                ///The torque on the orbit in the reference frame of the
                ///outermost zone of body 1.
                const Eigen::Vector3d &global_orbit_torque,

                ///The rate of change of the orbital parameters (see paramateres
                ///argument of binary_differential_equations().
                ///
                ///The size must be sufficient to hold all rates.
                double *evolution_rates) const;
                
        ///The differential equations for a system with both bodies present.
        int binary_differential_equations(
            ///On output is set to the rates of change of the evolution
            ///variables. See differintal_equations() for details.
            double *differential_equations);


        ///\brief Adds the derivatives of a rate by which the orbit is changing
        ///due to tidal interactions with a single body to an array.
        ///
        ///With the appropriate arguments this can either add the rate of change
        ///of the orbital energy or the torque on the orbit due to the tides on
        ///one body to the given array.
        template<typename VALUE_TYPE>
        void add_body_rate_deriv(
                ///The body whose contribution we wish to add.
                const DissipatingBody &body,

                ///The quantity we are collecting. Should be either
                ///body.tidal_orbit_energy_gain or body.tidal_orbit_torque.
                VALUE_TYPE (DissipatingBody::*func)(Dissipation::Derivative,
                                                unsigned,
                                                const Eigen::VectorXd &) const,
                
                ///The array to add to. The indices match the parameters being
                ///evolved.
                std::valarray<VALUE_TYPE> &orbit_rate_deriv,

                ///If body is __body1 this should be zero, if it is __body2, it
                ///should be the number of zones in __body1.
                unsigned offset) const;

        ///\brief Computes the derivatives w.r.t. the evolution quantities of the
        ///orbit energy gain.
        void fill_orbit_energy_gain_deriv(
                ///Location to fill.
                std::valarray<double> &orbit_energy_gain_deriv) const;

        ///\brief Computes the derivatives w.r.t. the evolution quantities of the
        ///orbit angular momentum gain.
        void fill_orbit_angmom_gain_deriv(
                ///Location to fill.
                std::valarray<double> &orbit_angmom_gain_deriv) const;

        ///Computes the row of the jacobian corresponding to the semimajor axis.
        void semimajor_jacobian(
                ///The derivatives of the orbit energy gain w.r.t. the evolution
                ///variables and age (last entry).
                const std::valarray<double> &orbit_energy_gain_deriv,

                ///Is the first variable \f$a^{6.5}\f$ instead of a?
                bool a6p5,
                
                ///The location to fill with the deriavtives of the semimajor
                ///evolution equation w.r.t. the evolution variables.
                double *param_derivs, 
                
                ///The location to set to the age derivative of the semimajor
                ///evolution equation.
                double &age_deriv) const;

        ///Computes the row of the jacobian corresponding to the eccentricity.
        void eccentricity_jacobian(
                ///The derivatives of the orbit energy gain w.r.t. the evolution
                ///variables and age (last entry).
                const std::valarray<double> &orbit_energy_gain_deriv,

                ///The derivatives of the orbit angular momentum gain w.r.t. the
                ///evolution variables and age (last entry).
                const std::valarray<double> &orbit_angmom_gain_deriv,

                ///Is the first variable \f$a^{6.5}\f$ instead of a?
                bool a6p5,
                
                ///The location to fill with the deriavtives of the eccentricity
                ///evolution equation w.r.t. the evolution variables.
                double *param_derivs, 
                
                ///The location to set to the age derivative of the eccentricity
                ///evolution equation.
                double &age_deriv) const;

        ///\brief Computes the partial age derivatives of the inclination or
        ///periapsis evolutiono equations.
        void angle_evolution_age_deriv(
                ///The body whose zone to calculate the age derivative for.
                DissipatingBody &body,

                ///The index of the zone in body to calculate the age derivative
                ///for.
                unsigned zone_ind,

                ///The sin of the inclination between the zone and the orbit.
                double sin_inc, 
                
                ///The cos of the inclination between the zone and the orbit.
                double cos_inc, 

                ///The index of the current zone in the list of locked zones
                ///(used only if zone is locked, ignored otherwise)
                unsigned locked_zone_ind,
                
                ///Destination overwritten by the age derivative of the
                ///inclination evolution rate.
                double &inclination,
                
                ///Destination overwritten by the age derivative of the periapsis
                ///evolution rate.
                double &periapsis) const;

        ///\brief Computes the partial derivatives w.r.t. semimamjor axis or
        ///eccentricity of the inclination and periapsis evolution equations.
        void angle_evolution_orbit_deriv(
                ///The derivative to compute, should be either
                ///Dissipation::SEMIMAJOR or Dissipation::ECCENTRICITY
                Dissipation::Derivative deriv, 

                ///The derivative of the orbital angular momentum w.r.t. deriv.
                double angmom_deriv,

                ///The body whose zone's evolution we are differentiating.
                DissipatingBody &body,

                ///The index of the zone whose evolution we are differentiating.
                unsigned zone_ind,

                ///The sin of the inclination between the zone and the orbit.
                double sin_inc, 
                
                ///The cos of the inclination between the zone and the orbit.
                double cos_inc, 

                ///The index of the current zone in the list of locked zones
                ///(used only if zone is locked, ignored otherwise)
                unsigned locked_zone_ind,
                
                ///Overwritten by the derivative of the inclination evolution
                ///rate.
                double &inclination,
                
                ///Overwritten by the derivative of the periapsis evolution rate.
                double &periapsis) const;

        void fill_orbit_torque_deriv(
                ///The derivatives to compute. Sholud be one of:
                ///Dissipation::INCLINATION, Dissipation::PERIAPSIS,
                ///Dissipation::MOMENT_OF_INERTIA, Dissipation::SPIN_ANGMOM
                Dissipation::Derivative deriv, 

                ///The body whose zone's coordinate system we are expressing the
                ///torque in.
                DissipatingBody &body,

                ///The index of the zone whose coordinate system we are
                ///expressing the torque in.
                unsigned zone_ind,

                ///On output gets filled with the derivatives of the torque on
                ///the orbit in the coordinate system of the zone identified by
                ///body and zone_ind w.r.t. deriv. Should be indexed by zone,
                ///with __body1 zones first.
                std::valarray<Eigen::Vector3d> &orbit_torque_deriv) const;

        ///\brief Calculates the derivatives of the torque on a zone w.r.t. a
        ///zone-specific quantity.
        void fill_zone_torque_deriv(
                ///The derivatives to compute. Sholud be one of:
                ///Dissipation::INCLINATION, Dissipation::PERIAPSIS,
                ///Dissipation::MOMENT_OF_INERTIA, Dissipation::SPIN_ANGMOM
                Dissipation::Derivative deriv, 

                ///The body whose zone's evolution we are differentiating.
                DissipatingBody &body,

                ///The index of the zone whose evolution we are differentiating.
                unsigned zone_ind,

                ///Overwritten with the result. If the zone is not locked this
                ///contains the derivatives w.r.t. to the quantity for the zone
                ///above, the zone and the zone below. If the zone is locked, the
                ///first three quantities are for the torque below the lock and
                ///an additional quantity is added at the end givin the derivative
                ///of the torque above the lock w.r.t. to the quantity of this
                ///zone.
                std::valarray<Eigen::Vector3d> &zone_torque_deriv) const;

        ///\brief Calculates the derivatives of an inclination evolution equation
        ///w.r.t. a zone specific quentity.
        void inclination_evolution_zone_derivs(
                ///The derivatives to compute. Sholud be one of:
                ///Dissipation::INCLINATION, Dissipation::PERIAPSIS,
                ///Dissipation::MOMENT_OF_INERTIA, Dissipation::SPIN_ANGMOM
                Dissipation::Derivative deriv, 

                ///The body whose zone's evolution we are differentiating.
                DissipatingBody &body,

                ///The index of the zone whose evolution we are differentiating.
                unsigned zone_ind,

                ///The x component of the torque on the zone (above a potential
                ///spin orbit lock). Ignored if the zone is not locked.
                double zone_x_torque_above, 

                ///The x component of the torque on the zone (below a potential
                ///spin orbit lock).
                double zone_x_torque_below,

                ///To be computed by fill_zone_torque_deriv().
                const std::valarray<Eigen::Vector3d> &zone_torque_deriv,

                ///The torque on the orbit in the coordinate system of the zone
                ///defined by the above 2 arguments.
                const Eigen::Vector3d &orbit_torque,

                ///The result of fill_orbit_torque_deriv.
                const std::valarray<Eigen::Vector3d> &orbit_torque_deriv,

                ///The derivatives of the above_lock fractions.
                const std::valarray<Eigen::VectorXd> &above_frac_deriv,

                ///The sin of the inclination between the zone and the orbit.
                double sin_inc, 
                
                ///The cos of the inclination between the zone and the orbit.
                double cos_inc, 

                ///The index of the current zone in the list of locked zones
                ///(used only if zone is locked, ignored otherwise)
                unsigned locked_zone_ind,

                ///Overwritten by the result.
                double *result) const;

        ///\brief Computes the derivatives of the periapsis evolution w.r.t. a
        ///zone specific quantity for all zones.
        void periapsis_evolution_zone_derivs(
                ///The derivative to compute.
                Dissipation::Derivative deriv, 

                ///The body whose zone's periapsis evolution to differentiate.
                DissipatingBody &body,

                ///The zone whose periapsis evolution to differentiate.
                unsigned zone_ind, 

                ///The y component of the torque on the zone above a lock if
                ///locked.
                double zone_y_torque_above, 

                ///The y component of the torque on the zone (assuming below a
                ///lock if locked).
                double zone_y_torque_below,

                ///The derivative of the zone torque.
                const std::valarray<Eigen::Vector3d> &zone_torque_deriv,

                ///The y torque on the orbit.
                double orbit_y_torque,

                ///The derivative of the orbit torque.
                const std::valarray<Eigen::Vector3d> &orbit_torque_deriv,

                ///The derivatives of the above lock fractions
                const std::valarray<Eigen::VectorXd> &above_frac_deriv,

                ///The sin of the inclination between the zone and the orbit.
                double sin_inc, 

                ///The cos of the inclination between the zone and the orbit.
                double cos_inc,

                ///The index in the list of locked zones of this zone if it is
                ///locked.
                unsigned locked_zone_ind, 

                ///Overwritten by the result. Should have the correct size.
                double *result) const;

        void spin_angmom_evolution_zone_derivs(
                ///The derivative to compute.
                Dissipation::Derivative deriv,

                ///The body whose zone's periapsis evolution to differentiate.
                DissipatingBody &body,

                ///The zone whose periapsis evolution to differentiate.
                unsigned zone_ind,

                ///The z component of the torque on the zone above a lock if
                ///locked.
                double zone_z_torque_above,

                ///The z component of the torque on the zone (assuming below a
                ///lock if locked).
                double zone_z_torque_below,

                ///The derivative of the zone torque.
                const std::valarray<Eigen::Vector3d> &zone_torque_deriv,

                ///The derivatives of the above lock fractions
                const std::valarray<Eigen::VectorXd> &above_frac_deriv,

                ///The index in the list of locked zones of this zone if it is
                ///locked.
                unsigned locked_zone_ind,

                ///Overwritten by the result. Should have the correct size.
                double *result) const;
                
        ///Calculates the jacobian for the evolution of a system of two bodies.
        void binary_jacobian(
                ///On output is set to the Jacobian.
                double *param_derivs,
                
                ///On output is set to the partial age derivatives of the
                ///evolution equations.
                double *age_derivs) const;

        ///Implements fill_orbit() for LOCKED_SURFACE_SPIN evolution mode.
        void fill_locked_surface_orbit(std::valarray<double> &orbit) const;

        ///Implements fill_orbit() for BINARY evolution mode.
        void fill_binary_orbit(std::valarray<double> &orbit) const;

        ///Implements fill_orbit() for SINGLE evolution mode.
        void fill_single_orbit(std::valarray<double> &orbit) const;

    public:
        ///Construct a binary system.
        BinarySystem(
                ///The first body in the system. Assumed to always be there, so
                ///for a star-planet system this should be the star.
                DissipatingBody &body1, 

                ///The second body in the system, initially may not be there and
                ///later may be engulfed by the first body.
                DissipatingBody &body2,

                ///The name of the system.
                const std::string &system_name="")
            : __name(system_name),
            __above_lock_fractions(Dissipation::NUM_DERIVATIVES), __body1(body1),
            __body2(body2) {}

        ///Returns the name of the system.
        const std::string get_name() const {return __name;}

        ///Sets the current state of the system.
        virtual int configure(
            ///Is this the first time configure() is invoked?
            bool initialize,

            ///The age to set the system to.
            double age,

            ///The semimajor axis of the orbit.
            double semimajor,

            ///The eccentricity of the orbit.
            double eccentricity,

            ///The spin angular momenta of the zones of the bodies (body 1 
            ///first, outermost zone to innermost, followed by body 2).
            const double *spin_angmom,

            ///The inclinations of the zones of the bodies (same order as
            ///spin_angmom). The surface zone inclination must be omitted for
            ///single body systems.
            const double *inclination,

            ///The arguments of periapsis of the zones of the bodies (same
            ///order as spin_angmom, but not including the surface zone of
            ///the first body).
            const double *periapsis,

            ///The evolution mode to assume.
            Core::EvolModeType evolution_mode
        );

        ///Sets the current state of the system directly from the evolutino
        ///variables.
        int configure(
            ///Is this the first time configure() is invoked?
            bool initialize,

            ///The current age of the system.
            double age,

            ///The parameters being evolved. See differential_equations().
            const double *parameters, 

            ///The evolution mode to assume for the system.
            Core::EvolModeType evolution_mode
        );

        ///Returns the present age of the system in Gyr.
        double age() const {return __age;}

        ///Returns the primary body in the system (const).
        const DissipatingBody &primary() const {return __body1;}

        ///Returns the secondary body in the system (const).
        const DissipatingBody &secondary() const {return __body2;}

        ///The total number of zones in both system bodies.
        unsigned number_zones() const
        {return (__evolution_mode==Core::BINARY
                 ? __body1.number_zones() + __body2.number_zones()
                 : __body1.number_zones());}

        ///How many zones on either body are currently locked.
        unsigned number_locked_zones() const
        {return __body1.number_locked_zones() + __body2.number_locked_zones();}

        ///The current semimajor axis of the system.
        double semimajor() const {return __semimajor;}

        ///\brief Fills an array with the parameters expected by 
        ///differential_equations() and jacobian(), returning the evolution mode.
        ///
        ///The system must be appropriately configure() -ed already.
        Core::EvolModeType fill_orbit(
                ///The orbit to fill (resized as necessary).
                std::valarray<double> &orbit) const;


        ///\brief The fraction of an infinitesimal timestep that a zone spends
        ///spinning faster than the lock it is in.
        double above_lock_fraction(
                ///The index within the list of locked zones of the zone we need.
                unsigned locked_zone_index,

                ///What derivative to return. Only non-zone specific quantities
                ///and, Dissipation::INCLINATION, Dissipation::PERIAPSIS,
                ///Dissipation::SPIN_ANGMOM and Dissipation::MOMENT_OF_INERTIA
                ///are allowed.
                Dissipation::Derivative deriv=Dissipation::NO_DERIV,

                ///The index of the zone whose quantity we want the derivative
                ///w.r.t. for zone specific quantities.
                unsigned deriv_zone_index=0,
                
                ///If deriv==Dissipation::RADIUS this argument determines the
                ///body whose radius the derivative is w.r.t.
                bool secondary_radius=false);

        ///\brief The differential equation and jacobian for the evolution of the
        ///system.
        ///
        ///Calls configure().
        int differential_equations(
                ///The current age of the system.
                double age,

                ///Contains the variables being evolved. The expected content
                ///depends on the values of evolution_mode and mhether any
                ///spin-orbit locks are held.
                ///
                ///If evolution_mode is LOCKED_SURFACE_SPIN:
                /// -# \f$S^0_i,\ i=1\ldots\f$ - all zones' spin axes are
                ///assumed aligned
                ///
                ///If evolution mode is SINGLE
                /// -# \f$theta^0_i,\ i=1\ldots\f$ (\f$\theta^0_0\equiv0\f$)
                /// -# \f$\omega^0_i,\ i=1\ldots\f$ (\f$\omega^0_0\equiv0\f$)
                /// -# \f$S^0_i,\ i=0,1,\ldots\f$
                ///
                ///The TABULATION evolution mode is not allowed.
                ///
                ///If evolution mode is BINARY
                /// -# a^6.5 if no zones are locked, a otherwise
                /// -# e
                /// -# \f$\theta^0_i,\ i=0,\ldots\f$
                /// -# \f$\theta^1_i,\ i=0,\ldots\f$
                /// -# \f$\omega^0_i,\ i=1,\ldots\f$ (\f$\omega^0_0\equiv0\f$)
                /// -# \f$\omega^1_i,\ i=0,\ldots\f$
                /// -# \f$S^0_i,\ i=0\ldots\f$ with zones in spin-orbit lock
                ///omitted
                /// -# \f$S^1_i,\ i=0\ldots\f$ with zones in spin-orbit lock
                ///omitted
                const double *parameters,

                ///The evolution mode to assume for the system.
                Core::EvolModeType evolution_mode,

                ///On outputs gets filled  with the rates at which the entries in
                ///parameters evolve. It is assumed that sufficient space has
                ///been allocated to hold the results.
                double *differential_equations);

        ///The jacobian of the evolution equations.
        ///
        ///Calls configure()!
        int jacobian(
                ///The age of the system.
                double age, 

                ///See differential_equations()
                const double *parameters, 

                ///The evolution mode to assume for the system.
                Core::EvolModeType evolution_mode,
                
                ///The matrix of partial derivatives of the evolution rates for
                ///the parameters w.r.t. each of the parameters.
                double *param_derivs,

                ///The partialderivatives of the evolution rates for the
                ///parameters w.r.t. age.
                double *age_derivs);
        
        ///\brief Check if a spin-orbit lock can be held and updates the system
        ///as necessary to calculate subsequent evolution.
        void check_for_lock(
                ///The multiplier of the orbital frequency at the lock
                int orbital_freq_mult, 
                
                ///The multiplier of the spin frequency of the potentially locked
                ///zone at the lock.
                int spin_freq_mult,

                ///The body whose zone is being locked.
                unsigned short body_index,

                ///The zone being unlocked.
                unsigned zone_index,
                
                ///The direction in which the lock will be left if it does not
                ///hold.
                short direction);

        ///\brief Smallest semimajor axis at which the secondary can survive for 
        ///the latest system configuration.
        ///
        ///By default returns the larger of:
        /// - primary radius.
        /// - 2.44*(secondary radius)*((primary mass)/(secondary mass))^(1/3)
        virtual double minimum_semimajor(
                ///If true the rate of change (per Gyr) is returned.
                bool deriv=false) const;

        ///The evolution mode of last call to configure().
        Core::EvolModeType evolution_mode() {return __evolution_mode;}

        ///\brief Update the system to account for the death of the secondary.
        ///
        ///Adds the entire orbital angular momentum to the surface zone of the
        ///primary and enters single body evolution mode.
        virtual void secondary_died();

        ///Releases the lock to one of the locked zones.
        virtual void release_lock(
                ///The index within the list of locked zones of the zone to
                ///unlock.
                unsigned locked_zone_index,
                
                ///The direction in which the zone's spin will evolve in the
                ///future relative to the lock.
                short direction);

        ///Appends the state defined by last configure(), to the evolution.
        virtual void add_to_evolution();

        ///Resets the evolution of the system.
        virtual void reset_evolution();

        ///Discards the last steps from the evolution.
        virtual void rewind_evolution(
                ///How many steps of evolution to discard.
                unsigned nsteps);

        ///\brief Conditions detecting the next possible doscontinuity in the
        ///evolution.
        ///
        ///Must be deleted when no longer necessary.
        virtual CombinedStoppingCondition *stopping_conditions();

        ///\brief Change the system as necessary at the given age.
        ///
        ///Handles things like the disk dissipating, the planet forming and
        ///interpolation discontinuities. 
        virtual void reached_critical_age(double age);

        ///\brief The next age when the evolution needs to be stopped for a
        ///system change
        virtual double next_stop_age() const;

        ///The tabulated evolution of the semimajor axis so far.
        const std::list<double> &semimajor_evolution() const
        {return __semimajor_evolution;}

        ///The tabulated evolution of the eccentricity so far.
        const std::list<double> &eccentricity_evolution() const
        {return __eccentricity_evolution;}

    }; //End BinarySystem class.

} //End Evolve namespace.

#endif
