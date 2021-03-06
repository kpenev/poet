/**\file
 *
 * \brief Declaration of enumerations of dissipation quantities and
 * derivatives.
 *
 * \ingroup Evolve
 */

#ifndef __DISSIPATION_QUANTITIES_H
#define __DISSIPATION_QUANTITIES_H

#include "../Core/SharedLibraryExportMacros.h"
#include <ostream>
#include <cassert>

namespace Evolve {

    ///\brief Isolates constants related to the tidal dissipation.
    ///
    ///\ingroup StellarSystem_group
    namespace Dissipation {

        ///The quantities which evolve due to tidal dissipation
        enum LIB_LOCAL Quantity {
            ///The rate at which energy is deposited into the body.
            ///Units: \f$\frac{M_\odot R_\odot^2 rad^2}{day^2\,Gyr}\f$
            POWER,

            ///The torque exerted on the body in the x direction.
            ///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
            TORQUEX,

            ///The torque exerted on the body in the y direction.
            ///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
            TORQUEY,

            ///The torque exerted on the body in the z direction.
            ///Units: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
            TORQUEZ,

            ///Minus the rate of change of the semimajor axis in AU/Gyr.
            SEMIMAJOR_DECAY,

            ///The rate of change of the orbital frequency in rad/day/Gyr.
            ORBIT_SPINUP,

            ///The rate of decrease of the the angle between the spin angular
            ///momentum of the body and the orbital angular momentum. In units of
            /// \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
            INCLINATION_DECAY,

            ///The rate at which the eccentricity decays per Gyr.
            ECCENTRICITY_DECAY,

            ///The total number of dissipation quantitise supported.
            NUM_QUANTITIES
        }; //End Quantity enumeration.

        ///All evolving quantities have a number of entries each (see below).
        enum LIB_LOCAL QuantityEntry {

            ///The quantity itself, undifferentiated.
            NO_DERIV,

            ///\brief The derivative w.r.t. age, excluding the dependence through
            ///the body's radius and the moments of inertia, but including all
            ///else.
            ///
            ///The full age derivative can be reconstructed later by adding the
            ///radius derivative times the rate of change of the radius and
            ///similarly for the moments of inertia of the zones involved.
            AGE,

            ///\brief The derivative w.r.t. the spin frequency of a dissipating
            ///zone.
            ///
            ///Holding the moment of inertia constant.
            SPIN_FREQUENCY,

            ///The derivative w.r.t. the orbital frequency.
            ORBITAL_FREQUENCY,

            ///The above derivatives exist for modified phase lags, below do not.
            END_PHASE_LAG_DERIV,

            ///The derivative w.r.t. the angle between the spin angular momentum
            ///of the body and the orbital angular momentum.
            INCLINATION = END_PHASE_LAG_DERIV,

            ///The derivative w.r.t. the eccentricity.
            ECCENTRICITY,

            ///All derivatives above are available for dimensionless powers and
            ///torques and those below are not.
            END_DIMENSIONLESS_DERIV,

            ///The derivative w.r.t. the argument of periapsis of the orbit in
            ///the equatorial plane a given zone.
            PERIAPSIS = END_DIMENSIONLESS_DERIV,

            ///The derivative w.r.t. the radius of the body in \f$R_\odot\f$.
            RADIUS,

            ///The derivative w.r.t. the moment of inertia of the zone in
            /// \f$M_\odot R_\odot^2\f$.
            ///
            ///Holding the angular momentum constant.
            MOMENT_OF_INERTIA,

            ///\brief The derivative w.r.t. the spin angular momentum in
            /// \f$M_\odot R_\odot^2 rad/day\f$.
            ///
            ///Holding the moment of inertia constant but not the spin frequency.
            SPIN_ANGMOM,

            ///The derivative w.r.t. the semimajor axis in AU.
            SEMIMAJOR,

            ///The total number of derivatives supported
            NUM_DERIVATIVES,

            ///The error due to truncating the eccentricity series to finite
            ///order.
            EXPANSION_ERROR = NUM_DERIVATIVES,

            ///The total number of Entries
            NUM_ENTRIES
        }; //End QuantityEntries enumeration.

    } //End Dissipation namespace.

    LIB_LOCAL inline bool zone_specific(Dissipation::QuantityEntry entry)
    {
        return (entry == Dissipation::SPIN_FREQUENCY
                ||
                entry == Dissipation::INCLINATION
                ||
                entry == Dissipation::PERIAPSIS
                ||
                entry == Dissipation::MOMENT_OF_INERTIA
                ||
                entry == Dissipation::SPIN_ANGMOM);
    }

    ///More civilized output for Dissipation::Quantity variables.
    LIB_LOCAL std::ostream &operator<<(std::ostream &os,
            const Dissipation::Quantity &quantity);

    ///More civilized output for Dissipation::QuantityError variables.
    LIB_LOCAL std::ostream &operator<<(std::ostream &os,
                                       Dissipation::QuantityEntry entry);

} //End Evolve namespace.

#endif
