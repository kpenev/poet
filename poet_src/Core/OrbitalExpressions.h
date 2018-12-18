/**\file
 *
 * \brief A collection of functions which calculate various quantities for
 * two body orbits.
 *
 */

#ifndef __ORBITAL_EXPRESSIONS_H
#define __ORBITAL_EXPRESSIONS_H

#include "../Core/SharedLibraryExportMacros.h"
#include <cmath>
#include "AstronomicalConstants.h"

namespace Core {

    ///Returns the orbital angular velocity of the given orbit in rad/day.
    LIB_PUBLIC double orbital_angular_velocity(
        ///The mass of the first body in \f$M_\odot\f$.
        double m1,

        ///The mass of the second body in \f$M_\odot\f$.
        double m2,

        ///The semimajor axis in \f$R_\odot\f$
        double semimajor,

        ///Whether to return the derivative with respect to the semimajor
        ///axis instead of the value.
        bool deriv=false
    );

    ///\brief The energy of the orbit (assuming 0 gravitational potential at
    ///infinity) in \f$\frac{M_\odot R_\odot^2 rad^2}{day^2}\f$.
    LIB_PUBLIC double orbital_energy(
        ///The mass of the first body in \f$M_\odot\f$.
        double m1,

        ///The mass of the second body in \f$M_\odot\f$.
        double m2,

        ///The semimajor axis in AU
        double semimajor,

        ///The order of the derivative w.r.t. a required.
        unsigned deriv_order=0
    );

    ///\brief The angular momentum of the orbit in 
    /// \f$\frac{M_\odot R_\odot^2 rad}{day}\f$.
    LIB_PUBLIC double orbital_angular_momentum(
        ///The mass of the first body in \f$M_\odot\f$.
        double m1,

        ///The mass of the second body in \f$M_\odot\f$.
        double m2,

        ///The semimajor axis in AU
        double semimajor,

        ///The eccentricity.
        double eccentricity
    );

    ///\brief Return the semiamjor axis in solar radii required to have the
    ///given masses orbit with the given period.
    LIB_PUBLIC double semimajor_from_period(
        ///The mass of the first body in \f$M_\odot\f$.
        double m1,

        ///The mass of the second body in \f$M_\odot\f$.
        double m2,

        ///The orbital period in days
        double period
    );

}//End Core namespace.

#endif
