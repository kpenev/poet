/**\file
 *
 * \brief Defines constant to serve as identifier of the various quantities
 * handled by the stellar evolution interpoltaion.
 *
 * \ingroup StellarEvolution_group.
 */

#ifndef __STELLAR_EVOLUTION_QUANTITIES_H
#define __STELLAR_EVOLUTION_QUANTITIES_H

#include <vector>
#include <string>

namespace StellarEvolution {

    ///Defines the quantities tracked by stellar evolution and their order.
    enum QuantityID {
        ///The radius of the star for in \f$R_\odot\f$.
        RADIUS = 0,

        ///\brief The convective zone moment of inertia in
        /// \f$M_\odot \cdot R_\odot^2\f$.
        ICONV,

        ///The luminosity in \f$L_\odot\f$.
        LUM,

        ///The index of the first core quantity.
        FIRST_CORE_QUANTITY,

        ///\brief The radiative zone moment of inertia in
        /// \f$M_\odot \cdot R_\odot^2\f$.
        IRAD = FIRST_CORE_QUANTITY,

        ///The radiative zone mass in \f$M_\odot\f$.
        MRAD,

        ///The convective-radiative boundary in \f$R_\odot\f$.
        RRAD,

        ///The number of stellar evolution quantities tracked.
        NUM_QUANTITIES
    };

    static const std::vector<std::string> QUANTITY_NAME {
        "R*",
        "Iconv",
        "L*",
        "Irad",
        "Mrad",
        "Rrad",
    };


}//End StellarEvolution namespace.

#endif
