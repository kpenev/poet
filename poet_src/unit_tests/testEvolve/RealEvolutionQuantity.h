/**\file
 *
 * \brief Declare an inumeration for the various quantities tracked by the
 * evolution and checked.
 *
 * \ingroup UnitTests_group
 */

#ifndef __REAL_EVOLUTION_QUANTITY_H
#define __REAL_EVOLUTION_QUANTITY_H

namespace Evolve {

    ///Define identifiers for the quantities whose evolution we check.
    enum RealEvolutionQuantity {
        SEMIMAJOR,
        ECCENTRICITY,
        CONV_INCLINATION,
        RAD_INCLINATION,
        CONV_PERIAPSIS,
        RAD_PERIAPSIS,
        CONV_ANGMOM,
        RAD_ANGMOM,
        AGE,
        NUM_REAL_QUANTITIES
    };

}//End Evolve namespace.

#endif
