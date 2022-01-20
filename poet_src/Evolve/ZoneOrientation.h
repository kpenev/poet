#ifndef __ZONE_ORIENTATION_H
#define __ZONE_ORIENTATION_H

/**\file
 *
 * \brief Declares a class for orientations of zones of DissipatingBody
 * objects.
 *
 * \ingroup StellarSystem_group
 */

/**\brief Orientations of zones of bodies in a binary system.
 *
 * \ingroup StellarSystem_group
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "DissipationQuantities.h"
#include "../Core/Common.h"
#include "../Core/IncludeEigen.h"

namespace Evolve {

    class LIB_PUBLIC ZoneOrientation {
    private:

        double
            ///The inclination of the zone relative to the orbit.
            __inclination,

            ///The rate at which the inclination of the zone is evolving.
            __inclination_rate,

            ///The periapsis of the orbit in the equatorial frame of the zone.
            __periapsis,

            ///The rate at which the periapsis of the zone is evolving.
            __periapsis_rate;
    public:
        ZoneOrientation(double inclination = Core::NaN,
                        double periapsis = Core::NaN) :
            __inclination(inclination),
            __periapsis(periapsis)
        {}

        ///Changes the zone orientation.
        void configure(double inclination, double periapsis)
        {__inclination = inclination; __periapsis = periapsis;}

        void set_evolution_rates(double inclination, double periapsis)
        {__inclination_rate = inclination; __periapsis_rate = periapsis;}

        ///The angle between the angular momenta of the zone and the orbit.
        double inclination(bool evolution_rate=false) const {
            return (evolution_rate ? __inclination_rate: __inclination);
        }

        ///The argument of periapsis of this zone minus the reference zone's
        double periapsis(bool evolution_rate=false) const {
            return (evolution_rate ? __periapsis_rate : __periapsis);
        }

    }; //End ZoneOrientation class.

    ///Transforms a vector betwen the coordinates systems of two zones.
    LIB_PUBLIC Eigen::Vector3d zone_to_zone_transform(
            ///The zone whose coordinate system the vectors are currently in.
            const ZoneOrientation &from_zone,

            ///The zone whose coordinate system we want to transform the vectors
            ///to.
            const ZoneOrientation &to_zone,

            ///The vector to transform.
            const Eigen::Vector3d &vector,

            ///Derivatives with respect to inclination and periapsis can be
            ///computed (assuming vector does not depend on these), in addition
            ///to the regular transform. It is an error to request another
            ///derivative.
            Dissipation::QuantityEntry deriv=Dissipation::NO_DERIV,

            ///If deriv is not NO_DERIV, derivatives can be computed with respect
            ///to quantities of the from_zone (if this argument is true) or the
            ///to_zone (if false).
            bool with_respect_to_from=false
    );

    /*
    ///Transforms the inclination and periapsis of a zone between references.
    void transform_zone_orientation(
            ///The zone whose orientation we wish to transform
            const ZoneOrientation &zone,

            ///The reference frame in which we want to express the zone's
            ///orientation in the old reference frame.
            const ZoneOrientation &reference,

            ///Overwritten by the inclination of zone in the new reference frame.
            double &inclination,

            ///Overwritten by the periapsis of zone in the new reference frame.
            double &periapsis
    );*/

}//End Evolve namespace.

#endif
