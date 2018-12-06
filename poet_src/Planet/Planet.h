/**\file
 *
 * \brief Declares a class for planets that are always locked to the orbit.
 *
 * \ingroup Planet_group
 */

#ifndef __LOCKED_PLANET_H
#define __LOCKED_PLANET_H

#include "../Core/SharedLibraryExportMacros.h"
#include "PlanetZone.h"
#include "../Evolve/DissipatingBody.h"

namespace Planet {

    ///\brief Single zone non-evolving planets with huge dissipation, so they
    ///always remain locked to the disk.
    ///
    ///\ingroup Planet_group
    class LIB_PUBLIC Planet : virtual public Evolve::DissipatingBody {
    private:
        ///The only zone of the planet.
        PlanetZone __zone;
    public:
        ///Create a planet with a constant mass and radius.
        Planet(double mass,
                     double radius,
                     double inertia_factor=1e-6) :
            __zone(mass, radius, inertia_factor)
        {};

        ///The number of zones the body consists of.
        unsigned number_zones() const {return 1;}

        ///Returns the only zone.
        const Evolve::DissipatingZone &zone(unsigned zone_index) const
        {
            assert(zone_index == 0);
            return __zone;
        }

        ///Returns the only zone.
        Evolve::DissipatingZone &zone(unsigned zone_index)
        {
            assert(zone_index == 0);
            return __zone;
        }

        ///Returns the only zone.
        const PlanetZone &zone() const {return __zone;}

        ///Returns the only zone.
        PlanetZone &zone() {return __zone;}

        ///Should never be called.
        Eigen::Vector3d angular_momentum_coupling(
            unsigned,
            Evolve::Dissipation::QuantityEntry=Evolve::Dissipation::NO_DERIV,
            bool = false
        ) const
        {
            throw Core::Error::Runtime("Request for the angular momentum "
                                       "coupling of a Planet!");
        }

        ///Always zero.
        double angular_momentum_loss(
            Evolve::Dissipation::QuantityEntry=Evolve::Dissipation::NO_DERIV
        ) const
        {return 0;}

        ///No critical ages for non-evolving non-dissipating planets.
        void reached_critical_age(double) {}

    }; //End LockedPlanet class.

}//End Planet namespace.

#endif
