/**\file
 * 
 * \brief Declares a class for planets that are always locked to the orbit.
 *
 * \ingroup Planet_group
 */

#ifndef __LOCKED_PLANET_H
#define __LOCKED_PLANET_H

#include "LockedPlanetZone.h"
#include "../Evolve/DissipatingBody.h"

namespace Planet {

    /**\brief Single zone non-evolving planets with huge dissipation, so they
     * always remain locked to the disk.
     */
    class LockedPlanet : virtual public Evolve::DissipatingBody {
    private:
        ///The only zone of the planet.
        LockedPlanetZone __zone;
    public:
        ///Create a planet with a constant mass and radius.
        LockedPlanet(double mass, double radius) : __zone(mass, radius) {};

        ///The number of zones the body consists of.
        unsigned number_zones() const {return 1;}

        ///Returns the only zone.
        const Evolve::DissipatingZone &zone(unsigned) const {return __zone;}

        ///Returns the only zone.
        Evolve::DissipatingZone &zone(unsigned) {return __zone;}

        ///Should never be called.
        Eigen::Vector3d angular_momentum_coupling(
            unsigned,
            Evolve::Dissipation::Derivative = Evolve::Dissipation::NO_DERIV,
            bool = false
        ) const
        {
            throw Core::Error::Runtime("Request for the angular momentum "
                                       "coupling of a LockedPlanet!");
        }

        ///Always zero.
        double angular_momentum_loss(
            Evolve::Dissipation::Derivative =Evolve::Dissipation::NO_DERIV
        ) const
        {return 0;}
    }; //End LockedPlanet class.

}//End Planet namespace.

#endif
