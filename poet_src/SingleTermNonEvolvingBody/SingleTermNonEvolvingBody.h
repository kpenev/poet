/**\file
 *
 * \brief Declares a class for planets that are always locked to the orbit.
 *
 * \ingroup Planet_group
 */

#ifndef __SINGLE_TERM_NON_EVOLVING_BODY_H
#define __SINGLE_TERM_NON_EVOLVING_BODY_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SingleTermNonEvolvingZone.h"
#include "../Evolve/DissipatingBody.h"

namespace SingleTermNonEvolvingBody {

    ///\brief Single zone non-evolving planets with huge dissipation, so they
    ///always remain locked to the disk.
    ///
    ///\ingroup Planet_group
    class LIB_PUBLIC SingleTermNonEvolvingBody :
        virtual public Evolve::DissipatingBody {
    private:
        ///The only zone of the planet.
        SingleTermNonEvolvingZone __zone;
    public:
        ///Create a planet with a constant mass and radius.
        SingleTermNonEvolvingBody(double mass,
                                  double radius,
                                  double inertia_factor=0.3,
                                  int orbital_frequency_multiplier=0,
                                  int spin_frequency_multiplier=0,
                                  double phase_lag=0):
            __zone(mass,
                   radius,
                   inertia_factor,
                   orbital_frequency_multiplier,
                   spin_frequency_multiplier,
                   phase_lag)
        {};

        ///The number of zones the body consists of.
        unsigned number_zones() const {return 1;}

        ///Returns the only zone.
        const Evolve::DissipatingZone &zone(
            unsigned
#ifndef NDEBUG
            zone_index
#endif
        ) const
        {
            assert(zone_index == 0);
            return __zone;
        }

        ///Returns the only zone.
        Evolve::DissipatingZone &zone(
            unsigned
#ifndef NDEBUG
            zone_index
#endif
        )
        {
            assert(zone_index == 0);
            return __zone;
        }

        ///Returns the only zone.
        const SingleTermNonEvolvingZone &zone() const {return __zone;}

        ///Returns the only zone.
        SingleTermNonEvolvingZone &zone() {return __zone;}

        ///Should never be called.
        Eigen::Vector3d angular_momentum_coupling(
            unsigned,
            Evolve::Dissipation::QuantityEntry=Evolve::Dissipation::NO_DERIV,
            bool = false
        ) const
        {
            throw Core::Error::Runtime("Request for the angular momentum "
                                       "coupling of a SingleTermNonEvolving!");
        }

        ///Always zero.
        double angular_momentum_loss(
            Evolve::Dissipation::QuantityEntry=Evolve::Dissipation::NO_DERIV
        ) const
        {return 0;}

        ///No critical ages for non-evolving non-dissipating planets.
        void reached_critical_age(double) {}

    }; //End SingleTermNonEvolvingBody class.

}//End SingleTermNonEvolvingBody namespace.

#endif
