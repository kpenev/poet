/**\file
 *
 * \brief Definition of the Dissipation::Quantity and
 * Dissipation::Derivative output.
 *
 * \ingroup Evolve
 */

#define BUILDING_LIBRARY
#include "DissipationQuantities.h"

namespace Evolve {

    std::ostream &operator<<(std::ostream &os, 
            const Dissipation::Quantity &quantity)
    {
        switch(quantity) {
            case Dissipation::POWER : os << "POWER"; break;
            case Dissipation::TORQUEX : os << "TORQUEX"; break;
            case Dissipation::TORQUEZ : os << "TORQUEZ"; break;
            case Dissipation::SEMIMAJOR_DECAY : os << "SEMIMAJOR_DECAY"; break;
            case Dissipation::ORBIT_SPINUP : os << "ORBIT_SPINUP"; break;
            case Dissipation::INCLINATION_DECAY : os << "INCLINATION_DECAY";
                                                  break;
            default : assert(false);
        };
        return os;
    }

    ///More civilized output for Dissipation::Derivative variables.
    std::ostream &operator<<(std::ostream &os,
            const Dissipation::Derivative &deriv)
    {
        switch(deriv) {
            case Dissipation::NO_DERIV : os << "NO_DERIV";
                                         break;
            case Dissipation::AGE : os << "AGE";
                                    break;
            case Dissipation::SPIN_FREQUENCY : os << "SPIN_FREQUENCY";
                                               break;
            case Dissipation::ORBITAL_FREQUENCY: os << "ORBITAL_FREQUENCY";
                                                 break;
            case Dissipation::INCLINATION : os << "INCLINATION";
                                            break;
            case Dissipation::ECCENTRICITY: os << "ECCENTRICITY";
                                            break;
            case Dissipation::PERIAPSIS: os << "PERIAPSIS";
                                         break;
            case Dissipation::RADIUS : os << "RADIUS";
                                       break;
            case Dissipation::MOMENT_OF_INERTIA : os << "MOMENT_OF_INERTIA";
                                                  break;
            case Dissipation::SPIN_ANGMOM : os << "SPIN_ANGMOM";
                                            break;
            case Dissipation::SEMIMAJOR : os << "SEMIMAJOR";
                                          break;
            default : assert(false);
        };
        return os;
    }
} //End Evolve namespace.
