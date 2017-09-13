/**\file
 * 
 * \brief Definitions of some of the methods of ExponentialDecayDiffRotBody.
 *
 * \ingroup Star_group
 */

#include "ExponentialDecayDiffRotBody.h"

namespace Star {

    void ExponentialDecayDiffRotBody::reset_torque()
    {
        if(__torque.size()!=number_zones()-1)
            __torque.resize(number_zones()-1);
        for(unsigned i=0; i<number_zones()-1; ++i)
            __torque[i].resize(11, Eigen::Vector3d(Core::NaN,
                                                   Core::NaN,
                                                   Core::NaN));
    }

    Eigen::Vector3d &ExponentialDecayDiffRotBody::torque_entry(
            unsigned top_zone_index,
            Evolve::Dissipation::Derivative deriv,
            bool with_respect_to_top
    ) const
    {
        assert(deriv != Evolve::Dissipation::ORBITAL_FREQUENCY);
        assert(deriv != Evolve::Dissipation::ECCENTRICITY);
        assert(deriv != Evolve::Dissipation::SEMIMAJOR);
        assert(top_zone_index < number_zones()-1);

        std::valarray<Eigen::Vector3d> &zone_torque=__torque[top_zone_index];
        switch(deriv) {
            case Evolve::Dissipation::NO_DERIV :
                return zone_torque[0];
            case Evolve::Dissipation::SPIN_FREQUENCY : 
                return zone_torque[(with_respect_to_top ? 1 : 2)];
            case Evolve::Dissipation::INCLINATION :
                return zone_torque[(with_respect_to_top ? 3 : 4)];
            case Evolve::Dissipation::PERIAPSIS :
                return zone_torque[(with_respect_to_top ? 5 : 6)];
            case Evolve::Dissipation::MOMENT_OF_INERTIA :
                return zone_torque[(with_respect_to_top ? 7 : 8)];
            case Evolve::Dissipation::SPIN_ANGMOM :
                return zone_torque[(with_respect_to_top ? 9 : 10)];
            default:
                throw Core::Error::BadFunctionArguments(
                    "Unsupported derivative in "
                    "ExponentialDecayDiffRotBody::torque_entry"
                );
        }
    }

    void ExponentialDecayDiffRotBody::configure(bool initialize,
                                                double age,
                                                double companion_mass,
                                                double semimajor,
                                                double eccentricity,
                                                const double *spin_angmom,
                                                const double *inclination,
                                                const double *periapsis,
                                                bool locked_surface,
                                                bool zero_outer_inclination,
                                                bool zero_outer_periapsis)
    {
        if(age!=__current_age) {
            __current_age=age;
            reset_torque();
        }
        Evolve::DissipatingBody::configure(initialize,
                                           age,
                                           companion_mass,
                                           semimajor,
                                           eccentricity,
                                           spin_angmom,
                                           inclination,
                                           periapsis,
                                           locked_surface,
                                           zero_outer_inclination,
                                           zero_outer_periapsis);
    }

    Eigen::Vector3d ExponentialDecayDiffRotBody::angular_momentum_coupling(
        unsigned top_zone_index,
        Evolve::Dissipation::Derivative deriv,
        bool with_respect_to_top
    ) const
    {
        assert(top_zone_index<number_zones()-1);

        if(
            deriv == Evolve::Dissipation::ORBITAL_FREQUENCY
            || deriv == Evolve::Dissipation::ECCENTRICITY
            || deriv == Evolve::Dissipation::SEMIMAJOR
            || deriv == Evolve::Dissipation::AGE
            || deriv == Evolve::Dissipation::RADIUS
        )
            return Eigen::Vector3d(0, 0, 0);
        Eigen::Vector3d &result=torque_entry(top_zone_index,
                                             deriv,
                                             with_respect_to_top);
        if(!std::isnan(result[0])) return result;
        const Evolve::DissipatingZone &zone1 = zone(top_zone_index),
                                      &zone2 = zone(top_zone_index + 1);
        double i1 = zone1.moment_of_inertia(),
               i2 = zone2.moment_of_inertia();
        if(i1 == 0 || i2 == 0) {
            result.setZero();
            return result;
        }
        if(
            deriv == Evolve::Dissipation::INCLINATION
            ||
            deriv == Evolve::Dissipation::PERIAPSIS
        )
            result = zone_to_zone_transform(
                    zone2,
                    zone1,
                    Eigen::Vector3d(0, 0, zone2.spin_frequency()),
                    deriv,
                    !with_respect_to_top
            );
        else if(
            (
                deriv == Evolve::Dissipation::SPIN_FREQUENCY
                || 
                deriv==Evolve::Dissipation::SPIN_ANGMOM
            )
            &&
            !with_respect_to_top
        )
            result = Eigen::Vector3d(
                0,
                0,
                -(deriv == Evolve::Dissipation::SPIN_ANGMOM ? 1.0 / i1 : 1.0)
            );
        else {
            result = zone_to_zone_transform(zone2,
                                            zone1,
                                            Eigen::Vector3d(0, 0, 1));
            if(deriv == Evolve::Dissipation::SPIN_ANGMOM) result /= i2;
            else if(deriv != Evolve::Dissipation::SPIN_FREQUENCY)
                result *= zone2.spin_frequency();
            if(
                deriv == Evolve::Dissipation::NO_DERIV
                || 
                deriv == Evolve::Dissipation::MOMENT_OF_INERTIA
            )
                result[2] -= zone1.spin_frequency();
        }
        if(deriv != Evolve::Dissipation::MOMENT_OF_INERTIA) 
            result *= i1 * i2;
        else if(with_respect_to_top)
            result *= i2 * (1 - i1 / (i1 + i2));
        else result *= i1 * (1 - i2 / (i1 + i2));
        result/=__timescale * (i1 + i2);
        return result;
    }

}//End Star namespace.
