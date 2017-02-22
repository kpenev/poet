#include "BreakLockCondition.h"
#include "BinarySystem.h"

namespace Evolve {

    std::valarray<double> BreakLockCondition::operator()(
#ifdef DEBUG
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
#else
            Core::EvolModeType,
            const std::valarray<double> &,
#endif
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv) const
    {
#ifdef DEBUG
        assert(evol_mode==Core::BINARY);
        assert(orbit.size()==1 + 3*__system.number_zones() -
                             __system.number_locked_zones());
        assert(orbit.size()==derivatives.size());
        if(__system.number_locked_zones()==0)
            assert(std::pow(__system.semimajor(), 6.5)==orbit[0]);
        else assert(__system.semimajor()==orbit[0]);
#endif
        double frac=__system.above_lock_fraction(__locked_zone_index),
               dfrac_dt=__system.above_lock_fraction(__locked_zone_index,
                                                     Dissipation::AGE)
                        +
                        __system.above_lock_fraction(__locked_zone_index,
                                                     Dissipation::RADIUS, 0,
                                                     false)
                        *__system.primary().radius(1)
                        +
                        __system.above_lock_fraction(__locked_zone_index,
                                                     Dissipation::RADIUS, 0,
                                                     true)
                        *__system.secondary().radius(1);
        unsigned deriv_zone_ind=0,
                 unlocked_zone_ind=0,
                 num_zones=__system.number_zones();
        unsigned inclination_offset=2,
                   periapsis_offset=inclination_offset+num_zones-1,
                   angmom_offset=periapsis_offset+num_zones;
        for(int body_ind=0; body_ind<2; ++body_ind) {
            const DissipatingBody &body=(body_ind==0 ? __system.primary()
                                                     : __system.secondary());
            for(unsigned zone_ind=0; zone_ind<body.number_zones(); ++zone_ind) {
                double dangmom_dt;
                const DissipatingZone &zone=body.zone(zone_ind);
                if(!zone.locked())
                    dangmom_dt=derivatives[angmom_offset+unlocked_zone_ind++];
                else {
                    double dworb_dt=Core::orbital_angular_velocity(
                            __system.primary().mass(),
                            __system.secondary().mass(), __system.semimajor(),
                            true)
                        *derivatives[0];
                    dangmom_dt=zone.moment_of_inertia(1)*zone.spin_frequency()
                               +
                               zone.moment_of_inertia()
                               *zone.lock_held().spin(dworb_dt);
                }
                dfrac_dt+=__system.above_lock_fraction(
                        __locked_zone_index, Dissipation::MOMENT_OF_INERTIA,
                        deriv_zone_ind)*zone.moment_of_inertia(1)
                    +
                    __system.above_lock_fraction(
                            __locked_zone_index, Dissipation::INCLINATION,
                            deriv_zone_ind)
                    *derivatives[inclination_offset+deriv_zone_ind]
                    +
                    __system.above_lock_fraction(
                            __locked_zone_index, Dissipation::SPIN_ANGMOM,
                            deriv_zone_ind)
                    *derivatives[inclination_offset+deriv_zone_ind]
                    *dangmom_dt;
                if(deriv_zone_ind) 
                    dfrac_dt+=__system.above_lock_fraction(
                            __locked_zone_index, Dissipation::PERIAPSIS,
                            deriv_zone_ind)
                        *derivatives[periapsis_offset+deriv_zone_ind];
                ++deriv_zone_ind;
            }

        }
        dfrac_dt=Core::NaN;
        stop_deriv.resize(2, dfrac_dt);
        std::valarray<double> result(2);
        result[0]=frac;
        result[1]=frac-1.0;
        return result;
    }

    void BreakLockCondition::reached(short deriv_sign, unsigned index)
    {
#ifdef DEBUG
        if(index==0) assert(deriv_sign == -1);
        else {
            assert(index == 1);
            assert(deriv_sign == 1);
        }
#endif
        StoppingCondition::reached(deriv_sign, index);
        __system.release_lock(__locked_zone_index, deriv_sign);
    }

}//End Evolve namespace.
