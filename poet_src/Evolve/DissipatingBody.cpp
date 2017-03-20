#include "DissipatingBody.h"

namespace Evolve {

    double DissipatingBody::normalize_torques(double companion_mass,
            double semimajor, double orbital_frequency)
    {
        double torque_norm = (
            std::pow(companion_mass / std::pow(semimajor, 3), 2)
            *
            std::pow(radius(), 5)
            *
            Core::AstroConst::G
            *
            Core::AstroConst::day
            *
            Core::AstroConst::Gyr
            *
            Core::AstroConst::solar_mass
            /
            std::pow(Core::AstroConst::solar_radius, 3)
        );
        double dangular_velocity_da = -1.5 * orbital_frequency / semimajor;
        for(unsigned zone_index = 0; zone_index < number_zones(); ++zone_index) {
            bool above = false;
            const DissipatingZone &current_zone = zone(zone_index);
            do {
                std::valarray<Eigen::Vector3d> &tidal_torque=
                    (above ? __tidal_torques_above : __tidal_torques_below)
                    [zone_index];
                tidal_torque[Dissipation::ORBITAL_FREQUENCY] += (
                    4.0
                    /
                    orbital_frequency
                    *
                    tidal_torque[Dissipation::NO_DERIV]
                );
                tidal_torque[Dissipation::RADIUS] = (
                    5.0 / radius() * tidal_torque[Dissipation::NO_DERIV]
                );
                tidal_torque[Dissipation::MOMENT_OF_INERTIA] = (
                    -current_zone.spin_frequency()
                    /
                    current_zone.moment_of_inertia()
                    *
                    tidal_torque[Dissipation::SPIN_FREQUENCY]
                );
                tidal_torque[Dissipation::SPIN_ANGMOM] = (
                    tidal_torque[Dissipation::SPIN_FREQUENCY]
                    /
                    current_zone.moment_of_inertia()
                );
                tidal_torque[Dissipation::SEMIMAJOR] = (
                    dangular_velocity_da
                    *
                    tidal_torque[Dissipation::ORBITAL_FREQUENCY]
                );
                above = !above;
            } while (above);
        }

        assert(__tidal_torques_above.size() == __tidal_torques_below.size());
        for(unsigned i = 0; i < __tidal_torques_above.size(); ++i) {
            assert(__tidal_torques_above[i].size()
                   ==
                   __tidal_torques_below[i].size());
            for(unsigned j = 0; j < __tidal_torques_above[i].size(); ++j) {
                __tidal_torques_above[i][j] *= torque_norm;
                __tidal_torques_below[i][j] *= torque_norm;
            }
        }
        return torque_norm;
    }

    void DissipatingBody::collect_orbit_rates(double orbital_frequency,
            double torque_norm)
    {
        __power_norm=torque_norm*orbital_frequency;
        DissipatingZone &surface_zone=zone(0);
        unsigned invalid_ind=__orbit_deriv.size(),
                 no_deriv_ind=invalid_ind, 
                 orbital_freq_deriv_ind=invalid_ind,
                 radius_deriv_ind=invalid_ind,
                 semimajor_deriv_ind=invalid_ind;
        for(unsigned i=0; i<__orbit_deriv.size(); ++i)
            switch(__orbit_deriv[i]) {
                case Dissipation::NO_DERIV : no_deriv_ind=i; break;
                case Dissipation::ORBITAL_FREQUENCY : orbital_freq_deriv_ind=i;
                                                      break;
                case Dissipation::RADIUS : radius_deriv_ind=i; break;
                case Dissipation::SEMIMAJOR : semimajor_deriv_ind=i; break;
                default:;
            }
        __orbit_energy_gain=0;
        for(unsigned zone_index=0; zone_index<number_zones(); ++zone_index) {
            DissipatingZone &current_zone=zone(zone_index);
            std::valarray<Eigen::Vector3d> &tidal_torque=
                __tidal_torques_below[zone_index];
            for(unsigned deriv_ind=0; deriv_ind<__orbit_deriv.size()-1;
                    ++deriv_ind) {
                Dissipation::Derivative deriv=__orbit_deriv[deriv_ind];
                if(deriv<Dissipation::END_DIMENSIONLESS_DERIV)
                    __orbit_energy_gain[deriv_ind]-=
                        current_zone.tidal_power(false, deriv);
                if(zone_index) {
                    __orbit_torque[deriv]-=
                        zone_to_zone_transform(current_zone, surface_zone,
                                tidal_torque[deriv]);
                } else __orbit_torque[deriv]=-tidal_torque[deriv];
            }
            if(zone_index) {
                __orbit_torque[Dissipation::INCLINATION]-=
                    zone_to_zone_transform(current_zone, surface_zone,
                            tidal_torque[Dissipation::NO_DERIV],
                            Dissipation::INCLINATION, false);
                __orbit_torque[Dissipation::PERIAPSIS]-=
                    zone_to_zone_transform(current_zone, surface_zone,
                            tidal_torque[Dissipation::NO_DERIV],
                            Dissipation::PERIAPSIS, false);
            }
        }
        __orbit_energy_gain*=__power_norm;
        __orbit_energy_gain[orbital_freq_deriv_ind]+=
            5.0/orbital_frequency*__orbit_energy_gain[no_deriv_ind];
        __orbit_energy_gain[radius_deriv_ind]+=
            5.0/radius()*__orbit_energy_gain[no_deriv_ind];
        __orbit_energy_gain[semimajor_deriv_ind]=
            __dorbital_frequency_da*
            __orbit_energy_gain[orbital_freq_deriv_ind];
    }

    void DissipatingBody::calculate_orbit_rate_corrections()
    {
        __orbit_energy_gain_correction.resize(__num_locked_zones);
        __orbit_torque_correction.resize(__num_locked_zones);
        for(unsigned zone_index=0; zone_index<number_zones(); ++zone_index)
            if(zone(zone_index).locked()) {
                DissipatingZone &this_zone=zone(zone_index);
                __orbit_energy_gain_correction[this_zone.locked_zone_index()]=
                    tidal_power(zone_index, false)
                    -
                    tidal_power(zone_index, true);
                __orbit_torque_correction[this_zone.locked_zone_index()]=
                    zone_to_zone_transform(this_zone, zone(0),
                            tidal_torque(zone_index, false)
                            -
                            tidal_torque(zone_index, true));
            }
    }

    void DissipatingBody::angular_momentum_transfer(
            const DissipatingZone &outer_zone,
            const DissipatingZone &inner_zone,
            Eigen::Vector3d &outer_angmom_gain,
            Eigen::Vector3d &inner_angmom_gain,
            Dissipation::Derivative deriv, bool with_respect_to_outer) const
    {
#ifdef DEBUG
        assert(deriv==Dissipation::NO_DERIV || deriv==Dissipation::INCLINATION
                || deriv==Dissipation::PERIAPSIS);
#endif
        double dm_dt=inner_zone.outer_mass(1),
               lost_spin=(dm_dt>=0 ? outer_zone : inner_zone).spin_frequency(),
               angmom_transfer=-2.0/3.0*std::pow(inner_zone.outer_radius(), 2)*
                   lost_spin*dm_dt;
        if(dm_dt>0) {
            outer_angmom_gain[0]=outer_angmom_gain[1]=0;
            outer_angmom_gain[2]=
                (deriv==Dissipation::NO_DERIV ? angmom_transfer : 0);
            inner_angmom_gain=-zone_to_zone_transform(outer_zone, inner_zone,
                    Eigen::Vector3d(0, 0, angmom_transfer), deriv, 
                    with_respect_to_outer);
        } else {
            inner_angmom_gain[0]=inner_angmom_gain[1]=0;
            inner_angmom_gain[2]=
                (deriv==Dissipation::NO_DERIV ? -angmom_transfer : 0);
            outer_angmom_gain=zone_to_zone_transform(inner_zone, outer_zone,
                    Eigen::Vector3d(0, 0, angmom_transfer), deriv, 
                    !with_respect_to_outer);
        }
    }

    Eigen::Vector3d DissipatingBody::angular_momentum_transfer_from_top(
            unsigned zone_index, Dissipation::Derivative deriv, 
            bool with_respect_to_outer) const
    {
#ifdef DEBUG
        assert(zone_index>0);
#endif
        const DissipatingZone &this_zone=zone(zone_index),
                              &zone_above=zone(zone_index-1);
        double scaling=Core::NaN, dm_dt=this_zone.outer_mass(1);
        if(deriv==Dissipation::NO_DERIV) scaling=1;
        else if(deriv==Dissipation::AGE) {
            scaling=2.0*this_zone.outer_radius(1)/this_zone.outer_radius(0)
                    + this_zone.outer_mass(2)/dm_dt;
        } else if(deriv==Dissipation::SPIN_FREQUENCY || 
                  deriv==Dissipation::MOMENT_OF_INERTIA ||
                  deriv==Dissipation::SPIN_ANGMOM) {
            if((with_respect_to_outer && dm_dt<=0) ||
                    (!with_respect_to_outer && dm_dt>=0)) 
                return Eigen::Vector3d(0, 0, 0);
            else if(deriv==Dissipation::SPIN_FREQUENCY)
                scaling=1.0/(with_respect_to_outer ? zone_above.spin_frequency()
                                                   : this_zone.spin_frequency());
            else if(deriv==Dissipation::MOMENT_OF_INERTIA)
                scaling=-1.0/(with_respect_to_outer
                              ? zone_above.moment_of_inertia()
                              : this_zone.moment_of_inertia());
            else scaling=1.0/(with_respect_to_outer
                              ? zone_above.angular_momentum()
                              : this_zone.angular_momentum());
        } else if(deriv==Dissipation::INCLINATION || 
                  deriv==Dissipation::PERIAPSIS) {
            if(dm_dt<=0) return Eigen::Vector3d(0, 0, 0);
            else {
                Eigen::Vector3d dummy, result;
                angular_momentum_transfer(zone_above, this_zone, dummy, result,
                        deriv, with_respect_to_outer);
                return result;
            }
        }
#ifdef DEBUG
        else assert(false);
#endif
        return scaling*__angular_momentum_transfer[zone_index-1][1];
    }

    Eigen::Vector3d DissipatingBody::angular_momentum_transfer_from_bottom(
            unsigned zone_index, Dissipation::Derivative deriv, 
            bool with_respect_to_inner) const
    {
#ifdef DEBUG
        assert(zone_index<number_zones()-1);
#endif
        const DissipatingZone &this_zone=zone(zone_index),
                              &zone_below=zone(zone_index+1);
        double scaling=Core::NaN, dm_dt=zone_below.outer_mass(1);
        if(deriv==Dissipation::NO_DERIV) scaling=1;
        else if(deriv==Dissipation::AGE) {
            scaling=2.0*zone_below.outer_radius(1)/zone_below.outer_radius(0)
                    + zone_below.outer_mass(2)/dm_dt;
        } else if(deriv==Dissipation::SPIN_FREQUENCY || 
                  deriv==Dissipation::MOMENT_OF_INERTIA ||
                  deriv==Dissipation::SPIN_ANGMOM) {
            if((with_respect_to_inner && dm_dt>=0) ||
                    (!with_respect_to_inner && dm_dt<=0))
                return Eigen::Vector3d(0, 0, 0);
            else if(deriv==Dissipation::SPIN_FREQUENCY)
                scaling=1.0/(with_respect_to_inner ? zone_below.spin_frequency()
                                                   : this_zone.spin_frequency());
            else if(deriv==Dissipation::MOMENT_OF_INERTIA)
                scaling=-1.0/(with_respect_to_inner
                              ? zone_below.moment_of_inertia()
                              : this_zone.moment_of_inertia());
            else scaling=1.0/(with_respect_to_inner
                              ? zone_below.angular_momentum()
                              : this_zone.angular_momentum());
        } else if(deriv==Dissipation::INCLINATION ||
                  deriv==Dissipation::PERIAPSIS) {
            if(dm_dt>=0) return Eigen::Vector3d(0, 0, 0);
            else {
                Eigen::Vector3d dummy, result;
                angular_momentum_transfer(this_zone, zone_below, result, dummy,
                        deriv, !with_respect_to_inner);
            }
        }
#ifdef DEBUG
        else assert(false);
#endif
        return scaling*__angular_momentum_transfer[zone_index][0];
    }

    Eigen::Vector3d DissipatingBody::angular_momentum_transfer_to_zone(
                unsigned zone_index, Dissipation::Derivative deriv,
                int deriv_zone) const
    {
        if(deriv==Dissipation::ORBITAL_FREQUENCY
                || deriv==Dissipation::ECCENTRICITY
                || deriv==Dissipation::RADIUS
                || deriv==Dissipation::SEMIMAJOR) return Eigen::Vector3d(0,0,0);
        Eigen::Vector3d result(0, 0, 0);
        if(zone_index>0 && (deriv==Dissipation::NO_DERIV || 
                deriv==Dissipation::AGE || deriv_zone<=0)) 
            result=angular_momentum_transfer_from_top(zone_index, deriv, 
                    deriv_zone<0);
        if(zone_index<number_zones()-1 && (deriv==Dissipation::NO_DERIV || 
                deriv==Dissipation::AGE || deriv_zone>=0))
            result+=angular_momentum_transfer_from_bottom(zone_index, deriv, 
                    deriv_zone>0);
        return result;
    }

    DissipatingBody::DissipatingBody() : 
        __orbit_deriv(6),
        __orbit_energy_gain(__orbit_deriv.size()),
        __orbit_torque(Dissipation::NUM_DERIVATIVES),
        __num_locked_zones(0)
    {
        __orbit_deriv[0]=Dissipation::NO_DERIV;
        __orbit_deriv[1]=Dissipation::AGE;
        __orbit_deriv[2]=Dissipation::ORBITAL_FREQUENCY;
        __orbit_deriv[3]=Dissipation::ECCENTRICITY;
        __orbit_deriv[4]=Dissipation::RADIUS;
        __orbit_deriv[5]=Dissipation::SEMIMAJOR;
    }

    void DissipatingBody::correct_orbit_energy_gain(
            Eigen::VectorXd &above_lock_fractions_age_deriv,
            Eigen::VectorXd &above_lock_fractions_semimajor_deriv,
            Eigen::VectorXd &above_lock_fractions_eccentricity_deriv,
            Eigen::VectorXd &above_lock_fractions_radius_deriv)
    {
        for(unsigned zone_index=0; zone_index<number_zones(); ++zone_index) {
            if(zone(zone_index).locked()) {
                unsigned locked_zone_index=zone(zone_index).locked_zone_index();
                __orbit_energy_gain_correction[locked_zone_index]=
                    tidal_power(zone_index, false)
                    -
                    tidal_power(zone_index, true);
                for(unsigned deriv_ind=0; deriv_ind<__orbit_deriv.size();
                        ++deriv_ind) {
                    Dissipation::Derivative deriv=__orbit_deriv[deriv_ind];
                    __orbit_energy_gain[deriv_ind]+=
                        __above_lock_fractions[Dissipation::NO_DERIV]
                        [locked_zone_index]*(
                                tidal_power(zone_index, false, deriv)
                                -
                                tidal_power(zone_index, true, deriv));
                    if(deriv!=Dissipation::NO_DERIV) {
                        double frac_deriv=Core::NaN;
                        if(deriv==Dissipation::AGE) frac_deriv=
                            above_lock_fractions_age_deriv[locked_zone_index];
                        else if(deriv==Dissipation::ORBITAL_FREQUENCY)
                            frac_deriv=above_lock_fractions_semimajor_deriv
                                [locked_zone_index]/__dorbital_frequency_da;
                        else if(deriv==Dissipation::ECCENTRICITY) 
                            frac_deriv=above_lock_fractions_eccentricity_deriv
                                [locked_zone_index]/__dorbital_frequency_da;
                        else if(deriv==Dissipation::RADIUS) frac_deriv=
                            above_lock_fractions_radius_deriv[locked_zone_index];
                        else if(deriv==Dissipation::SEMIMAJOR) 
                            frac_deriv=above_lock_fractions_semimajor_deriv
                                [locked_zone_index];
#ifdef DEBUG
                        else assert(false);
#endif
                        __orbit_energy_gain[deriv_ind]+=frac_deriv*
                            __orbit_energy_gain_correction[locked_zone_index];
                    }
                }
            }
        }
    }

    void DissipatingBody::correct_orbit_torque(
                std::valarray<Eigen::VectorXd> &above_lock_fractions)
    {
        DissipatingZone &surface_zone=zone(0);
        for(unsigned zone_index=0; zone_index<number_zones(); ++zone_index) {
            if(zone(zone_index).locked()) {
                unsigned locked_zone_index=zone(zone_index).locked_zone_index();
                DissipatingZone &this_zone=zone(zone_index);
                for(int deriv=Dissipation::NO_DERIV;
                        deriv<Dissipation::NUM_DERIVATIVES; ++deriv) {
                    Eigen::Vector3d correction=
                        above_lock_fractions[Dissipation::NO_DERIV]
                        [locked_zone_index]
                        *(tidal_torque(zone_index, false,
                                    static_cast<Dissipation::Derivative>(deriv))
                          -
                          tidal_torque(zone_index, true,
                                    static_cast<Dissipation::Derivative>(deriv))
                         );
                    if(zone_index) {
                        __orbit_torque[deriv]+=
                            zone_to_zone_transform(this_zone, surface_zone,
                                    correction);
                        if(deriv==Dissipation::INCLINATION || 
                                deriv==Dissipation::PERIAPSIS)
                            __orbit_torque[deriv]+=zone_to_zone_transform(
                                    this_zone, surface_zone,
                                    above_lock_fractions[Dissipation::NO_DERIV]
                                                        [locked_zone_index]
                                    *
                                    __orbit_torque_correction[locked_zone_index],
                                    static_cast<Dissipation::Derivative>(deriv));
                    } else __orbit_torque[deriv]+=correction;
                    if(deriv!=Dissipation::NO_DERIV)
                        __orbit_torque[deriv]+=
                            above_lock_fractions[deriv][locked_zone_index]
                            *__orbit_torque_correction[locked_zone_index];
                }
            }
        }

    }

    void DissipatingBody::configure(bool initialize,
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
        double orbital_angmom=Core::orbital_angular_momentum(
            companion_mass,
            mass(),
            semimajor,
            eccentricity
        );
        __orbital_frequency=Core::orbital_angular_velocity(
            companion_mass,
            mass(),
            semimajor,
            false
        );
        __dorbital_frequency_da=Core::orbital_angular_velocity(
            companion_mass,
            mass(),
            semimajor,
            true
        );

        __tidal_torques_above.resize(number_zones());
        __tidal_torques_below.resize(number_zones());
        __angular_momentum_transfer.resize(number_zones()-1);
        unsigned angmom_offset=(locked_surface ? 1 : 0);
        for(unsigned zone_index=0; zone_index<number_zones(); ++zone_index) {
            DissipatingZone &current_zone=zone(zone_index);
            double zone_inclination, zone_periapsis, zone_spin;
            if(!inclination) zone_inclination=0;
            else if(zero_outer_inclination) 
                zone_inclination=(zone_index
                                  ? inclination[zone_index-1]
                                  : 0);
            else zone_inclination=inclination[zone_index];
            if(!periapsis) zone_periapsis=0;
            else if(zero_outer_periapsis)
                zone_periapsis=(zone_index ? periapsis[zone_index-1] : 0);
            else zone_periapsis=periapsis[zone_index];
            if(locked_surface && zone_index==0)
                zone_spin=surface_lock_frequency();
            else if(current_zone.locked()) {
                zone_spin=Core::NaN;
                ++angmom_offset;
            } else zone_spin=spin_angmom[zone_index-angmom_offset];
            current_zone.configure(initialize,
                                   age,
                                   __orbital_frequency,
                                   eccentricity,
                                   orbital_angmom,
                                   zone_spin,
                                   zone_inclination,
                                   zone_periapsis,
                                   locked_surface && zone_index==0);
        }
        for(unsigned zone_index=0; zone_index<number_zones(); ++zone_index) {
            DissipatingZone &current_zone=zone(zone_index);
            if(zone_index<number_zones()-1) {
                __angular_momentum_transfer[zone_index].resize(2);
                angular_momentum_transfer(current_zone, zone(zone_index+1),
                        __angular_momentum_transfer[zone_index][0],
                        __angular_momentum_transfer[zone_index][1]);
            }
            bool above=false;
            do {
                std::valarray<Eigen::Vector3d> &tidal_torque=
                    (above ? __tidal_torques_above : __tidal_torques_below)
                    [zone_index];
                tidal_torque.resize(Dissipation::NUM_DERIVATIVES);
                for(int deriv=Dissipation::NO_DERIV;
                    deriv<Dissipation::END_DIMENSIONLESS_DERIV; ++deriv) {
                    tidal_torque[deriv][0]=
                        current_zone.tidal_torque_x(above,
                                static_cast<Dissipation::Derivative>(deriv));
                    tidal_torque[deriv][1]=
                        current_zone.tidal_torque_y(above,
                                static_cast<Dissipation::Derivative>(deriv));
                    tidal_torque[deriv][2]=
                        current_zone.tidal_torque_z(above,
                                static_cast<Dissipation::Derivative>(deriv));
                }
                above=!above;
            } while(above);
        }
        collect_orbit_rates(__orbital_frequency, 
            normalize_torques(companion_mass,
                              semimajor,
                              __orbital_frequency));
        calculate_orbit_rate_corrections();
        __above_lock_fractions.resize(0);
    }

    Eigen::Vector3d DissipatingBody::nontidal_torque(unsigned zone_index,
            Dissipation::Derivative deriv, int deriv_zone) const
    {
#ifdef DEBUG
        assert(zone_index<number_zones());
        assert(static_cast<int>(zone_index)+deriv_zone>=0);
        assert(static_cast<int>(zone_index)+deriv_zone
               <
               static_cast<int>(number_zones()));
#endif
        const DissipatingZone &this_zone=zone(zone_index);
        Eigen::Vector3d result(0, 0, 0);
        if(zone_index==0 && deriv_zone==0)
            result[2]=-angular_momentum_loss(deriv);
        result+=angular_momentum_transfer_to_zone(zone_index, deriv, deriv_zone);
        if(zone_index<number_zones()-1 && 
                (!zone_specific(deriv) || deriv_zone>=0))
            result+=angular_momentum_coupling(zone_index, deriv, deriv_zone==0);
        if(zone_index>0 && (!zone_specific(deriv) || deriv_zone<=0)) {
            result-=zone_to_zone_transform(
                    zone(zone_index-1), this_zone,
                    angular_momentum_coupling(zone_index-1, deriv,
                                              deriv_zone<0));
            if(deriv==Dissipation::INCLINATION || deriv==Dissipation::PERIAPSIS)
                result-=zone_to_zone_transform(
                        zone(zone_index-1), this_zone,
                        angular_momentum_coupling(zone_index-1), deriv, 
                        deriv_zone<0);
        }
        return result;
    }

    double DissipatingBody::tidal_power(unsigned zone_index,
            bool above, Dissipation::Derivative deriv) const
    {
#ifdef DEBUG
        assert(zone_index<number_zones());
#endif
        const DissipatingZone &this_zone=zone(zone_index);
        double result=(deriv<Dissipation::END_DIMENSIONLESS_DERIV
                       ? __power_norm*this_zone.tidal_power(above, deriv)
                       : 0.0);
        if(deriv==Dissipation::ORBITAL_FREQUENCY || 
                deriv==Dissipation::SEMIMAJOR) {
            result+=5.0/__orbital_frequency*__power_norm*
                this_zone.tidal_power(above);
            if(deriv==Dissipation::SEMIMAJOR) result*=__dorbital_frequency_da;
        }
        else if(deriv==Dissipation::RADIUS)
            result+=5.0/radius()*this_zone.tidal_power(above);
        return result;
    }

    void DissipatingBody::set_above_lock_fractions(
            std::valarray<Eigen::VectorXd> &above_lock_fractions)
    {
#ifdef DEBUG
        assert(above_lock_fractions.size()==Dissipation::NUM_DERIVATIVES);
#endif
        __above_lock_fractions.resize(Dissipation::NUM_DERIVATIVES);
        for(unsigned i=0; i<Dissipation::NUM_DERIVATIVES; ++i) {
#ifdef DEBUG
            assert(above_lock_fractions[i].size()>=__num_locked_zones);
#endif
            __above_lock_fractions[i].resize(above_lock_fractions[i].size());
        }
        __above_lock_fractions=above_lock_fractions;
        correct_orbit_energy_gain(above_lock_fractions[Dissipation::AGE],
                above_lock_fractions[Dissipation::SEMIMAJOR],
                above_lock_fractions[Dissipation::ECCENTRICITY],
                above_lock_fractions[Dissipation::RADIUS]);
        correct_orbit_torque(above_lock_fractions);
    }

    double DissipatingBody::tidal_orbit_energy_gain(
            Dissipation::Derivative deriv, unsigned deriv_zone_index,
            const Eigen::VectorXd &above_lock_fraction_deriv) const
    {
        const DissipatingZone &deriv_zone=zone(deriv_zone_index);
        double result=Core::NaN;
        if(!zone_specific(deriv) || !deriv_zone.locked() 
                || __above_lock_fractions.size()==0) {
            if(zone_specific(deriv))
                result=tidal_power(deriv_zone_index, false, deriv);
            else {
                for(unsigned i=0; i<__orbit_deriv.size(); ++i)
                    if(deriv==__orbit_deriv[i]) 
                        return __orbit_energy_gain[i];
#ifdef DEBUG
                assert(false);
#endif
            }
            if(__above_lock_fractions.size()==0) return result;
        } else {
            double above_frac=
                __above_lock_fractions[Dissipation::NO_DERIV]
                                    [zone(deriv_zone_index).locked_zone_index()];
            result=above_frac
                   *tidal_power(deriv_zone_index, true, deriv)
                   +
                   (1.0-above_frac)
                   *tidal_power(deriv_zone_index, false, deriv);
        }
#ifdef DEBUG
        assert(above_lock_fraction_deriv.size()==
                __above_lock_fractions[Dissipation::NO_DERIV].size());
#endif
        for(unsigned i=0; i<above_lock_fraction_deriv.size(); ++i)
            result+=
                above_lock_fraction_deriv[i]*__orbit_energy_gain_correction[i];
        return result;
    }

    Eigen::Vector3d DissipatingBody::tidal_orbit_torque(
            Dissipation::Derivative deriv, unsigned deriv_zone_index,
            const Eigen::VectorXd &above_lock_fraction_deriv) const
    {
        if(zone_specific(deriv) && deriv_zone_index!=0) {
            const DissipatingZone &deriv_zone=zone(deriv_zone_index);
            double above_frac=
                __above_lock_fractions[Dissipation::NO_DERIV]
                                                [deriv_zone.locked_zone_index()];
            Eigen::Vector3d result=zone_to_zone_transform(deriv_zone, zone(0),
                    (above_frac-1.0)
                    *__tidal_torques_below[deriv_zone_index][deriv]
                    -
                    above_frac*__tidal_torques_above[deriv_zone_index][deriv]);
#ifdef DEBUG
            assert(above_lock_fraction_deriv.size()==
                    __above_lock_fractions[Dissipation::NO_DERIV].size());
#endif
            for(unsigned i=0; i<above_lock_fraction_deriv.size(); ++i)
                result+=
                    above_lock_fraction_deriv[i]*__orbit_torque_correction[i];
            if(deriv==Dissipation::INCLINATION ||
                    deriv==Dissipation::PERIAPSIS)
                result+=zone_to_zone_transform(deriv_zone, zone(0),
                        above_frac*__tidal_torques_above[deriv_zone_index]
                                                        [Dissipation::NO_DERIV]
                        +
                        (1.0-above_frac)
                        *__tidal_torques_below[deriv_zone_index]
                                              [Dissipation::NO_DERIV],
                        deriv, true);
            return result;
        } else return __orbit_torque[deriv];
    }

    Eigen::Vector3d DissipatingBody::tidal_orbit_torque(
            const DissipatingZone &reference_zone,
            Dissipation::Derivative deriv, 
            unsigned deriv_zone_index,
            const Eigen::VectorXd &above_lock_fraction_deriv) const
    {
        Eigen::Vector3d result=zone_to_zone_transform(
                zone(0), reference_zone,
                tidal_orbit_torque(deriv, deriv_zone_index,
                                   above_lock_fraction_deriv));
        if((deriv==Dissipation::INCLINATION || deriv==Dissipation::PERIAPSIS)
                &&
                (deriv_zone_index==0 || 
                 &zone(deriv_zone_index)==&reference_zone)) {
            result+=zone_to_zone_transform(zone(0), reference_zone,
                    tidal_orbit_torque(), deriv, deriv_zone_index==0);
        } 
        return result;
    }

    void DissipatingBody::add_to_evolution()
    {
        for(unsigned zone_ind=0; zone_ind<number_zones(); ++zone_ind)
            zone(zone_ind).add_to_evolution();
    }

    void DissipatingBody::rewind_evolution(unsigned nsteps)
    {
        for(unsigned zone_ind=0; zone_ind<number_zones(); ++zone_ind)
            zone(zone_ind).rewind_evolution(nsteps);
    }

    void DissipatingBody::reset_evolution()
    {
        for(unsigned zone_ind=0; zone_ind<number_zones(); ++zone_ind)
            zone(zone_ind).reset_evolution();
    }

    CombinedStoppingCondition *DissipatingBody::stopping_conditions(
            BinarySystem &system, bool primary)
    {
        CombinedStoppingCondition *result=new CombinedStoppingCondition();
        for(unsigned zone_ind=0; zone_ind<number_zones(); ++zone_ind)
            (*result)|=zone(zone_ind).stopping_conditions(system, primary,
                                                          zone_ind);
        return result;
    }

} //End Envolve namespace.
