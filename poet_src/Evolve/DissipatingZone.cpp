#define BUILDING_LIBRARY
#include "DissipatingZone.h"
#include "BinarySystem.h"

namespace Evolve {

    std::ostream &operator<<(std::ostream &os,
                             const ZoneEvolutionQuantities &evol_var)
    {
        switch(evol_var) {
            case ANGULAR_MOMENTUM : os << "AGE"; break;
            case INCLINATION : os << "INCLINATION"; break;
            case PERIAPSIS : os << "PERIAPSIS"; break;
            case MOMENT_OF_INERTIA : os << "MOMENT_OF_INERTIA"; break;
            case MOMENT_OF_INERTIA_FIRST_DERIV :
                                     os << "DMOMENT_OF_INERTIA_DT"; break;
            case MOMENT_OF_INERTIA_SECOND_DERIV :
                                     os << "D2MOMENT_OF_INERTIA_DT2"; break;
            case OUTER_RADIUS : os << "OUTER_RADIUS"; break;
            case OUTER_RADIUS_FIRST_DERIV : os << "DOUTER_RADIUS_DT"; break;
            case OUTER_RADIUS_SECOND_DERIV : os << "D2OUTER_RADIUS_DT2"; break;
            case OUTER_MASS : os << "OUTER_MASS"; break;
            case OUTER_MASS_DERIV : os << "DOUTER_MASS_DT"; break;
            case E_ORDER : os << "E_ORDER"; break;
            case ORBITAL_FREQ_MULTIPLIER : os << "ORBITAL_FREQ_MULTIPLIER";
                                           break;
            case SPIN_FREQ_MULTIPLIER : os << "SPIN_FREQ_MULTIPLIER"; break;
            default : assert(false);
        };
        return os;
    }

    const double DissipatingZone::__torque_x_plus_coef[]={1.0, 
                                                          std::sqrt(1.5),
                                                          std::sqrt(1.5),
                                                          1.0,
                                                          0.0};

    const double DissipatingZone::__torque_x_minus_coef[]={0.0,           //m=-2
                                                           1.0,           //m=-1
                                                           std::sqrt(1.5),//m=0
                                                           std::sqrt(1.5),//m=1
                                                           1.0};          //m=2

    void DissipatingZone::fix_forcing_frequency(
        const SpinOrbitLockInfo &limit, 
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        double &forcing_frequency
    ) const
    {
        if(__initializing) return;
        assert(limit.spin_frequency_multiplier()==1
               ||
               limit.spin_frequency_multiplier()==2);

        if(orbital_frequency_multiplier > 5)
            std::cerr << "ofm" << std::endl;
        if(spin_frequency_multiplier > 5)
            std::cerr << "sfm" << std::endl;

        if(
            limit.term(orbital_frequency_multiplier,
                       spin_frequency_multiplier)
        ) {
            if(
                spin_frequency_multiplier
                *
                limit.lock_direction()
                *
                forcing_frequency
                >
                0
            ) 
                forcing_frequency = (
                    (
                        spin_frequency_multiplier*limit.lock_direction() > 0
                        ? -1
                        : 1
                    )
                    *
                    std::numeric_limits<double>::epsilon()
                );
            return;
        }
        int expected_sign = (
            limit.spin_frequency_multiplier()
            *
            (
                orbital_frequency_multiplier
                *
                limit.spin_frequency_multiplier()
                -
                limit.orbital_frequency_multiplier()
                *
                spin_frequency_multiplier
            )
        );
        if(expected_sign * limit.lock_direction() > 0) return;
        if(forcing_frequency * expected_sign > 0) return;

        assert(limit.lock_direction());

        forcing_frequency = (std::numeric_limits<double>::epsilon()
                             *
                             limit.lock_direction());
    }

#ifndef NDEBUG
    void DissipatingZone::check_locks_consistency() const
    {
        if(__initializing) return;
        int max_abs_orb_mult = static_cast<int>(__e_order + 2);
        assert(
            __lock
            || 
            (__lock.lock_direction() * __other_lock.lock_direction( )== -1)
            ||
            (
                __lock.spin_frequency_multiplier() == 1
                && 
                (
                    __lock.orbital_frequency_multiplier()
                    ==
                    __lock.lock_direction() * max_abs_orb_mult
                )
                &&
                __other_lock.orbital_frequency_multiplier() == 1
                && 
                __other_lock.spin_frequency_multiplier() == 0
                &&
                __other_lock.lock_direction() == 1)
        );
        assert(__lock.spin_frequency_multiplier() == 1
               ||
               __lock.spin_frequency_multiplier() == 2);
        assert(__other_lock.spin_frequency_multiplier() >= 0
               && 
               __other_lock.spin_frequency_multiplier() <= 2);
        if(__lock) return;
        return;//<++>
        assert(
            (
                __lock.lock_direction()
                *
                __lock.spin_frequency_multiplier()
                *
                spin_frequency()
                +
                1.0e-5 * __orbital_frequency
            )
            >=
            (
                __lock.lock_direction()
                *
                __lock.orbital_frequency_multiplier()
                *
                __orbital_frequency
            )
        );
        if(__other_lock.spin_frequency_multiplier()) {
            assert(
                (
                    __other_lock.lock_direction()
                    *
                    __other_lock.spin_frequency_multiplier()
                    *
                    spin_frequency()
                )
                >
                (
                    __other_lock.lock_direction()
                    *
                    __other_lock.orbital_frequency_multiplier()
                    *
                    __orbital_frequency
                )
            );
        } else assert(
            __lock.lock_direction() * __lock.orbital_frequency_multiplier()
            >
            0
        );
    }
#endif

    void DissipatingZone::update_lock_to_lower_e_order(SpinOrbitLockInfo &lock)
    {
#ifndef NDEBUG
        assert(lock.lock_direction());
        check_locks_consistency();
#endif
        if(static_cast<unsigned>(std::abs(lock.orbital_frequency_multiplier()))
                >__e_order+2
                && 
                lock.spin_frequency_multiplier()==2)
            lock.set_lock((lock.orbital_frequency_multiplier()
                           -
                           lock.lock_direction())/2,
                          1, lock.lock_direction());
        if(lock.orbital_frequency_multiplier()>static_cast<int>(__e_order+2)) {
            if(lock.lock_direction()>0)
                lock.set_lock(__e_order+2, 1, 1);
            else lock.set_lock(1, 0, 1);
        } else if(lock.orbital_frequency_multiplier()<
                -static_cast<int>(__e_order)-2) {
            if(lock.lock_direction()>0) lock.set_lock(1, 0, 1);
            else lock.set_lock(-static_cast<int>(__e_order)-2, 1, -1);
        }
#ifndef NDEBUG
        check_locks_consistency();
#endif
    }

    void DissipatingZone::update_locks_to_higher_e_order(unsigned new_e_order)
    {
        if(__other_lock.spin_frequency_multiplier()==0) {
            int orb_mult=std::ceil(std::abs(2*spin_frequency())
                    /__orbital_frequency);
            if(orb_mult%2==1 && static_cast<int>(new_e_order)+2>=orb_mult)
                __other_lock.set_lock(__lock.lock_direction()*orb_mult, 2,
                        -__lock.lock_direction());
            else if(static_cast<int>(new_e_order)+2>=orb_mult/2)
                __other_lock.set_lock(__lock.lock_direction()*orb_mult/2, 1,
                        -__lock.lock_direction());
        } else if(__lock.spin_frequency_multiplier()==1
                && __other_lock.spin_frequency_multiplier()==1
                && __lock.orbital_frequency_multiplier()
                +
                __other_lock.orbital_frequency_multiplier()
                <=
                static_cast<int>(new_e_order)+2) {
            int mid_mult=__lock.orbital_frequency_multiplier()
                +
                __other_lock.orbital_frequency_multiplier();
            if(2*spin_frequency()<mid_mult*__orbital_frequency)
                (__lock.lock_direction()>0 ? __other_lock : __lock)
                    .set_lock(mid_mult, 2, -1);
            else (__lock.lock_direction()>0 ? __lock : __other_lock)
                .set_lock(mid_mult, 2, 1);
        }
    }

    void DissipatingZone::initialize_locks()
    {
#ifndef NDEBUG
        std::cerr << "Initializing locks" << std::endl;
#endif
        int below_orb_mult=std::floor(2.0*__spin_frequency/__orbital_frequency),
            max_abs_orb_mult=static_cast<int>(__e_order+2);
        if(below_orb_mult%2) {
            if(static_cast<unsigned>(std::abs(below_orb_mult))<=__e_order+2) {
                __lock.set_lock(below_orb_mult, 2, 1);
                __other_lock.set_lock((below_orb_mult+1)/2, 1, -1);
            } else if(static_cast<unsigned>(std::abs((below_orb_mult-1)/2))
                      <=__e_order+2) {
                __lock.set_lock((below_orb_mult-1)/2, 1, 1);
                if((below_orb_mult+1)/2>max_abs_orb_mult)
                    __other_lock.set_lock(1, 0, 1);
                else __other_lock.set_lock((below_orb_mult+1)/2, 1, -1);
            } else {
                if(__spin_frequency>0) __lock.set_lock(max_abs_orb_mult, 1, 1);
                else __lock.set_lock(-max_abs_orb_mult, 1, -1);
                __other_lock.set_lock(1, 0, 1);
            }
        } else if(std::abs(below_orb_mult/2)<=max_abs_orb_mult) {
            __lock.set_lock(below_orb_mult/2, 1, 1);
            if(std::abs(below_orb_mult+1)<=max_abs_orb_mult)
                __other_lock.set_lock(below_orb_mult+1, 2, -1);
            else if(std::abs(below_orb_mult/2+1)<=max_abs_orb_mult)
                __other_lock.set_lock(below_orb_mult/2+1, 1, -1);
            else __other_lock.set_lock(1, 0, 1);
        } else {
            if(__spin_frequency>0) __lock.set_lock(max_abs_orb_mult, 1, 1);
            else __lock.set_lock(-max_abs_orb_mult, 1, -1);
            __other_lock.set_lock(1, 0, 1);
        }
#ifndef NDEBUG
        check_locks_consistency();
#endif
    }

    void DissipatingZone::add_tidal_term(int m,
                                         int mp,
                                         double tidal_frequency,
                                         const TidalTermTriplet &U_value,
                                         const TidalTermTriplet &U_i_deriv,
                                         const TidalTermTriplet &U_e_deriv,
                                         const TidalTermTriplet &U_error)
    {
        int m_ind = m + 2;

        bool locked_term = locked(mp, m);

        for(
            int deriv = Dissipation::NO_DERIV; 
            (
                (m != 0 || mp != 0)
                && 
                deriv < Dissipation::END_DIMENSIONLESS_DERIV
            );
            ++deriv
        ) {
            Dissipation::Derivative phase_lag_deriv =
                (deriv < Dissipation::END_PHASE_LAG_DERIV
                 ? static_cast<Dissipation::Derivative>(deriv)
                 : Dissipation::NO_DERIV);
            double mod_phase_lag_above,
                   mod_phase_lag_below = modified_phase_lag(
                       mp,
                       m,
                       tidal_frequency,
                       phase_lag_deriv,
                       mod_phase_lag_above
                   ),
                   love_coef = love_coefficient(
                       mp,
                       m,
                       (
                           phase_lag_deriv == Dissipation::AGE
                           ? Dissipation::AGE
                           : Dissipation::NO_DERIV
                       )
                   );

            TidalTermTriplet U;
            if(deriv < Dissipation::END_PHASE_LAG_DERIV)
                U = U_value;
            else if(deriv == Dissipation::INCLINATION)
                U = U_i_deriv;
            else
                U = U_e_deriv;

            double U_mmp_squared = std::pow(U.m, 2),
                   term_power = U_mmp_squared*mp,
                   term_torque_z = U_mmp_squared*m,
                   term_torque_x = U.m * (
                       __torque_x_minus_coef[m_ind] * U.m_minus_one
                       +
                       __torque_x_plus_coef[m_ind] * U.m_plus_one
                   );
            if(
                !locked_term
                &&
                (tidal_frequency != 0 || __lock.lock_direction() < 0)
            )
                mod_phase_lag_above = mod_phase_lag_below;
            else if(
                !locked_term
                &&
                tidal_frequency == 0
                &&
                __lock.lock_direction() > 0
            )
                mod_phase_lag_below = mod_phase_lag_above;
            int deriv_ind = 2 * deriv;
            __power[deriv_ind] += term_power * mod_phase_lag_below;
            __torque_z[deriv_ind] += (term_torque_z
                                      *
                                      mod_phase_lag_below);
            __torque_x[deriv_ind] += (term_torque_x
                                      *
                                      mod_phase_lag_below);
            __torque_y[deriv_ind + 1] = -(
                __torque_y[deriv_ind] -= term_torque_x*love_coef
            );
            __power[deriv_ind + 1] += term_power * mod_phase_lag_above;
            __torque_z[deriv_ind + 1] += (term_torque_z
                                          *
                                          mod_phase_lag_above);
            __torque_x[deriv_ind + 1] += (term_torque_x
                                          *
                                          mod_phase_lag_above);
#ifdef VERBOSE_DEBUG
            if(deriv == Dissipation::NO_DERIV)
                std::cerr << ", Wzone = "
                    << spin_frequency()
                    << ", U(" << m << ", " << mp << ") = "
                    << U_mmp
                    << ", term_power="
                    << term_power 
                    << ", mod_phase_lag(above="
                    << mod_phase_lag_above
                    << ", below="
                    << mod_phase_lag_below
                    << ")";
#endif
        }

    }

    void DissipatingZone::configure_spin(double spin,
                                         bool spin_is_frequency)
    {
        if(spin_is_frequency) {
            __angular_momentum = spin * moment_of_inertia(); 
            __spin_frequency = spin;
        } else {
            __angular_momentum = spin; 
            if(spin == 0 && moment_of_inertia() == 0) __spin_frequency = 0;
            else __spin_frequency = spin / moment_of_inertia();
        }
    }

    DissipatingZone::DissipatingZone() :
        __e_order(0),
        __power(0.0, 2*Dissipation::END_DIMENSIONLESS_DERIV),
        __torque_x(0.0, 2*Dissipation::END_DIMENSIONLESS_DERIV),
        __torque_y(0.0, 2*Dissipation::END_DIMENSIONLESS_DERIV),
        __torque_z(0.0, 2*Dissipation::END_DIMENSIONLESS_DERIV),
        __evolution_real(NUM_REAL_EVOL_QUANTITIES),
        __evolution_integer(NUM_EVOL_QUANTITIES - NUM_REAL_EVOL_QUANTITIES),
        __initializing(false)
    {}

    void DissipatingZone::configure(bool initialize,
                                    double
#ifndef NDEBUG
                                    age
#endif
                                    ,
                                    double orbital_frequency,
                                    double eccentricity,
                                    double orbital_angmom,
                                    double spin,
                                    double inclination,
                                    double periapsis,
                                    bool spin_is_frequency)
    {
        assert(age >= 0);

        if(initialize) {
            __initializing = true;
#ifndef NDEBUG
            std::cerr << "Initializing DissipatingZone" << std::endl;
#endif
        }
        ZoneOrientation::configure(inclination, periapsis);
        __orbital_angmom = orbital_angmom;
        __orbital_frequency = orbital_frequency;
        if(__lock && !initialize) {
            __spin_frequency = __lock.spin(orbital_frequency);
            __angular_momentum = __spin_frequency*moment_of_inertia();
        } else 
            configure_spin(spin, spin_is_frequency);

        if(initialize) {
            initialize_locks();
            __initializing = false;
        }
        if(std::isnan(orbital_frequency)) return;

        __potential_term.configure(inclination);
        __power = 0;
        __torque_x = 0;
        __torque_y = 0;
        __torque_z = 0;

        if(!dissipative()) return;

        for(
            int mp = -static_cast<int>(__e_order) - 2;
            mp <= static_cast<int>(__e_order) + 2; 
            ++mp
        ) {
            TidalTermTriplet U_value,
                             U_error,
                             U_i_deriv,
                             U_e_deriv;
            __potential_term(eccentricity,
                             -2,
                             mp,
                             U_value.m_plus_one,
                             U_i_deriv.m_plus_one,
                             U_e_deriv.m_plus_one,
                             U_error.m_plus_one);

            for(int m = -2; m <= 2; ++m) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Term: m' = "
                          << mp
                          << ", m = "
                          << m;
#endif

                U_value.m = U_value.m_plus_one;
                U_error.m = U_error.m_plus_one;
                U_i_deriv.m = U_i_deriv.m_plus_one;
                U_e_deriv.m = U_e_deriv.m_plus_one;
                if(m < 2)
                    __potential_term(eccentricity,
                                     m + 1,
                                     mp,
                                     U_value.m_plus_one,
                                     U_i_deriv.m_plus_one,
                                     U_e_deriv.m_plus_one,
                                     U_error.m_plus_one);
                else {
                    U_value.m_plus_one = 0;
                    U_error.m_plus_one = 0;
                    U_i_deriv.m_plus_one = 0;
                    U_e_deriv.m_plus_one = 0;
                }

                add_tidal_term(m,
                               mp, 
                               forcing_frequency(mp,
                                                 m,
                                                 orbital_frequency),
                               U_value,
                               U_i_deriv,
                               U_e_deriv,
                               U_error);

                U_value.m_minus_one = U_value.m;
                U_i_deriv.m_minus_one = U_i_deriv.m;
                U_e_deriv.m_minus_one = U_e_deriv.m;
#ifdef VERBOSE_DEBUG
                std::cerr << std::endl;
#endif
            }
        }
    }

    double DissipatingZone::forcing_frequency(
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        double orbital_frequency
    ) const
    {
#ifndef NDEBUG
        check_locks_consistency();
#endif
        if(__lock(orbital_frequency_multiplier, spin_frequency_multiplier))
            return 0;
        double forcing_freq = (
            orbital_frequency_multiplier * orbital_frequency
            -
            spin_frequency_multiplier * spin_frequency()
        );

#ifdef VERBOSE_DEBUG
        std::cerr << ", Wtide = " << forcing_freq << " -> ";
#endif

        fix_forcing_frequency(__lock,
                              orbital_frequency_multiplier,
                              spin_frequency_multiplier,
                              forcing_freq);
        if(!__lock && __other_lock.spin_frequency_multiplier()!=0)
            fix_forcing_frequency(__other_lock, orbital_frequency_multiplier,
                    spin_frequency_multiplier, forcing_freq);
#ifdef VERBOSE_DEBUG
        std::cerr << forcing_freq;
#endif
        return forcing_freq;
    }



    double DissipatingZone::periapsis_evolution(
            const Eigen::Vector3d &orbit_torque,
            const Eigen::Vector3d &zone_torque,
            Dissipation::Derivative deriv,
            const Eigen::Vector3d &orbit_torque_deriv,
            const Eigen::Vector3d &zone_torque_deriv)
    {
        double sin_inc = std::sin(inclination()),
               cos_inc = std::cos(inclination()),
               zone_y_torque,
               orbit_y_torque;
        if(deriv == Dissipation::NO_DERIV) {
            orbit_y_torque = orbit_torque[1];
            zone_y_torque = zone_torque[1];
        } else {
            orbit_y_torque = orbit_torque_deriv[1];
            zone_y_torque = zone_torque_deriv[1];
        }
#ifndef NDEBUG
        if(sin_inc == 0) {
            assert(orbit_y_torque == 0 || std::isnan(orbit_y_torque));
            assert(zone_y_torque == 0 || std::isnan(zone_y_torque));
        } else {
            assert(!std::isnan(orbit_y_torque));
            assert(!std::isnan(zone_y_torque));
        }
#endif
        double result = (
            sin_inc == 0
            ? 0
            : (
                -orbit_y_torque * cos_inc / (__orbital_angmom * sin_inc)
                +
                zone_y_torque / (angular_momentum() * sin_inc)
            )
        );
        assert(!std::isnan(result));

        if(
            deriv == Dissipation::NO_DERIV 
            || deriv == Dissipation::AGE 
            || deriv == Dissipation::ECCENTRICITY
            || deriv == Dissipation::PERIAPSIS
            || deriv == Dissipation::RADIUS
            || deriv == Dissipation::MOMENT_OF_INERTIA
            || deriv == Dissipation::SEMIMAJOR
        )
            return result;
        else if(
            deriv == Dissipation::SPIN_FREQUENCY ||
            deriv == Dissipation::SPIN_ANGMOM
        ) {
            if(sin_inc == 0) return 0.0;
            else return (
                result
                -
                zone_torque[1]
                /
                (std::pow(angular_momentum(), 2) * sin_inc)
                *
                (deriv == Dissipation::SPIN_FREQUENCY ? moment_of_inertia() : 1)
            );
        } else if(deriv == Dissipation::INCLINATION) {
            if(sin_inc == 0) return 0.0;
            else return (
                result
                -
                (
                    orbit_torque[1] / __orbital_angmom
                    +
                    zone_torque[1] * cos_inc / angular_momentum()
                )
                /
                std::pow(sin_inc, 2)
            );
        } else {
            assert(false);
        }

        return Core::NaN;
    }

    double DissipatingZone::inclination_evolution(
            const Eigen::Vector3d &orbit_torque,
            const Eigen::Vector3d &zone_torque,
            Dissipation::Derivative deriv,
            const Eigen::Vector3d &orbit_torque_deriv,
            const Eigen::Vector3d &zone_torque_deriv)
    {
        double sin_inc=std::sin(inclination()),
               cos_inc=std::cos(inclination()),
               zone_x_torque,
               orbit_x_torque,
               orbit_z_torque;
        if(deriv==Dissipation::NO_DERIV) {
            orbit_x_torque=orbit_torque[0];
            orbit_z_torque=orbit_torque[2];
            zone_x_torque=zone_torque[0];
        } else {
            orbit_x_torque=orbit_torque_deriv[0];
            orbit_z_torque=orbit_torque_deriv[2];
            zone_x_torque=zone_torque_deriv[0];
        }
        double result;
        if(orbit_x_torque == 0 && orbit_z_torque == 0)
            result = 0.0;
        else 
            result = ((orbit_x_torque * cos_inc - orbit_z_torque*sin_inc)
                      /
                      __orbital_angmom);

        if(zone_x_torque!=0 && moment_of_inertia()!=0)
            result -= zone_x_torque / __angular_momentum;
        if(		deriv==Dissipation::NO_DERIV 
                || deriv==Dissipation::AGE 
                || deriv==Dissipation::ECCENTRICITY
                || deriv==Dissipation::PERIAPSIS
                || deriv==Dissipation::RADIUS
                || deriv==Dissipation::MOMENT_OF_INERTIA
                || deriv==Dissipation::SEMIMAJOR)
            return result;
        else if(deriv==Dissipation::SPIN_FREQUENCY ||
                deriv==Dissipation::SPIN_ANGMOM)
            return result
                   +
                   zone_torque[0]/std::pow(angular_momentum(), 2)
                   *(deriv==Dissipation::SPIN_FREQUENCY ? moment_of_inertia():1);
        else if(deriv==Dissipation::INCLINATION)
            return result
                   +
                   (orbit_torque[2]*cos_inc + orbit_torque[0]*sin_inc)
                   /angular_momentum();
        else assert(false);
        return Core::NaN;
    }

    void DissipatingZone::release_lock()
    {
        if(__lock.spin_frequency_multiplier()==2) {
            assert(__lock.orbital_frequency_multiplier()%2==1);
            __other_lock.set_lock((__lock.orbital_frequency_multiplier()+1)/2,
                    1, -1);
            __lock.set_lock((__lock.orbital_frequency_multiplier()-1)/2, 1, -1);
        }
        update_lock_to_lower_e_order(__lock);
        update_lock_to_lower_e_order(__other_lock);
        if(__lock.spin_frequency_multiplier()==0) {
            __lock=__other_lock;
            __other_lock.set_lock(1, 0, 1);
        }
    }

    void DissipatingZone::release_lock(short direction)
    {
        assert(__lock);
        assert(direction == 1 || direction == -1);
        __lock.lock_direction(direction);
        int orbit_mult = (
            (__lock.spin_frequency_multiplier() == 2 ? 1 : 2)
            *
            __lock.orbital_frequency_multiplier()
            +
            direction
        );
        if(orbit_mult % 2) __other_lock.set_lock(orbit_mult, 2, -direction);
        else __other_lock.set_lock(orbit_mult/2, 1, -direction);
        update_lock_to_lower_e_order(__other_lock);
    }

    void DissipatingZone::change_e_order(
        unsigned new_e_order,
        BinarySystem &, 
        bool ,
        unsigned
    )
    {
        __potential_term.change_e_order(new_e_order);
        if(__lock.spin_frequency_multiplier()==0) {
            __e_order=new_e_order;
            return;
        }
        if(__lock) {
           __e_order=new_e_order;
           if(__lock.orbital_frequency_multiplier()>
                   static_cast<int>(__e_order)+2) release_lock();
           return;
        }
#ifndef NDEBUG
        check_locks_consistency();
#endif
        if(new_e_order>__e_order) {
            update_locks_to_higher_e_order(new_e_order);
            __e_order=new_e_order;
        } else {
            __e_order=new_e_order;
            update_lock_to_lower_e_order(__lock);
            update_lock_to_lower_e_order(__other_lock);
            if(__lock.spin_frequency_multiplier()==0) {
                __lock=__other_lock;
                __other_lock.set_lock(1, 0, 1);
            }
        }
#ifndef NDEBUG
        check_locks_consistency();
#endif
    }

    void DissipatingZone::add_to_evolution()
    {
        __evolution_real[ANGULAR_MOMENTUM].push_back(__angular_momentum);

        __evolution_real[INCLINATION].push_back(inclination());

        __evolution_real[PERIAPSIS].push_back(periapsis());

        __evolution_real[MOMENT_OF_INERTIA].push_back(moment_of_inertia());

        __evolution_real[MOMENT_OF_INERTIA_FIRST_DERIV].push_back(
                moment_of_inertia(1)
        );

        __evolution_real[MOMENT_OF_INERTIA_SECOND_DERIV].push_back(
                moment_of_inertia(2)
        );

        __evolution_real[OUTER_RADIUS].push_back(outer_radius());

        __evolution_real[OUTER_RADIUS_FIRST_DERIV].push_back(outer_radius(1));

        __evolution_real[OUTER_RADIUS_SECOND_DERIV].push_back(outer_radius(2));

        __evolution_real[OUTER_MASS].push_back(outer_mass());

        __evolution_real[OUTER_MASS_DERIV].push_back(outer_mass(1));

        __evolution_integer[E_ORDER-NUM_REAL_EVOL_QUANTITIES].push_back(
                __e_order
        );

        if(__lock) {
            __evolution_integer[ORBITAL_FREQ_MULTIPLIER-NUM_REAL_EVOL_QUANTITIES]
                .push_back(__lock.orbital_frequency_multiplier());
            __evolution_integer[SPIN_FREQ_MULTIPLIER-NUM_REAL_EVOL_QUANTITIES]
                .push_back(__lock.spin_frequency_multiplier());
        } else {
            __evolution_integer[ORBITAL_FREQ_MULTIPLIER-NUM_REAL_EVOL_QUANTITIES]
                .push_back(0);
            __evolution_integer[SPIN_FREQ_MULTIPLIER-NUM_REAL_EVOL_QUANTITIES]
                .push_back(0);
        }
    }

    void DissipatingZone::reset_evolution()
    {
        for(unsigned i=0; i<__evolution_real.size(); ++i)
            __evolution_real[i].clear();
        for(unsigned i=0; i<__evolution_integer.size(); ++i)
            __evolution_integer[i].clear();
    }

    void DissipatingZone::rewind_evolution(unsigned nsteps)
    {
        for(unsigned i=0; i<nsteps; ++i) {
            for(unsigned i=0; i<__evolution_real.size(); ++i)
                __evolution_real[i].pop_back();
            for(unsigned i=0; i<__evolution_integer.size(); ++i)
                __evolution_integer[i].pop_back();
        }
    }

    CombinedStoppingCondition *DissipatingZone::stopping_conditions(
        BinarySystem &system, 
        bool primary,
        unsigned zone_index
    )
    {
        CombinedStoppingCondition *result = new CombinedStoppingCondition();
        if(!dissipative()) return result;
        if(__lock)
            (*result) |= new BreakLockCondition(system, __locked_zone_index);
        else if(system.evolution_mode() == Core::BINARY) {
            (*result) |= new SynchronizedCondition(
                __lock.orbital_frequency_multiplier(),
                __lock.spin_frequency_multiplier(),
                __lock.lock_direction(),
                primary,
                zone_index,
                system
            );
            (*result) |= new SynchronizedCondition(
                __other_lock.orbital_frequency_multiplier(),
                __other_lock.spin_frequency_multiplier(),
                __other_lock.lock_direction(),
                primary,
                zone_index,
                system
            );
        }
        return result;
    }

} //End Evolve namespace.
