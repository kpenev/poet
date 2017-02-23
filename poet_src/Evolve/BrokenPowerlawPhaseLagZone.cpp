#include "BrokenPowerlawPhaseLagZone.h"
#include <iostream>

namespace Evolve {

    void BrokenPowerlawPhaseLagZone::fill_tidal_frequency_conditions(
        BinarySystem &system, 
        bool primary,
        unsigned zone_index
    )
    {
        if(eccentricity_order() < __tidal_frequency_conditions.size()) {
            for(
                unsigned e_order = __tidal_frequency_conditions.size() - 1;
                e_order > eccentricity_order();
                --e_order
            ) {
                delete __tidal_frequency_conditions.back();
                __tidal_frequency_conditions.pop_back();
            }
            return;
        }

        const DissipatingBody 
            &this_body = (primary ? system.primary() : system.secondary()),
            &other_body = (primary ? system.secondary() : system.primary());
        double orbital_frequency = Core::orbital_angular_velocity(
            this_body.mass(),
            other_body.mass(),
            system.semimajor()
        );

        for(
            int e_order = __tidal_frequency_conditions.size(); 
            e_order <= static_cast<int>(eccentricity_order());
            ++e_order
        )
        {
            CombinedStoppingCondition *condition = 
                new CombinedStoppingCondition();
            int mp_step = (e_order == 0 ? 1 : 4 + 2 * e_order);

            for(int mp = -e_order - 2; mp <= e_order + 2;  mp += mp_step)
                for(int m = -2; m <= 2; ++m)
                    (*condition) |= new CriticalForcingFrequencyCondition(
                        this_body,
                        other_body,
                        primary,
                        zone_index,
                        mp,
                        m,
                        __tidal_frequency_breaks,
                        orbital_frequency
                    );
            __tidal_frequency_conditions.push_back(condition);
        }

    }

    void BrokenPowerlawPhaseLagZone::print_configuration()
    {
        std::clog << "Tidal breaks: ";
        for(
            std::vector<double>::const_iterator
                i = __tidal_frequency_breaks.begin();
            i != __tidal_frequency_breaks.end();
            ++i
        )
            std::clog << *i << " ";
        std::clog << std::endl;

        std::clog << "Tidal powers: ";
        for(
            std::vector<double>::const_iterator
                i = __tidal_frequency_powers.begin();
            i != __tidal_frequency_powers.end();
            ++i
        )
            std::clog << *i << " ";
        std::clog << std::endl;

        std::clog << "Spin breaks: ";
        for(
            std::vector<double>::const_iterator
                i = __spin_frequency_breaks.begin();
            i != __spin_frequency_breaks.end();
            ++i
        )
            std::clog << *i << " ";
        std::clog << std::endl;

        std::clog << "Spin powers: ";
        for(
            std::vector<double>::const_iterator
                i = __spin_frequency_powers.begin();
            i != __spin_frequency_powers.end();
            ++i
        )
            std::clog << *i << " ";
        std::clog << std::endl;

        std::clog << "Break lags: ";
        for(
            std::vector<double>::const_iterator
                i = __break_phase_lags.begin();
            i != __break_phase_lags.end();
            ++i
        )
            std::clog << *i << " ";
        std::clog << std::endl;
    }

    void BrokenPowerlawPhaseLagZone::setup(
        std::vector<double> tidal_frequency_breaks,
        std::vector<double> spin_frequency_breaks,
        std::vector<double> tidal_frequency_powers,
        std::vector<double> spin_frequency_powers,
        double reference_phase_lag
    )
    {
        assert(__spin_frequency_breaks.size() == 0);
        assert(__tidal_frequency_breaks.size() == 0);
        assert(__spin_frequency_breaks.size() == 0);
        assert(__spin_frequency_powers.size() == 0);
        assert(__break_phase_lags.size() == 0);
        assert(tidal_frequency_powers.size()
               == 
               tidal_frequency_breaks.size() + 1);
        assert(spin_frequency_powers.size()
               ==
               spin_frequency_breaks.size() + 1);
        assert(tidal_frequency_breaks.size() > 0
               ||
               tidal_frequency_powers.front() == 0);
        assert(spin_frequency_breaks.size() > 0
               ||
               spin_frequency_powers.front() == 0);

        __tidal_frequency_breaks = tidal_frequency_breaks;
        __spin_frequency_breaks = spin_frequency_breaks;
        __tidal_frequency_powers = tidal_frequency_powers;
        __spin_frequency_powers = spin_frequency_powers;
        __break_phase_lags.resize(
            std::max(int(spin_frequency_breaks.size()), 1)
            *
            std::max(int(tidal_frequency_breaks.size()), 1)
        );

        unsigned break_lag_i = 0;
        for(
            int spin_break_i = 0;
            spin_break_i < std::max(int(spin_frequency_breaks.size()), 1);
            ++spin_break_i
        ) {
            if(break_lag_i == 0)
                __break_phase_lags[break_lag_i] = reference_phase_lag;
            else 
                __break_phase_lags[break_lag_i] = (
                    __break_phase_lags[
                        break_lag_i - tidal_frequency_breaks.size()
                    ]
                    *
                    (
                        spin_frequency_powers[spin_break_i] == 0
                        ? 1.0
                        : std::pow(
                            (
                                spin_frequency_breaks[spin_break_i]
                                /
                                spin_frequency_breaks[spin_break_i - 1]
                            ),
                            spin_frequency_powers[spin_break_i]
                        )
                    )
                );
            ++break_lag_i;
            for(
                int tidal_break_i = 1;
                tidal_break_i < std::max(int(tidal_frequency_breaks.size()),
                                         1);
                ++tidal_break_i
            ) {
                __break_phase_lags[break_lag_i] = (
                    __break_phase_lags[break_lag_i - 1]
                    *
                    (
                        tidal_frequency_powers[spin_break_i] == 0
                        ? 1.0
                        : std::pow(
                            (
                                tidal_frequency_breaks[tidal_break_i]
                                /
                                tidal_frequency_breaks[tidal_break_i - 1]
                            ),
                            tidal_frequency_powers[spin_break_i]
                        )
                    )
                );
                ++break_lag_i;
            }
        }

#ifndef NDEBUG
        print_configuration();
#endif
    }//End BrokenPowerlawPhaseLagZone::BrokenPowerlawPhaseLagZone definition.

    double BrokenPowerlawPhaseLagZone::modified_phase_lag(
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        double forcing_frequency,
        Dissipation::Derivative deriv,
        double &above_lock_value
    ) const
    {
        if(deriv == Dissipation::AGE) return 0;

        double abs_forcing_frequency = std::abs(forcing_frequency),
               abs_spin_frequency = std::abs(spin_frequency());

        std::vector<double>::const_iterator 
            tidal_break_iter = std::lower_bound(
                __tidal_frequency_breaks.begin(),
                __tidal_frequency_breaks.end(),
                abs_forcing_frequency
            ),
            spin_break_iter = std::lower_bound(
                __spin_frequency_breaks.begin(),
                __spin_frequency_breaks.end(),
                abs_spin_frequency
            );
        if(tidal_break_iter == __tidal_frequency_breaks.end())
            --tidal_break_iter;
        if(spin_break_iter == __spin_frequency_breaks.end())
            --spin_break_iter;
        unsigned tidal_index = (tidal_break_iter
                                -
                                __tidal_frequency_breaks.begin()),
                 spin_index = (spin_break_iter
                               -
                               __spin_frequency_breaks.begin()),
                 lag_index = (spin_index * __tidal_frequency_breaks.size()
                              +
                              tidal_index);
        double 
            tidal_power = __tidal_frequency_powers[tidal_index],
            spin_power = __spin_frequency_powers[spin_index],
            result = (
                __break_phase_lags[lag_index]
                *
                (
                    tidal_power == 0
                    ? 1.0
                    : std::pow(abs_forcing_frequency / *tidal_break_iter,
                               tidal_power)
                )
                *
                (
                    spin_power == 0
                    ? 1.0
                    : std::pow(abs_spin_frequency / *spin_break_iter,
                               spin_power)
                )
            );
        std::clog << "lag index: " << lag_index << std::endl;
        std::clog << "base lag: " 
                  << __break_phase_lags[lag_index]
                  << std::endl;
        std::clog << "scaled lag: " << result << std::endl;
        switch(deriv) {
            case Dissipation::SPIN_FREQUENCY :
                result *= (
                    spin_power / spin_frequency()
                    -
                    spin_frequency_multiplier
                    *
                    tidal_power
                    /
                    forcing_frequency
                );
                break;
            case Dissipation::ORBITAL_FREQUENCY :
                result *= (
                    orbital_frequency_multiplier
                    * 
                    tidal_power
                    /
                    forcing_frequency
                );
                break;
            default :
                assert(deriv == Dissipation::NO_DERIV);
        }

        if(forcing_frequency == 0) above_lock_value  = -result;
        return (forcing_frequency >= 0 ? result : -result);

    }//End BrokenPowerlawPhaseLagZone::modified_phase_lag definition.

    CombinedStoppingCondition *
        BrokenPowerlawPhaseLagZone::stopping_conditions(
            BinarySystem &system, 
            bool primary,
            unsigned zone_index
        )
    {
        CombinedStoppingCondition *result = 
            DissipatingZone::stopping_conditions(system,
                                                 primary,
                                                 zone_index);
        if(system.evolution_mode() != Core::BINARY) return result;

        fill_tidal_frequency_conditions(system, primary, zone_index);

        if(__spin_condition == NULL)
            __spin_condition = new CriticalSpinCondition(
                (primary ? system.primary() : system.secondary()),
                (primary ? system.secondary() : system.primary()),
                primary,
                zone_index,
                __spin_frequency_breaks
            );

        (*result) |= __spin_condition;
        for(
            std::list<CombinedStoppingCondition *>::const_iterator 
                tidal_cond_iter = __tidal_frequency_conditions.begin();
            tidal_cond_iter != __tidal_frequency_conditions.end();
            ++tidal_cond_iter
        )
            (*result) |= *tidal_cond_iter;

        return result;
    }//End BrokenPowerlawPhaseLagZone::stopping_conditions definition.

    void BrokenPowerlawPhaseLagZone::change_e_order(
        unsigned new_e_order,
        BinarySystem &system, 
        bool primary,
        unsigned zone_index
    )
    {
        assert(eccentricity_order() + 1 == __tidal_frequency_conditions.size());

        DissipatingZone::change_e_order(new_e_order,
                                        system,
                                        primary,
                                        zone_index);
        fill_tidal_frequency_conditions(system, primary, zone_index);
    }

    BrokenPowerlawPhaseLagZone::~BrokenPowerlawPhaseLagZone()
    {
        delete __spin_condition;
        for(
            std::list<CombinedStoppingCondition *>::iterator
                condition = __tidal_frequency_conditions.begin();
            condition != __tidal_frequency_conditions.end();
            ++condition
        )
            delete *condition;
    }

} //End BinarySystem namespace.
