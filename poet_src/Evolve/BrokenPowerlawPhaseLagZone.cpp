#define BUILDING_LIBRARY
#include "BrokenPowerlawPhaseLagZone.h"
#include <iostream>

namespace Evolve {

    double BrokenPowerlawPhaseLagZone::get_orbital_frequency(
        const BinarySystem &system
    ) const
    {
        return Core::orbital_angular_velocity(
            system.primary().mass(),
            system.secondary().mass(),
            system.semimajor()
        );
    }

    void BrokenPowerlawPhaseLagZone::reset()
    {
        __tidal_frequency_breaks.clear();
        __spin_frequency_breaks.clear();
        __tidal_frequency_powers.clear();
        __spin_frequency_powers.clear();
        __break_phase_lags.clear();
    }

    void BrokenPowerlawPhaseLagZone::set_spin_index()
    {
        if(__spin_frequency_breaks.size() != 0) {
            double abs_spin_frequency = std::abs(spin_frequency());
            __spin_index = (
                std::lower_bound(
                    __spin_frequency_breaks.begin(),
                    __spin_frequency_breaks.end(),
                    abs_spin_frequency
                )
                -
                __spin_frequency_breaks.begin()
            );
            if(
                __spin_index < __spin_frequency_breaks.size()
                &&
                __spin_frequency_breaks[__spin_index] == abs_spin_frequency
            )
                throw Core::Error::BadFunctionArguments(
                    "Starting evolution from exactly a critical spin "
                    "frequency is not currently supported."
                );
        } else
            __spin_index = 0;


    }

    std::vector<double>::size_type 
        BrokenPowerlawPhaseLagZone::get_tidal_index(
            double abs_forcing_frequency
        )
        {
            assert(abs_forcing_frequency >= 0);
            if(__tidal_frequency_breaks.size() != 0)
                return (
                    std::lower_bound(
                        __tidal_frequency_breaks.begin(),
                        __tidal_frequency_breaks.end(),
                        abs_forcing_frequency
                    )
                    -
                    __tidal_frequency_breaks.begin()
                );
            else
                return 0;
        }

    void BrokenPowerlawPhaseLagZone::add_tidal_frequency_conditions(
        BinarySystem &system, 
        bool primary,
        unsigned,
        CombinedStoppingCondition &result
    )
    {
        const DissipatingBody 
            &this_body = (primary ? system.primary() : system.secondary()),
            &other_body = (primary ? system.secondary() : system.primary());

        for(
            int e_order = 0; 
            e_order <= static_cast<int>(eccentricity_order());
            ++e_order
        )
        {

            int mp_step = (e_order == 0 ? 1 : 4 + 2 * e_order);

            for(int mp = 0; mp <= e_order + 2;  mp += mp_step)
                for(int m = -2; m <= 2; ++m)
                    result |= new LagForcingFrequencyBreakCondition(
                        *this,
                        this_body,
                        other_body,
                        mp,
                        m
                    );
        }

    }

    void BrokenPowerlawPhaseLagZone::print_configuration(std::ostream &out_stream)
    {
        return;
        out_stream << "Tidal breaks: ";
        for(
            std::vector<double>::const_iterator
                i = __tidal_frequency_breaks.begin();
            i != __tidal_frequency_breaks.end();
            ++i
        )
            out_stream << *i << " ";
        out_stream << std::endl;

        out_stream << "Tidal powers: ";
        for(
            std::vector<double>::const_iterator
                i = __tidal_frequency_powers.begin();
            i != __tidal_frequency_powers.end();
            ++i
        )
            out_stream << *i << " ";
        out_stream << std::endl;

        out_stream << "Spin breaks: ";
        for(
            std::vector<double>::const_iterator
                i = __spin_frequency_breaks.begin();
            i != __spin_frequency_breaks.end();
            ++i
        )
            out_stream << *i << " ";
        out_stream << std::endl;

        out_stream << "Spin powers: ";
        for(
            std::vector<double>::const_iterator
                i = __spin_frequency_powers.begin();
            i != __spin_frequency_powers.end();
            ++i
        )
            out_stream << *i << " ";
        out_stream << std::endl;

        out_stream << "Break lags: ";
        for(
            std::vector<double>::const_iterator
                i = __break_phase_lags.begin();
            i != __break_phase_lags.end();
            ++i
        )
            out_stream << *i << " ";
        out_stream << std::endl;
    }

    void BrokenPowerlawPhaseLagZone::setup(
        const std::vector<double> &tidal_frequency_breaks,
        const std::vector<double> &spin_frequency_breaks,
        const std::vector<double> &tidal_frequency_powers,
        const std::vector<double> &spin_frequency_powers,
        double reference_phase_lag
    )
    {
        reset();
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
                        tidal_frequency_powers[tidal_break_i] == 0
                        ? 1.0
                        : std::pow(
                            (
                                tidal_frequency_breaks[tidal_break_i]
                                /
                                tidal_frequency_breaks[tidal_break_i - 1]
                            ),
                            tidal_frequency_powers[tidal_break_i]
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

    void BrokenPowerlawPhaseLagZone::configure(
        bool initialize,
        double age,
        double orbital_frequency,
        double eccentricity,
        double orbital_angmom,
        double spin,
        double inclination,
        double periapsis,
        bool spin_is_frequency
    )
    {

        if(initialize && !std::isnan(orbital_frequency)) {
            initializing(true);

#ifndef NDEBUG
            std::cerr << "Initializing broken powerlaw lag zone"
                      << (initialize ? " for the first time." : ".")
                      << std::endl;
#endif

            configure_spin(spin, spin_is_frequency);

            set_spin_index();

            __tidal_indices.resize(5 * (eccentricity_order() + 3));
#ifndef NDEBUG
            std::cerr << "__tidal_indices size = "
                      << __tidal_indices.size()
                      << std::endl;
#endif

            std::vector< std::vector<double>::size_type >::iterator
                destination = __tidal_indices.begin();

            for(
                int mp = 0;
                mp <= static_cast<int>(eccentricity_order()) + 2;
                ++mp
            ) {
                for(int m = -2; m <= 2; ++m) {
                    double abs_forcing_frequency = std::abs(
                        forcing_frequency(mp, m, orbital_frequency)
                    );
                    *destination = get_tidal_index(abs_forcing_frequency);
                    if(
                        *destination < __tidal_frequency_breaks.size()
                        &&
                        (
                            __tidal_frequency_breaks[*destination]
                            ==
                            abs_forcing_frequency
                        )
                    )
                        throw Core::Error::BadFunctionArguments(
                            "Starting evolution from exactly a critical tidal "
                            "forcing frequency is not currently supported."
                        );
                    ++destination;
                }
            }

            initializing(false);
        }

        DissipatingZone::configure(initialize,
                                   age,
                                   orbital_frequency,
                                   eccentricity,
                                   orbital_angmom,
                                   spin,
                                   inclination,
                                   periapsis,
                                   spin_is_frequency);

    }

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

        std::vector<double>::size_type tidal_index
            = __tidal_indices[tidal_term_index(orbital_frequency_multiplier,
                                               spin_frequency_multiplier)];

        double tidal_power = __tidal_frequency_powers[tidal_index],
               spin_power = __spin_frequency_powers[__spin_index];

        std::vector<double>::size_type tidal_break_index = tidal_index,
                                       spin_break_index = __spin_index;
        if(
            spin_break_index > 0
            &&
            spin_break_index >= __spin_frequency_breaks.size()
        )
            --spin_break_index;

        if(
            tidal_break_index > 0
            &&
            tidal_break_index >= __tidal_frequency_breaks.size()
        )
            --tidal_break_index;

        std::vector<double>::size_type 
            lag_index = (spin_break_index * __tidal_frequency_breaks.size()
                         +
                         tidal_break_index);
        double result = (
            __break_phase_lags[lag_index]
            *
            (
                tidal_power == 0
                ? 1.0
                : std::pow(abs_forcing_frequency
                           /
                           __tidal_frequency_breaks[tidal_break_index]
                           ,
                           tidal_power)
            )
            *
            (
                spin_power == 0
                ? 1.0
                : std::pow(
                    (
                        abs_spin_frequency
                        /
                        __spin_frequency_breaks[spin_break_index]
                    ),
                    spin_power
                )
            )
        );
        switch(deriv) {
            case Dissipation::SPIN_FREQUENCY :
                result *= (
                    (spin_power ? spin_power / spin_frequency() : 0.0)
                    -
                    (
                        tidal_power
                        ? (spin_frequency_multiplier * tidal_power
                           /
                           forcing_frequency)
                        : 0.0
                    )
                );
                break;
            case Dissipation::ORBITAL_FREQUENCY :
                result *= (
                    tidal_power 
                    ? (orbital_frequency_multiplier * tidal_power
                       /
                       forcing_frequency)
                    : 0.0
                );
                break;
            default :
                assert(deriv == Dissipation::NO_DERIV);
        }
        
        if(forcing_frequency == 0) {
            if(spin_frequency_multiplier >= 0) {
                above_lock_value = -result;
                return result;
            } else {
                above_lock_value = result;
                return -result;
            }
        } else {
            return (forcing_frequency > 0 ? result : -result);
        }

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

        if(__spin_frequency_breaks.size() != 0) {
            *result |= new LagSpinBreakCondition(
                *this,
                (primary ? system.primary() : system.secondary()),
                (primary ? system.secondary() : system.primary()),
                primary,
                zone_index
            );
        }

        if(__tidal_frequency_breaks.size() != 0)
            add_tidal_frequency_conditions(system,
                                           primary,
                                           zone_index, 
                                           *result);

        return result;
    }//End BrokenPowerlawPhaseLagZone::stopping_conditions definition.

    void BrokenPowerlawPhaseLagZone::change_e_order(
        unsigned new_e_order,
        BinarySystem &system, 
        bool primary,
        unsigned zone_index
    )
    {
        __tidal_indices.resize(5 * (new_e_order + 3));

        double orbital_frequency = get_orbital_frequency(system);

        std::vector< std::vector<double>::size_type >::iterator
            destination = (__tidal_indices.begin()
                           +
                           5 * (eccentricity_order() + 3));
        for(
            int mp = static_cast<int>(eccentricity_order()) + 3;
            mp <= static_cast<int>(new_e_order) + 2;
            ++mp
        ) {
            for(int m = -2; m <= 2; ++m) {
                *destination = get_tidal_index(
                    forcing_frequency(mp, m, orbital_frequency)
                );
                ++destination;
            }
        }

        DissipatingZone::change_e_order(new_e_order,
                                        system,
                                        primary,
                                        zone_index);
    }

} //End BinarySystem namespace.
