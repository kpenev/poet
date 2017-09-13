/**\file
 * 
 * \brief Definitions of some of the methods of the 
 * LagForcingFrequencyBreakCondition class.
 *
 * \ingroup Evolve_group
 */

#include "LagForcingFrequencyBreakCondition.h"
#include "BrokenPowerlawPhaseLagZone.h"

namespace Evolve {

    void LagForcingFrequencyBreakCondition::set_num_subconditions()
    {
        if(
            __powerlaw_index == 0
            ||
            __powerlaw_index == __zone.__tidal_frequency_breaks.size()
        )
            __num_subconditions = 1;
        else
            __num_subconditions = 2;
    }

    LagForcingFrequencyBreakCondition::LagForcingFrequencyBreakCondition(
        BrokenPowerlawPhaseLagZone &zone,
        const DissipatingBody &body,
        const DissipatingBody &other_body,
        bool primary,
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier
    ) :
        __orbital_frequency_multiplier(orbital_frequency_multiplier),
        __spin_frequency_multiplier(spin_frequency_multiplier),
        __zone(zone),
        __body(body),
        __other_body(other_body),
        __primary(primary),
        __term_index(
            __zone.tidal_term_index(orbital_frequency_multiplier,
                                    spin_frequency_multiplier)
        ),
        __powerlaw_index(__zone.__tidal_indices[__term_index])
    {
        set_num_subconditions();
    }

    std::valarray<double> LagForcingFrequencyBreakCondition::operator()(
        Core::EvolModeType evol_mode,
        const std::valarray<double> &orbit,
        const std::valarray<double> &,
        std::valarray<double> &stop_deriv
    ) const
    {
        assert(evol_mode == Core::BINARY);
        double semimajor;
        if(
            __body.number_locked_zones()
            ||
            __other_body.number_locked_zones()
        )
            semimajor = orbit[0];
        else if(orbit[0] < 0)
            throw Core::Error::Runtime("Negitave a^6.5 encountered!");
        else semimajor = std::pow(orbit[0], 1.0 / 6.5);

        double orbital_frequency = Core::orbital_angular_velocity(
            __body.mass(),
            __other_body.mass(),
            semimajor
        );

        double abs_forcing_frequency = std::abs(
            __zone.forcing_frequency(
                __orbital_frequency_multiplier,
                __spin_frequency_multiplier,
                orbital_frequency
            )
        );
        std::valarray<double> result(__num_subconditions);
        size_t above_index;
        if(__powerlaw_index == 0)
            above_index = 0;
        else {
            double critical_frequency =
                __zone.__tidal_frequency_breaks[__powerlaw_index - 1];

            result[0] = ((abs_forcing_frequency - critical_frequency)
                         /
                         critical_frequency);
            above_index = 1;
        }
        if(__powerlaw_index < __zone.__tidal_frequency_breaks.size()) {
            double critical_frequency = 
                __zone.__tidal_frequency_breaks[__powerlaw_index];

            result[above_index] = (
                (abs_forcing_frequency - critical_frequency)
                /
                critical_frequency
            );
        }
        stop_deriv.resize(__num_subconditions, Core::NaN);
        return result;
    }

    void LagForcingFrequencyBreakCondition::reached(short deriv_sign,
                                                    unsigned index)
    {
        assert(index < __num_subconditions);

        if(__powerlaw_index > 0 && index == 0) {
            assert(deriv_sign == -1);
            --__powerlaw_index;
            --__zone.__tidal_indices[__term_index];
        } else {
            assert(deriv_sign == 1);
            assert(__powerlaw_index < __zone.__tidal_frequency_breaks.size());

            ++__powerlaw_index;
            ++__zone.__tidal_indices[__term_index];

            assert(__zone.__tidal_indices[__term_index]
                   <
                   __zone.__tidal_frequency_powers.size());
        }

        set_num_subconditions();

    }

    short LagForcingFrequencyBreakCondition::expected_crossing_deriv_sign(
        unsigned index
    ) const
    {
        if(index == 1 || __powerlaw_index == 0) {
            assert(__powerlaw_index
                   <
                   __zone.__tidal_frequency_breaks.size());
            return 1;
        } else {
            return -1;
        }
    }

    std::string LagForcingFrequencyBreakCondition::describe(int index) const
    {
        std::ostringstream description;
        description << "Forcing frequency leaving the interval ";
        if(__powerlaw_index > 0 && index <= 0) {
            description
                << __zone.__tidal_frequency_breaks[__powerlaw_index - 1]
                << " < ";
        }
        description << "|"
                    << __orbital_frequency_multiplier
                    << "Worb"
                    << (__spin_frequency_multiplier > 0 ? " - " : " + ")
                    << std::abs(__spin_frequency_multiplier)
                    << "Wspin|";
        if(
            __powerlaw_index < __zone.__tidal_frequency_breaks.size()
            &&
            (
                index < 0
                || 
                index == __num_subconditions - 1
            )
        ) {
            description << " < "
                        << __zone.__tidal_frequency_breaks[__powerlaw_index];
        }
        return description.str();
    }

}//End Evolve namespace.
