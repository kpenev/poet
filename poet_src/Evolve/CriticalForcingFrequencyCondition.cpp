/**\file
 * 
 * \brief Definitions of some of the methods of the 
 * CriticalForcingFrequencyCondition class.
 *
 * \ingroup Evolve_group
 */

#include "CriticalForcingFrequencyCondition.h"

namespace Evolve {

    CriticalForcingFrequencyCondition::CriticalForcingFrequencyCondition(
        const DissipatingBody &body,
        const DissipatingBody &other_body,
        bool primary,
        unsigned zone_index,
        int orbital_frequency_multiplier,
        int spin_frequency_multiplier,
        std::vector<double> critical_frequencies,
        double orbital_frequency
    ) :
        __critical_frequencies(critical_frequencies),
        __zone(body.zone(zone_index)),
        __critical_above_iter(
            std::lower_bound(
                __critical_frequencies.begin(),
                __critical_frequencies.end(),
                __zone.forcing_frequency(orbital_frequency_multiplier,
                                         spin_frequency_multiplier,
                                         orbital_frequency)
            )
        ),
        __critical_below_iter(
            __critical_above_iter == __critical_spins.begin()
            ? __critical_spins.end()
            : __critical_above_iter - 1
        ),
        __body(body),
        __other_body(other_body),
        __primary(primary),
        __zone_index(zone_index)
    {
        if(
            __zone.forcing_frequency(
                orbital_frequency_multiplier,
                spin_frequency_multiplier,
                orbital_frequency
            ) == *critical_above_iter
        )
            throw Core::Error::BadFunctionArguments(
                "Starting evolution from exactly a critical spin frequency "
                "is not currently supported."
            )
    }

    std::valarray<double> CriticalForcingFrequencyCondition::operator()(
        Core::EvolModeType evol_mode,
        const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
        double semimajor;
        if(
            __body.number_locked_zones()
            ||
            __other_body.number_locked_zones()
        )
            semimajor = parameters[0];
        else if(parameters[0] < 0)
            throw Core::Error::Runtime("Negitave a^6.5 encountered!");
        else semimajor = std::pow(parameters[0], 1.0 / 6.5);

        double orbital_frequency = orbital_angular_velocity(
            body.mass(),
            other_body.mass(),
            semimajor
        );

        double forcing_frequency = __zone.forcing_frequency(
            __orbital_frequency_multiplier,
            __spin_frequency_multiplier,
            orbital_frequency
        );
        std::valarray<double> result(2);
        result[0] = ((forcing_frequency - *__critical_below_iter)
                     /
                     *__critical_below_iter);
        result[1] = ((forcing_frequency - *__critical_above_iter)
                     /
                     *__critical_above_iter);
        stop_deriv.resize(2, Core::NaN);
        return result;
    }
}//End Evolve namespace.
