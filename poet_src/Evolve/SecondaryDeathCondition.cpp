#define BUILDING_LIBRARY
#include "SecondaryDeathCondition.h"
#include "BinarySystem.h"

namespace Evolve {

    std::valarray<double> SecondaryDeathCondition::operator()(
            Core::EvolModeType
#ifndef NDEBUG
            evol_mode
#endif
            ,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv) const
    {
        assert(evol_mode == Core::BINARY);
        assert(orbit.size() == (1
                                +
                                3 * __system.number_zones()
                                -
                                __system.number_locked_zones()));
        assert(orbit.size() == derivatives.size());

        double min_semimajor = __system.minimum_separation(),
               semimajor = __system.semimajor(),
               e = __system.eccentricity(),
               dsemimajor_dt = derivatives[0],
               de_dt = derivatives[1];
        if(__system.number_locked_zones() == 0)
            dsemimajor_dt *= semimajor/(6.5 * orbit[0]);
        stop_deriv.resize(
            1,
            (
                (1.0 - e) * (
                    dsemimajor_dt * min_semimajor
                    -
                    semimajor * __system.minimum_separation(true)
                )
                -
                semimajor * min_semimajor * de_dt
            )
            /
            std::pow(min_semimajor, 2));
#ifdef VERBOSE_DEBUG
        std::cerr << "amin = " << min_semimajor
                  << ", a = " << semimajor
                  << ", da/dt = " << dsemimajor_dt
                  << ", e = " << e
                  << ", de/dt = " << de_dt
                  << ", d(death cond)/dt: " << stop_deriv[0] << std::endl;
#endif
        return std::valarray<double>(
            (semimajor * (1 - e) - min_semimajor) / min_semimajor,
            1
        );
    }

    void SecondaryDeathCondition::reached(short deriv_sign,
                                          unsigned index)
    {
        assert(index == 0);
        assert(deriv_sign == -1);

        StoppingCondition::reached(deriv_sign, index);
        __system.secondary_died();
    }

    std::string SecondaryDeathCondition::describe(int ) const
    {
        std::ostringstream description;
        description << "Semimajor axis crossing death boundary of "
                    << __system.minimum_separation();
        return description.str();
    }

} //End Evolve namespace.
