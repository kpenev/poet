#define BUILDING_LIBRARY
#include "SynchronizedCondition.h"
#include "BinarySystem.h"
#include "DissipatingZone.h"

namespace Evolve {

    SynchronizedCondition::SynchronizedCondition(int orbital_freq_mult,
                                                 int spin_freq_mult,
                                                 short deriv_sign,
                                                 bool primary,
                                                 unsigned zone_index,
                                                 BinarySystem &system) :
        StoppingCondition(deriv_sign),
        __orbital_freq_mult(orbital_freq_mult),
        __spin_freq_mult(spin_freq_mult),
        __primary(primary),
        __zone_index(zone_index),
        __zone((primary
                ? system.primary()
                : system.secondary()).zone(zone_index)),
        __system(system)
    {}

    std::valarray<double> SynchronizedCondition::operator()(
        Core::EvolModeType
#ifndef NDEBUG
        evol_mode
#endif
        , const std::valarray<double> &orbit,
        const std::valarray<double> &derivatives,
        std::valarray<double> &stop_deriv
    ) const
    {
#ifndef NDEBUG
        assert(evol_mode == Core::BINARY);
        if(__system.number_locked_zones())
            assert(orbit[0] == __system.semimajor());
        else assert(std::pow(std::max(0.0, orbit[0]), 1.0 / 6.5)
                    ==
                    __system.semimajor());
        assert(orbit.size() == (1
                                +
                                3 * __system.number_zones()
                                -
                                __system.number_locked_zones()));
        assert(orbit.size() == derivatives.size());
#endif
        double
               m2 = __system.secondary().mass(),
               semimajor = __system.semimajor(),
               worb = __system.orbital_frequency(),
               wspin = __zone.spin_frequency(),
               dworb_dt = (__system.orbital_frequency(true)
                           *
                           derivatives[0]
                           *
                           Core::AstroConst::solar_radius);


        unsigned angmom_ind = 1 + 2 * __system.number_zones();
        if(!__primary) angmom_ind += (__system.primary().number_zones()
                                      -
                                      __system.primary().number_locked_zones());

        const DissipatingBody &body = (__primary
                                       ? __system.primary()
                                       : __system.secondary());
        for(unsigned i = 0; i < __zone_index; ++i)
            if(!body.zone(i).locked()) ++angmom_ind;

        double dwspin_dt = (
            (derivatives[angmom_ind] - __zone.moment_of_inertia(1) * wspin)
            /
            __zone.moment_of_inertia()
        );
        if(__system.number_locked_zones() == 0)
            dworb_dt /= 6.5 * orbit[0] / semimajor;

        stop_deriv.resize(
            1,
            (
                (wspin * dworb_dt - dwspin_dt * worb)
                /
                std::pow(worb, 2)
                *
                __spin_freq_mult
            )
        );

#ifdef VERBOSE_DEBUG
        std::cerr << describe() << " angmom index: " << angmom_ind
            << " worb = " << worb
            << " dworb_dt = " << dworb_dt
            << " wspin = " << wspin
            << " dwspin_dt = " << dwspin_dt
            << " adot = " << derivatives[0] / (6.5 * orbit[0] / semimajor)
            << " has deriv = " << stop_deriv[0]
            << std::endl;
#endif


#ifndef NDEBUG
        if(
            std::isnan((__orbital_freq_mult * worb - wspin * __spin_freq_mult)
                       /
                       worb)
        ) {
            std::cerr << "Synchronized value for zone " << __zone_index << "/"
                << __system.primary().number_zones() << " for worb=" << worb
                << ", wspin=" << wspin
                << ", m'=" << __orbital_freq_mult << ", m="
                << "__spin_freq_mult is NaN."
                << ", semimajor=" << semimajor << std::endl ;
            assert(false);
        }
#endif
        return std::valarray<double>(
            (__orbital_freq_mult * worb - wspin * __spin_freq_mult) / worb,
            1
        );
    }

    void SynchronizedCondition::reached(short deriv_sign, unsigned index)
    {
#ifndef NDEBUG
        std::cerr << "Synchronization reached: "
                  << "Expected sign: " << expected_crossing_deriv_sign()
                  << "sign: " << deriv_sign
                  << std::endl;
#endif
        StoppingCondition::reached(deriv_sign, index);
#ifndef NDEBUG
        std::cerr << "Now expected sign: " << expected_crossing_deriv_sign()
                  << std::endl;
#endif
    }

    std::string SynchronizedCondition::describe(int ) const
    {
        std::ostringstream description;
        description << (__primary ? "Primary" : "Secondary")
                    << " body, zone "
                    << __zone_index
                    << " satisfying "
                    << __orbital_freq_mult
                    << "(orbital frequency) = "
                    << __spin_freq_mult
                    << "(spin frequency)"
                    << ", expected crossing sign: "
                    << expected_crossing_deriv_sign();
        return description.str();
    }

} //End Evolve namespace.
