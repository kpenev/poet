#ifndef __EVOLVING_STAR_H
#define __EVOLVING_STAR_H

/**\file
 *
 * \brief Declares the class for stars that user pre-tabulated stellar
 * evolution tracks.
 *
 * \ingroup Star_group
 */

#include "SaturatingSkumanichWindBody.h"
#include "ExponentialDecayDiffRotBody.h"
#include "../StellarEvolution/Interpolator.h"
#include "EvolvingStellarCore.h"
#include "EvolvingStellarEnvelope.h"

namespace Star {

    class InterpolatedEvolutionStar : 
        public SaturatingSkumanichWindBody,
        public ExponentialDecayDiffRotBody {
    private:
        ///\brief The luminosity of the star in \f$L_\odot\f$ as a function of
        ///age in Gyr.
        ///
        ///Since the luminosity does not enter in the equations solved it
        ///might not be defined. If that is the case, this member should be
        ///NULL.
        const StellarEvolution::EvolvingStellarQuantity *__luminosity;

        double 
            ///The age at which the star leaves the main sequence in Gyrs.
            __lifetime;

        ///The surface zone of the star (the entire star if high mass).
        EvolvingStellarEnvelope __envelope;

        ///The core of the star (NULL if high mass).
        EvolvingStellarCore __core;

    public:
        InterpolatedEvolutionStar(
                ///Mass of the star
                double mass,

                ///The metallicity ([Fe/H]) of the star
                double metallicity,

                ///The strength of the wind.
                double wind_strength,

                ///The frequency at which the wind loss saturates in rad/day.
                double wind_saturation_frequency,

                ///The timescale for differential rotation coupling.
                double diff_rot_coupling_timescale,

                ///A StellarEvolution interpolator.
                const StellarEvolution::Interpolator &interpolator
        ) : 
            SaturatingSkumanichWindBody(wind_strength,
                                        wind_saturation_frequency),
            ExponentialDecayDiffRotBody(diff_rot_coupling_timescale),
            __luminosity(interpolator(StellarEvolution::LUM, mass, metallicity)),
            __lifetime(__luminosity->range_high()),
            __envelope(mass,
                       interpolator(StellarEvolution::RADIUS, mass, metallicity),
                       interpolator(StellarEvolution::ICONV, mass, metallicity)),
            __core(interpolator.core_formation_age(), 
                   interpolator(StellarEvolution::MRAD, mass, metallicity),
                   interpolator(StellarEvolution::RRAD, mass, metallicity),
                   interpolator(StellarEvolution::IRAD, mass, metallicity))
            {}

        ///The number of zones the body consists of.
        unsigned number_zones() const {return 2;}

        ///See Evolve::DissipatingBody::zone().
        Evolve::DissipatingZone &zone(unsigned zone_index)
        {
#ifdef DEBUG
            assert(zone_index<=1);
#endif
            if(zone_index == 0) return __envelope;
            else return __core;
        }

        ///The envelope of the star - inmodifiable.
        const EvolvingStellarEnvelope &envelope() const {return __envelope;}

        ///The envelope of the star - modifiable.
        EvolvingStellarEnvelope &envelope() {return __envelope;}

        ///The core of the star.
        const EvolvingStellarCore &core() const {return __core;}

        ///The core of the star.
        EvolvingStellarCore &core() {return __core;}

        ///See Evolve::DissipatingBody::zone().
        const Evolve::DissipatingZone &zone(unsigned zone_index) const
        {
#ifdef DEBUG
            assert(zone_index<=1);
#endif
            if(zone_index==0) return __envelope;
            else return __core;
        }

        ///The lifetime of the star (where tracks end).
        double lifetime() const {return __lifetime;}

        ///The luminosity of the star at the given age.
        double luminosity(double age) const {return (*__luminosity)(age);}

        ///Cleanup after the star.
        ~InterpolatedEvolutionStar() {delete __luminosity;}

        ///\brief Change the star as necessary at the given age.
        ///
        ///Handles things like interpolation discontinuities. 
        virtual void reached_critical_age(double age)
        {
            __core.reached_critical_age(age);
            __envelope.reached_critical_age(age);
        }

        ///\brief The next age when the evolution needs to be stopped for a
        ///change in one of the bodies.
        virtual double next_stop_age() const
        {return std::min(__core.next_stop_age(), __envelope.next_stop_age());}
    };//End InterpolatedEvolutionStar class

}//End Star namespace.

#endif
