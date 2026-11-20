/**\file
 *
 * \brief Declares the class for stars that user pre-tabulated stellar
 * evolution tracks.
 *
 * \ingroup Star_group
 */

#ifndef __EVOLVING_STAR_H
#define __EVOLVING_STAR_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SaturatingSkumanichWindBody.h"
#include "ExponentialDecayDiffRotBody.h"
#include "../StellarEvolution/Interpolator.h"
#include "EvolvingStellarCore.h"
#include "EvolvingStellarEnvelope.h"
#include <memory>

namespace Star {

    class LIB_PUBLIC InterpolatedEvolutionStar : 
        public SaturatingSkumanichWindBody,
        public ExponentialDecayDiffRotBody {
    private:
        ///\brief The luminosity of the star in \f$L_\odot\f$ as a function
        ///of age in Gyr.
        ///
        ///Since the luminosity does not enter in the equations solved it
        ///might not be defined. If that is the case, this member should be
        ///NULL.
        const StellarEvolution::EvolvingStellarQuantity *__luminosity;

        ///The age at which the star leaves the main sequence in Gyrs.
        double __lifetime;

        ///The surface zone of the star (the entire star if high mass).
        EvolvingStellarEnvelope __envelope;

        ///The core of the star (NULL if high mass).
        EvolvingStellarCore __core;

        // When present, the main track describes the envelope and radiative
        // shell, while this second track supplies the convective core.
        std::unique_ptr<EvolvingStellarCore> __convective_core;

    public:
        InterpolatedEvolutionStar(
            ///Mass of the star
            double mass,

            ///The [Fe/H] of the star
            double feh,

            ///The strength of the wind.
            double wind_strength,

            ///The frequency at which the wind loss saturates in rad/day.
            double wind_saturation_frequency,

            ///The timescale for differential rotation coupling.
            double diff_rot_coupling_timescale,

            ///A StellarEvolution interpolator.
            const StellarEvolution::Interpolator &interpolator,
            const StellarEvolution::Interpolator *core_interpolator = NULL
        ) : 
            SaturatingSkumanichWindBody(wind_strength,
                                        wind_saturation_frequency),
            ExponentialDecayDiffRotBody(diff_rot_coupling_timescale),
            __luminosity(interpolator(StellarEvolution::LUM, mass, feh)),
            __lifetime(__luminosity->range_high()),
            __envelope(mass,
                       interpolator(StellarEvolution::RADIUS,
                                    mass,
                                    feh),
                       interpolator(StellarEvolution::ICONV,
                                    mass,
                                    feh)),
            __core(std::max(__envelope.min_interp_age(),
                            interpolator.core_formation_age()),
                   interpolator(StellarEvolution::MRAD, mass, feh),
                   interpolator(StellarEvolution::RRAD, mass, feh),
                   interpolator(StellarEvolution::IRAD, mass, feh))
            {
                if(core_interpolator) {
                    __convective_core.reset(new EvolvingStellarCore(
                        std::max(__envelope.min_interp_age(),
                                 core_interpolator->core_formation_age()),
                        (*core_interpolator)(StellarEvolution::MRAD, mass, feh),
                        (*core_interpolator)(StellarEvolution::RRAD, mass, feh),
                        (*core_interpolator)(StellarEvolution::IRAD, mass, feh)
                    ));
                }
            }

        ///The number of zones the body consists of.
        unsigned number_zones() const {return __convective_core ? 3 : 2;}

        ///See Evolve::DissipatingBody::zone().
        Evolve::DissipatingZone &zone(unsigned zone_index)
        {
            assert(zone_index < number_zones());

            if(zone_index == 0) return __envelope;
            if(zone_index == 1) return __core;
            return *__convective_core;
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
            assert(zone_index < number_zones());

            if(zone_index==0) return __envelope;
            if(zone_index==1) return __core;
            return *__convective_core;
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
            Evolve::DissipatingBody::reached_critical_age(age);
        }

        ///\brief The next age when the evolution needs to be stopped for a
        ///change in one of the bodies.
        virtual double next_stop_age() const
        {return Evolve::DissipatingBody::next_stop_age();}

        ///\brief Prepare the stellar quantities for interpolation around the
        ///given age.
        ///
        ///After calling this method, requesting values or derivatives
        ///outside the range of the continuous region containing this age
        ///(see ::discontinuities) fails an assert.
        virtual void select_interpolation_region(double age) const
        {
            __core.select_interpolation_region(age);
            __envelope.select_interpolation_region(age);
            if(__convective_core)
                __convective_core->select_interpolation_region(age);
        }
    };//End InterpolatedEvolutionStar class

}//End Star namespace.

#endif
