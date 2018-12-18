/**\file
 *
 * \brief Declare a dissipative body class with a single SingleTidalTermZone.
 *
 * \ingroup UnitTests_group
 */

#ifndef __SINGLE_TIDAL_TERM_BODY_H
#define __SINGLE_TIDAL_TERM_BODY_H

#include "../../Star/SaturatingSkumanichWindBody.h"
#include "SingleTidalTermZone.h"

namespace Evolve {

    /**\brief A skumanich wind body with a single zone dissipative to only a
     * single tidal term.
     *
     * \ingroup UnitTests_group
     */
    class SingleTidalTermBody : public Star::SaturatingSkumanichWindBody {
    private:
        SingleTidalTermZone __zone;
    public:
        ///\brief Initialize the body's single zone and the wind with the given
        ///arguments.
        SingleTidalTermBody(
            ///The strength of the wind.
            double wind_strength,

            ///The frequency at which the wind loss saturates in rad/day.
            double saturation_frequency,

            ///The multiplier of the orbital frequency of the only dissipative
            ///term.
            int orbital_frequency_multiplier,

            ///The multiplier of the spin frequency of the only dissipative
            ///term.
            int spin_frequency_multiplier,

            ///The phase lag of the only dissipative term.
            double phase_lag,

            ///The coefficients of the polynomial for the mass evolution.
            const std::valarray<double> &mass_coefficients,

            ///The coefficients of the polynomial for the radius evolution.
            const std::valarray<double> &radius_coefficients,

            ///The coefficients of the polynomial for the moment of inertia
            ///evolution.
            const std::valarray<double> &inertia_coefficients
        ) :
            SaturatingSkumanichWindBody(wind_strength,
                                        saturation_frequency),
            __zone(orbital_frequency_multiplier,
                   spin_frequency_multiplier,
                   phase_lag,
                   mass_coefficients,
                   radius_coefficients,
                   inertia_coefficients)
        {}

        ///See DissipatingBody::number_zones().
        unsigned number_zones() const {return 1;}

        ///See DissipatingBody::zone(int) const.
        const DissipatingZone &zone(
            unsigned
#ifndef NDEBUG
            zone_index
#endif
        ) const
        {
            assert(zone_index == 0);
            return __zone;
        }

        ///See DissipatingBody::zone(int).
        DissipatingZone &zone(
            unsigned
#ifndef NDEBUG
            zone_index
#endif
        )
        {
            assert(zone_index == 0);
            return __zone;
        }

        ///See DissipatingBody::angular_momentum_coupling().
        Eigen::Vector3d angular_momentum_coupling(
            unsigned,
            Dissipation::QuantityEntry =Dissipation::NO_DERIV,
            bool =false
        ) const
        {assert(false);}
    };//End SingleTidalTermBody class.

}//End Evolve namespace.

#endif
