/**\file
 *
 * \brief Declare a dissipative zone class with polynomial evolution with only a
 * single tidel term having non-zero dissipation.
 *
 * \ingroup UnitTests_group
 */

#ifndef __SINGLE_TIDAL_TERM_ZONE_H
#define __SINGLE_TIDAL_TERM_ZONE_H

#include "PolynomialEvolutionZone.h"

namespace Evolve {

    /**\brief A zone dissipative to only a single tidal term.
     *
     * \ingroup UnitTests_group
     */
    class SingleTidalTermZone : public PolynomialEvolutionZone {
    private:
        int
            ///See orbital_frequency_multiplier argument to constructor.
            __orbital_frequency_multiplier,

            ///See spin_frequency_multiplier argument to constructor.
            __spin_frequency_multiplier;

        ///See phase_lag argument to constructor.
        double __phase_lag;
    public:
        SingleTidalTermZone(
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
            PolynomialEvolutionZone(mass_coefficients,
                                    radius_coefficients,
                                    inertia_coefficients),
            __orbital_frequency_multiplier(orbital_frequency_multiplier),
            __spin_frequency_multiplier(spin_frequency_multiplier),
            __phase_lag(phase_lag)
        {}

        ///See DissipatingZone::modified_phase_lag()
        double modified_phase_lag(int orbital_frequency_multiplier,
                                  int spin_frequency_multiplier,
                                  double forcing_frequency,
                                  Dissipation::QuantityEntry entry,
                                  double &above_lock_value) const
        {
            double result;

            if(entry != Dissipation::NO_DERIV)
                result = 0.0;
            else {
                if(
                    orbital_frequency_multiplier == __orbital_frequency_multiplier
                    &&
                    spin_frequency_multiplier == __spin_frequency_multiplier
                )
                    result = __phase_lag;
                else if(
                    orbital_frequency_multiplier == -__orbital_frequency_multiplier
                    &&
                    spin_frequency_multiplier == -__spin_frequency_multiplier

                )
                    result = -__phase_lag;
                else
                    result = 0.0;
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
                return result;
            }
        }

        ///See DissipatingZone::dissipative()
        bool dissipative() const {return true;}

        ///See DissipatingZone::can_lock()
        bool can_lock() const {return true;}

    };//End SingleTidalTermZone class.

}//End Evolve namespace.

#endif
