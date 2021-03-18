/**\file
 *
 * \brief Declare a dissipative zone class with polynomial evolution with a
 * constant phase lag
 *
 * \ingroup UnitTests_group
 */

#ifndef __CONST_PHASE_LAG_H
#define __CONST_PHASE_LAG_H

#include "PolynomialEvolutionZone.h"

namespace Evolve {

    /**\brief A zone with constant phase lag for all tidal terms.
     *
     * \ingroup UnitTests_group
     */
    class ConstPhaseLagZone : public PolynomialEvolutionZone {
    private:
        ///The constant value of the phase lag.
        double __phase_lag;
    public:
        ConstPhaseLagZone(
            ///The phase lag of all dissipative terms.
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
            __phase_lag(phase_lag)
        {}

        ///See DissipatingZone::modified_phase_lag()
        double modified_phase_lag(int,
                                  int spin_frequency_multiplier,
                                  double forcing_frequency,
                                  Dissipation::QuantityEntry entry,
                                  double &above_lock_value) const
        {
            double result;

            if(entry != Dissipation::NO_DERIV)
                result = 0.0;
            else
                result = __phase_lag;

            if(forcing_frequency == 0) {
                return 0.0;
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
        }

        ///See DissipatingZone::dissipative()
        bool dissipative() const {return true;}

        ///See DissipatingZone::can_lock()
        bool can_lock() const {return true;}
    };//End ConstPhaseLagZone class.
}//End Evolve namespace.

#endif
