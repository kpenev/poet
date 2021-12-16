/**\file
 *
 * \brief Declare an interface for evaluating the expansion of the tidal
 * potential.
 *
 * \ingroup Evolve_group
 */

#ifndef __TIDAL_POTENTIAL_TERMS_H
#define __TIDAL_POTENTIAL_TERMS_H

#include <cmath>
#include <valarray>
#include "EccentricityExpansionCoefficients.h"
#include "../Core/Common.h"

//TODO: No longer need to worry about e-order

namespace Evolve {
    class LIB_PUBLIC TidalPotentialTerms {
    private:
        ///The eccentricity expansion of \f$p_{m,s}\f$.
        static EccentricityExpansionCoefficients __pms;

        ///\brief The constant coefficiients in \f$\mathcal{U}_{m,m'}\f$ of Lai
        ///(2012).
        ///
        ///The first index is m+2 (since m starts from -2) and the second index
        ///is m'/2+1 since the only allowed values are -2, 0 and 1.
        static const double __Umm_coef[][3];

        ///The inclination with which __Ummp was last filled.
        double __Ummp_inclination;

        ///The argument of periaspsis set by the last call to configure().
        double __arg_of_periapsis;

        std::valarray< std::valarray<double> >
            ///The \f$\mathcal{U}_{m,m'}\f$ quantities defined in Lai (2012).
            __Ummp,

            ///\brief The derivatives of the \f$\mathcal{U}_{m,m'}\f$
            ///quantities w.r.t. the inclination.
            __Ummp_deriv;

    public:

        TidalPotentialTerms();

        ///See EccentricityExpansionCoefficients::prepare()
        static void prepare(const std::string &tabulated_pms_fname,
                            double precision,
                            bool pre_load)
        {__pms.prepare(tabulated_pms_fname, precision, pre_load);}


        ///\brief The maximum orbital frequency multiplier to include in the
        ///potential Fourier expansion in order to achive a specified precision.
        ///
        ///The return value (call it \f$O\f$) is such that \f$O p_{m, O+1} <
        //\mathrm{precision}\f$. The reasoning is that if \f$p_{m,s}\f$ will
        ///decay substantially as s doubles.
        static unsigned required_expansion_order(
            ///The eccentricity at which tidal potential needs to be evaluated.
            double e
        )
        {
            return std::max(
                std::max(
                    __pms.required_expansion_order(e, -2),
                    __pms.required_expansion_order(e, 0)
                ),
                __pms.required_expansion_order(e, 2)
            );
        }

        ///Return the expansion precision target for the expansion terms.
        inline double get_expansion_precision() const
        {return __pms.get_expansion_precision();}

        ///\brief Return the range of eccentricities (min, max) over which an
        ///expansion going up to given max m' is valid and required.
        ///
        ///For eccentricities below the minimum (first value of returned pair)
        ///terms with frequency \f$m\Omega_\star-m'\Omega_{orb}\f$ can be
        ///excluded from the expansion for all m, without violating the
        ///precision requirement. For eccentricities at or above maximum (second
        ///value of returned pair) at least
        ///one \f$m\Omega_\star-m'\Omega_{orb}\f$ term must be included in the
        ///series in order to satisfy the precision requirement.
        static std::pair<double,double> get_expansion_range(int max_mp);

        ///Set the inclination relative to the orbit.
        void configure(double inclination, double arg_of_periapsis = 0);

        ///\brief Calculates \f$\sum_s W_{2,s}D_{m,s}(\Theta)p_{s,m'}\f$ (see
        ///documentation) and its derivatives w.r.t. e and \f$\Theta\f$.
        ///
        ///configure() should already have been called with the appropriate
        ///inclination and argument of periapsis.
        void operator()(
            ///The eccentricity.
            double e,

            ///The m index (spin freuqency multiplier).
            int m,

            ///The m' index (orbital frequency multiplier).
            int mp,

            ///Set to the undifferentiated value.
            std::complex<double> &no_deriv,

            ///Set to the inclination derivative.
            std::complex<double> &inclination_deriv,

            ///Set to the eccentricity_derivative.
            std::complex<double> &eccentricity_deriv
        ) const;

        ///\brief Return only the real parts of the complex version of the
        ///operator, since only the real part enters the tidal torque and power.
        void operator()(
            ///The eccentricity.
            double e,

            ///The m index.
            int m,

            ///The m' index.
            int mp,

            ///Set to the undifferentiated value.
            double &no_deriv,

            ///Set to the inclination derivative.
            double &inclination_deriv,

            ///Set to the eccentricity_derivative.
            double &eccentricity_deriv
        ) const;

        ///\brief Reads the interpolation data for \f$p_{m,s}\f$.
        ///
        ///See EccentricityExpansionCoefficients::prepare for description of the
        ///arguments.
        static void prepare_eccentricity_expansion(const std::string &fname,
                                                   double precision,
                                                   bool pre_load)
        {__pms.prepare(fname, precision, pre_load);}

        ///\brief The maximum eccentricity expansion order (orbital frequency
        ///multiplier for which the expansion is known.
        static unsigned max_expansion_order()
        {return __pms.max_expansion_order();}

        ///Provide direct access to the eccentircity expansion coefficients
        static const EccentricityExpansionCoefficients &
            expansion_coefficient_evaluator()
            {return __pms;}
    }; //End TidalPotentialTerms class.
} //End Evolve namespace

#endif
