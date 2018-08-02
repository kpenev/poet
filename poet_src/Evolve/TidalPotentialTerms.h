/**\file
 *
 * \brief Declare an interface for evaluating the expansion of the tidal
 * potential.
 *
 * \ingroup Evolve_group
 */

#include <cmath>
#include <valarray>
#include "EccentricityExpansionCoefficients.h"
#include "../Core/Common.h"

namespace Evolve {
    class LIB_PUBLIC TidalPotentialTerms {
    private:
        ///The expansion order in eccentricity to use.
        unsigned __e_order;

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

        std::valarray< std::valarray<double> >
            ///The \f$\mathcal{U}_{m,m'}\f$ quantities defined in Lai (2012).
            __Ummp,

            ///\brief The derivatives of the \f$\mathcal{U}_{m,m'}\f$
            ///quantities w.r.t. the inclination.
            __Ummp_deriv;

    public:

        TidalPotentialTerms();

        ///Change the eccentricity expansion order.
        void change_e_order(unsigned new_e_order)
        {__e_order = new_e_order;}

        ///Set the inclination relative to the orbit.
        void configure(double inclination);

        ///\brief Calculates \f$\sum_s W_{2,s}D_{m,s}(\Theta)p_{s,m'}\f$ (see
        ///documentation) and its derivatives w.r.t. e and \f$\Theta\f$.
        ///
        ///fill_Umm should already have been called with the appropriate
        ///inclination.
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
        );

        ///\brief Reads the eccentricity expansion coefficients of \f$p_{m,s}\f$.
        ///
        ///The given file should have been generated by
        ///tabulate_eccentricity_expansion_coefficients.py.
        static void read_eccentricity_expansion(const std::string &fname)
        {__pms.read(fname);}
    }; //End TidalPotentialTerms class.
} //End Evolve namespace
