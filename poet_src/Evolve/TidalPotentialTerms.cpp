/**\file
 *
 * \brief Define the methods of TidalPotentialTerms.
 *
 * \ingroup Evolve_group
 */

#include "TidalPotentialTerms.h"

namespace Evolve {

    EccentricityExpansionCoefficients TidalPotentialTerms::__pms;

    const double TidalPotentialTerms::__Umm_coef[][3]={
        {
            std::sqrt(3.0 * M_PI / 10.0) / 4.0,
            -std::sqrt(6.0 * M_PI / 5.0) / 4.0,
            std::sqrt(3.0 * M_PI / 10.0) / 4.0
        },
        {
            -std::sqrt(3.0 * M_PI / 10.0) / 2.0,
            -std::sqrt(6.0 * M_PI / 5.0) / 2.0,
            std::sqrt(3.0 * M_PI / 10.0) / 2.0
        },
        {
            3.0 * std::sqrt(M_PI / 5.0) / 4.0,
            -std::sqrt(M_PI / 5.0) / 2.0,
            3.0 * std::sqrt(M_PI / 5.0) / 4.0
        },
        {
            -std::sqrt(3.0 * M_PI / 10.0) / 2.0,
            std::sqrt(6.0 * M_PI / 5.0) / 2.0,
            std::sqrt(3.0 * M_PI / 10.0) / 2.0
        },
        {
            std::sqrt(3.0 * M_PI / 10.0) / 4.0,
            -std::sqrt(6.0 * M_PI / 5.0) / 4.0,
            std::sqrt(3.0 * M_PI / 10.0) / 4.0
        }
    };

    void TidalPotentialTerms::configure(double inclination,
                                        double arg_of_periapsis)
    {
        __arg_of_periapsis = arg_of_periapsis;

        if(__Ummp_inclination==inclination) return;
        __Ummp_inclination=inclination;
        double c=std::cos(__Ummp_inclination), s=std::sin(__Ummp_inclination),
               s2=std::pow(s, 2), sc=s*c, cp1=c+1.0, cm1=c-1.0;

        __Ummp[0][0]=__Umm_coef[0][0]*std::pow(cp1, 2);
        __Ummp_deriv[0][0]=-__Umm_coef[0][0]*2.0*s*cp1;

        __Ummp[1][0]=__Umm_coef[1][0]*s*cp1;
        __Ummp_deriv[1][0]=__Umm_coef[1][0]*(cp1+2.0*s2);

        __Ummp[2][0]=__Umm_coef[2][0]*s2;
        __Ummp_deriv[2][0]=__Umm_coef[2][0]*2.0*sc;

        __Ummp[3][0]=-__Umm_coef[3][0]*s*cm1;
        __Ummp_deriv[3][0]=-__Umm_coef[3][0]*(c*cm1-s2);

        __Ummp[4][0]=__Umm_coef[4][0]*std::pow(cm1, 2);
        __Ummp_deriv[4][0]=-__Umm_coef[4][0]*2.0*s*cm1;



        __Ummp[0][1]=__Umm_coef[0][1]*s2;
        __Ummp_deriv[0][1]=__Umm_coef[0][1]*2.0*sc;

        __Ummp[1][1]=__Umm_coef[1][1]*sc;
        __Ummp_deriv[1][1]=__Umm_coef[1][1]*(1.0-2.0*s2);

        __Ummp[2][1]=__Umm_coef[2][1]*(2.0-3.0*s2);
        __Ummp_deriv[2][1]=-__Umm_coef[2][1]*6.0*sc;

        __Ummp[3][1]=__Umm_coef[3][1]*sc;
        __Ummp_deriv[3][1]=__Umm_coef[3][1]*(1.0-2.0*s2);

        __Ummp[4][1]=__Umm_coef[4][1]*s2;
        __Ummp_deriv[4][1]=__Umm_coef[4][1]*2.0*sc;



        __Ummp[0][2]=__Umm_coef[0][2]*std::pow(cm1, 2);
        __Ummp_deriv[0][2]=-__Umm_coef[0][2]*2.0*cm1*s;

        __Ummp[1][2]=-__Umm_coef[1][2]*s*cm1;
        __Ummp_deriv[1][2]=-__Umm_coef[1][2]*(c*cm1-s2);

        __Ummp[2][2]=__Umm_coef[2][2]*s2;
        __Ummp_deriv[2][2]=__Umm_coef[2][2]*2.0*sc;

        __Ummp[3][2]=__Umm_coef[3][2]*s*cp1;
        __Ummp_deriv[3][2]=__Umm_coef[3][2]*(c*cp1-s2);

        __Ummp[4][2]=__Umm_coef[4][2]*std::pow(cp1, 2);
        __Ummp_deriv[4][2]=-__Umm_coef[4][2]*2.0*cp1*s;
    }

    TidalPotentialTerms::TidalPotentialTerms() :
        __e_order(0),
        __Ummp_inclination(Core::NaN),
        __Ummp(5),
        __Ummp_deriv(5)
    {
        for(int i=0; i<5; ++i) {
            __Ummp[i].resize(3);
            __Ummp_deriv[i].resize(3);
        }
    }

    void TidalPotentialTerms::operator()(
        double e,
        int m,
        int mp, 
        std::complex<double> &no_deriv,
        std::complex<double> &inclination_deriv,
        std::complex<double> &eccentricity_deriv
    ) const
    {
        no_deriv=inclination_deriv=eccentricity_deriv=0;
        for(int i=0; i<3; ++i) {
            int s = 2 * (i - 1);
            double pms=__pms(s, mp, e, __e_order, false);
            std::complex<double> periapsis_factor(
                std::cos(s * __arg_of_periapsis),
                -std::sin(s * __arg_of_periapsis)
            );
            no_deriv += pms * __Ummp[m+2][i] * periapsis_factor;
            inclination_deriv += pms * __Ummp_deriv[m+2][i] * periapsis_factor;
            eccentricity_deriv += (
                __pms(2*(i-1), mp, e, __e_order, true)
                *
                __Ummp[m+2][i]
                *
                periapsis_factor
            );
        }
    }

    void TidalPotentialTerms::operator()(double e,
                                         int m,
                                         int mp, 
                                         double &no_deriv,
                                         double &inclination_deriv,
                                         double &eccentricity_deriv) const
    {
        std::complex<double> complex_no_deriv,
                             complex_inclination_deriv,
                             complex_eccentricity_deriv;
        operator()(e,
                   m,
                   mp,
                   complex_no_deriv,
                   complex_inclination_deriv,
                   complex_eccentricity_deriv);
        no_deriv = std::abs(complex_no_deriv);
        inclination_deriv = (
            complex_no_deriv.real() * complex_inclination_deriv.real()
            +
            complex_no_deriv.imag() * complex_inclination_deriv.imag()
        ) / no_deriv;
        eccentricity_deriv = (
            complex_no_deriv.real() * complex_eccentricity_deriv.real()
            +
            complex_no_deriv.imag() * complex_eccentricity_deriv.imag()
        ) / no_deriv;
    }
} //End Evolve namespace.
