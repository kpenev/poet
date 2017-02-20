/**\file
 *
 * \brief Declares a class which provides the
 * [\f$p_{m,s}\f$ coefficients]{@ref InclinationEccentricity_pms1}.
 */


#ifndef __ECCENTRICITY_EXPANSION_COEFFICIENTS_H
#define __ECCENTRICITY_EXPANSION_COEFFICIENTS_H

#include "../Core/Error.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>

namespace Evolve {

    ///\brief A class which reads-in and provides a convenient interface to the 
    /// \f$p_{m,s}\f$ coefficients
    class EccentricityExpansionCoefficients {
    private:
        ///Maximum eccentricity power with all necessary coefficients known.
        unsigned __max_e_power;

        std::vector< std::vector<double> >
            ///\brief The expansion coefficients of \f$p_{0,s}\f$
            ///(\f$\alpha\f$ along with the s factior in the
            /// [documentation]{@ref InclinationEccentricity_pms1}).
            ///
            ///The outer index is s+__max_e_power and the inner one is 
            ///min(n, n+s).
            __alpha,

            ///\brief The expansion coefficients of \f$p_{2,s}\f$
            ///(\f$\gamma_{s,n}^+(s/2)^{s+2n}\f$ in the
            /// [documentation]{@ref InclinationEccentricity_pms1}).
            ///
            ///The outer index is s+__max_e_power-2 and the inner one is 
            ///min(n+1, n+s-1).
            __gamma_plus,

            ///\brief The expansion coefficients of \f$p_{-2,s}\f$
            ///(\f$\gamma_{s,n}^-(s/2)^{s+2n}\f$ in the
            /// [documentation]{@ref InclinationEccentricity_pms1}).
            ///
            ///The outer index is s+__max_e_power+2 and the inner one is 
            ///min(n-1, n+s+1).
            __gamma_minus;

        ///Is the object ready to be used?
        bool __useable;

        ///\brief The inner index in the __alpha/gamma_plus/gamma_minus arrays
        //corresponding to the given term.
        int inner_index( 
                ///If -1 returns the index within __gamma_minus, if 0 returns the
                ///index within __alpha and if 1 returns the index within
                ///__gamma_plus.
                int msign,

                ///The multiplier of the orbital frequency of the desired term.
                int s,

                ///The power of the eccentricity in the desired term.
                int epower) const;

        ///Taylor series approximation of \f$p_{-2,s}(e)\f$.
        double p_m2s(
                ///The eccentricity
                double e, 

                ///The s index.
                int s, 
                
                ///Where to truncate the taylor series.
                unsigned max_e_power,
                
                ///If true the result is differentiated w.r.t. to the
                ///eccentricity.
                bool deriv) const;
        
        ///Taylor series approximation of \f$p_{0,s}(e)\f$.
        double p_0s(
                ///The eccentricity
                double e, 

                ///The s index.
                int s, 
                
                ///Where to truncate the taylor series.
                unsigned max_e_power,
                
                ///If true the result is differentiated w.r.t. to the
                ///eccentricity.
                bool deriv) const;

        ///Taylor series approximation of \f$p_{2,s}(e)\f$.
        double p_p2s(
                ///The eccentricity
                double e, 

                ///The s index.
                int s, 
                
                ///Where to truncate the taylor series.
                unsigned max_e_power,
                
                ///If true the result is differentiated w.r.t. to the
                ///eccentricity.
                bool deriv) const;

    public:
        ///Create an uninitialized object.
        EccentricityExpansionCoefficients() : __useable(false) {}

        ///Reads in tabulated expansion coefficients, making this object useable.
        void read(
                ///The name of the file to read pre-tabulated coefficients from.
                const std::string &tabulated_pms_fname="",

                ///Only reads the file until all coefficients necessary for
                ///expansions of up to this power are read-in. If this is
                ///negative, the full file is read.
                int max_e_power=-1);


        ///Maximum eccentricity power with all necessary coefficients known.
        unsigned max_e_power() const {return __max_e_power;}

        ///Taylor series approximation of \f$p_{m,s}\f$.
        double operator()(
                ///The first index (0 or +-2).
                int m, 

                ///The second index.
                int s,

                ///The value of the eccentricity to use.
                double e,

                ///The maximum eccentricity order to include in the Taylor
                ///series.
                unsigned max_e_power,
                
                ///If true the result is differentiated w.r.t. to the
                ///eccentricity.
                bool deriv) const;
    }; //End EccentricityExpansionCoefficients class.

} //End Evolve namespace.

#endif
