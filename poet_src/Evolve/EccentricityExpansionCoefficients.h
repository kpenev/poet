/**\file
 *
 * \brief Declares a class which provides the
 * [\f$p_{m,s}\f$ coefficients]{@ref InclinationEccentricity_pms1}.
 */


#ifndef __ECCENTRICITY_EXPANSION_COEFFICIENTS_H
#define __ECCENTRICITY_EXPANSION_COEFFICIENTS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Error.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <sqlite3>
#include <boost/math/special_functions/chebyshev.hpp>

namespace Evolve {

    ///\brief A class which reads-in and provides a convenient interface to the 
    /// \f$p_{m,s}\f$ coefficients
    class LIB_PUBLIC EccentricityExpansionCoefficients {
    private:
        ///Maximum eccentricity power with all necessary coefficients known.
        unsigned __max_e_power;

        std::vector< std::vector<double> >
            ///\brief The expansion coefficients for all \f$p_{m,s}\f$.
            ///
            ///Stored in the order m=-2,0,+2 for increasing s.
            __pms_expansions;
        
        //\\\
        std::vector<double> __max_precision;

        ///Is the object ready to be used?
        bool __useable;
        
        ///If you're seeing this, it means I haven't properly sorted
        ///out new documentation stuff or moved new variables/functions
        ///into a better place
        // The callback SQL function that updates the above values
        void get_expansions(sqlite3* db);
        void identify_expansions(sqlite3* db,double precision);
        /// Highest s available for requested precision
        /// It is possible for file to accommodate precision but only for s=0
        int get_max_s(sqlite3* db)
        int __max_s;
        
        bool __load_all; // Whether we load the whole database at the start or part of it as needed
        int __max_s_for_p2;
        int __max_s_for_0;
        int __max_s_for_m2;
        void change_frequency_order(double new_e);
        double __e;
        int current_largest_s(int m);
        std::vector< std::vector<double> >
            ///\brief The expansion coefficients for all \f$p_{m,s}\f$.
            ///
            ///Stored in the order m=-2,0,+2 for increasing s.
            __pms_metadata;
        std::vector<int> __db_index;
        std::vector<double> __min_e;
        std::vector<int> __step_num;
        std::vector<double> __max_e;
        std::vector<double> __accur;
        void load_metadata(sqlite3* db);
        void load_pms_on_demand();
        std::string __file_name;
        std::vector<double> load_coefficient(sqlite3* db,int m,int s)(sqlite3* db,int m, int s);
        std::vector<double> grab_specific_pms(int m, int s);
        std::vector<double> get_pms_boundary_values(int m,int s,double e);
        inline int local_index(int m, int s);

    public:
        ///Create an uninitialized object.
        EccentricityExpansionCoefficients() : __useable(false) {}

        ///Reads in tabulated expansion coefficients, making this object useable.
        void read(
                ///The name of the file to read pre-tabulated coefficients from.
                const std::string &tabulated_pms_fname="",

                ///Only reads the coefficients of the expansion which
                ///most closely achieves the indicated precision, erring on
                ///the side of being more precise.
                double precision=1,
                
                ///Which loading style to use. Setting to true means we keep
                ///a 9 GB database in memory.
                bool load_style=false);


        ///Maximum eccentricity power with all necessary coefficients known.
        unsigned max_e_power() const {return __max_e_power;}
        
        ///Maximum precision for a given expansion.
        double max_precision(int m, int s) const;

        ///\brief Taylor series approximation of \f$p_{m,s}\f$ and the
        ///contribution of the highest power eccentricity terms.
        double operator()(
                ///The first index (0 or +-2).
                int m, 

                ///The second index.
                int s,

                ///The value of the eccentricity to use.
                double e,

                ///Prevision, the maximum eccentricity order to include in the Taylor
                ///series. Currently does nothing.
                unsigned max_e_power,
                
                ///Previously, if true the result was differentiated w.r.t. to the
                ///eccentricity. Currently does nothing.
                bool deriv) const;
    }; //End EccentricityExpansionCoefficients class.

} //End Evolve namespace.

#endif
