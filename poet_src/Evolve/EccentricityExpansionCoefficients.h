/**\file
 *
 * \brief Declares a class which provides the
 * [\f$p_{m,s}\f$ coefficients]{@ref InclinationEccentricity_pms1}.
 */


#ifndef __ECCENTRICITY_EXPANSION_COEFFICIENTS_H
#define __ECCENTRICITY_EXPANSION_COEFFICIENTS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Error.h"
#include "../Core/Common.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <sqlite3.h>
//#include <boost/math/special_functions/chebyshev.hpp>

namespace Evolve {

    ///\brief A class which reads-in and provides a convenient interface to the
    /// \f$p_{m,s}\f$ coefficients
    class LIB_PUBLIC EccentricityExpansionCoefficients {
    private:
        ///Has this class been properly linked to a database with interpolation
        ///data?
        bool __useable;

        std::vector< std::vector<double> >
            ///\brief The expansion coefficients for all \f$p_{m,s}\f$.
            ///
            ///Stored in the order m=-2,0,+2 for increasing s (m changes
            ///faster).
            __pms_interp_data;


        ///If you're seeing this, it means I haven't properly sorted
        ///out new documentation stuff or moved new variables/functions
        ///into a better place

        ///\brief The largest s (second) index at which all three \f$p_{m,s}\f$
        ///coefficients are tabulated.
        int __max_s;

        ///The currently defined expansion precision (see prepare())
        double __expansion_precision;

        ///\brief Whether we load the whole database at the start (true) or part
        ///of it as needed (false)
        bool __load_all;

        ///\brief The identifiers of particular \f$p_{m,s}\f$ coefficients
        ///in the database.
        std::vector<int> __db_pms_id;

        ///The smallest value of e for which each \f$p_{m,s}\f$ coefficient has
        ///tabulated values.
        std::vector<double> __min_e;

        ///\brief The number of e values at which each \f$p_{m,s}\f$ coefficient
        ///is tabulated
        std::vector<int> __num_steps;

        ///\brief The maximum eccentricity at which each \f$p_{m,s}\f$
        ///coefficient can be reliably interpolated.
        std::vector<double> __max_e;


        ///\brief The guaranteed interpolation precision for each \f$p_{m,s}\f$ u
        ///coefficient.
        std::vector<double> __interp_precision;

        ///The name of the file contanining the interpolatiod sqlite database.
        std::string __file_name;

        ///\brief The maximum eccentricity at which a given \f$p_{m,s}\f$ can be
        ///ignored.
        ///
        ///The order is the same as __pms_interp_data.
        std::vector<double> __max_ignore_eccentricity;

        std::vector<double> load_coefficient(sqlite3* db,int m,int s);

        ///Return the highest s available for requested precision.
        /// TODO: what is this? It is possible for file to accommodate precision but only for s=0
        void set_max_s(sqlite3* db);

        ///\brief Read the metadata for the available \f$p_{m,s}\f$ coefficients
        ///from the databate.
        void load_metadata(sqlite3* db);

        ///Fill the __max_ignore_eccentricity member per the database.
        void load_max_ignore_eccentricity(sqlite3* db, double precision);

        double load_specific_e(int m,int s,int e_step) const;
        ///\brief The callback SQL function that updates the above values
        ///(__pms_interp_data)
        void get_interp_data(sqlite3* db);

        ///Return a single tabulated value of \f$p_{m,s}\f$.
        double get_specific_e(
            ///The m (first) index of the expansion coefficient to get.
            int m,

            ///The s (second) index of the expansion coefficient to get.
            int s,

            ///The index within the tabulated values to get
            int e_step
        ) const;

        std::vector<double> find_pms_boundary_values(int m,int s,double e) const;
        double return_known_e(int m, int s, double e) const;
        bool check_known_e(int m,int s,double e) const;
        inline int e_to_nearest_step(int m,int s,double e,bool flr) const;
        inline double step_to_e(int m,int s,int step) const;
        int current_largest_s(int m);
        inline int local_index(int m, int s) const;
        //void change_frequency_order(double new_e); // Has not been implemented yet but is new and is intended to be

        std::vector<double>* which_list(int m);
        int* which_order(int m);

    public:
        ///Create an uninitialized object.
        EccentricityExpansionCoefficients() : __useable(false) {}

        ///Reads in tabulated expansion coefficients, making this object useable.
        void prepare(
                ///The name of the file to read pre-tabulated coefficients from.
                const std::string &tabulated_pms_fname,

                ///Only reads the coefficients of the expansion which
                ///most closely achieves the indicated precision, erring on
                ///the side of being more precise.
                double precision,

                ///Should all data from the database be pre-loaded to avoid
                ///query each time a coefficient is evaluated. Setting to true
                ///means we keep a 9 GB database in memory.
                bool pre_load
        );

        ///Return the expansion precision set by the last call to prepare()
        inline double get_expansion_precision() const
        {return __expansion_precision;}

        ///The guaranteed interpolation precision for a given \f$p_{m,s}\f$.
        double interp_precision(int m, int s) const;

        ///\brief Return the smallest s value such that all \f$p_{m,s'}\f$ for
        //\f$s'>s\f$ can be ignored from the expansion without violating the
        ///precision requirements.
        int required_expansion_order(double e, int m) const;

        ///The largest s for which \f$p_{m,s'}\f$ is known.
        int max_expansion_order() const
        {return __max_s;}

        ///\brief Return the range of eccentricities (min, max) over which an
        ///expansion going up to given max s is valid and required.
        ///
        ///For eccentricities below the minimum (first value of returned pair)
        ///at least \f$p_{m,max_s}\f$ can be excluded from the expansion without
        ///violating the precision requirement. For eccentricities at or above
        ///maximum (second value of returned pair) at least \f$p_{m,max_s+1}\f$
        ///must be included in the series in order to satisfy the precision
        ///requirement.
        std::pair<double,double> get_expansion_range(int m, int max_s) const;

        ///\brief Approximate the value of \f$p_{m,s}(e)\f$
        double operator()(
                ///The first index (0 or +-2).
                int m,

                ///The second index.
                int s,

                ///The value of the eccentricity to use.
                double e,

                ///If true the result is differentiated w.r.t. to eccentricity.
                bool deriv
        ) const;
    }; //End EccentricityExpansionCoefficients class.

} //End Evolve namespace.

#endif
