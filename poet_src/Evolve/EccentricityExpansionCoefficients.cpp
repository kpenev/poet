#define BUILDING_LIBRARY
#include "EccentricityExpansionCoefficients.h"

namespace Evolve {

    std::vector<double> EccentricityExpansionCoefficients::load_coefficient(
        sqlite3* db,
        int m,
        int s
    )
    {
        bool error_flag = false;
        int ms_index = local_index(m, s);

        std::vector<double> pms;

        sqlite3_stmt *statement;
        std::string stmnt = "SELECT y_value, step_number";
        stmnt.append(" FROM interpolation_data WHERE p_id = ");
        stmnt.append(std::to_string(__db_pms_id[ms_index]));
        stmnt.append(" ORDER BY step_number DESC");
        const char *sql = stmnt.c_str();
        if(sqlite3_prepare_v2(db, sql, -1, &statement, NULL) == SQLITE_OK) {
            int rc = sqlite3_step(statement);
            int i_b = sqlite3_column_int(statement, 1);
            assert (i_b == __num_steps[ms_index] - 1);
            pms.resize(i_b + 1);
            while(rc == SQLITE_ROW) {
                pms[i_b] = sqlite3_column_double(statement, 0);
                --i_b;
                rc = sqlite3_step(statement);
            }

            if (rc != SQLITE_DONE)
                error_flag=true;
        }
        else
            error_flag = true;
        sqlite3_finalize(statement);

        if(error_flag)
            throw Core::Error::IO(
                "Unable to retrieve expansion id "
                +
                std::to_string(__db_pms_id[ms_index])
                +
                " in eccentricity expansion file during load_coefficient!"
            );

        return pms;
    }

    ///TODO: use smallest max s over m=+-2 and m=0
    void EccentricityExpansionCoefficients::set_max_s(sqlite3* db)
    {
        sqlite3_stmt *statement;
        std::string stmnt="SELECT MAX(s) FROM interpolations";
        const char *sql = stmnt.c_str();
        int result_s = 0;
        bool error_flag = false;

        int codeOne=sqlite3_prepare_v2(db,sql,-1,&statement,NULL);
        if(codeOne==SQLITE_OK) {
            int rc = sqlite3_step(statement);
            while(rc==SQLITE_ROW) {
                result_s = sqlite3_column_double(statement,0);
                rc=sqlite3_step(statement);
            }
            if (rc!=SQLITE_DONE) error_flag=true;
        }
        else error_flag=true;
        sqlite3_finalize(statement);

        if(error_flag==true)
        {
            std::string msg ="Eccentricity expansion file could not be read in";
            msg.append(" EccentricityExpansionCoefficients::set_max_s()!");
            throw Core::Error::IO(msg);
        }

        __max_s = result_s;
    }

    void EccentricityExpansionCoefficients::load_metadata(sqlite3* db)
    {
        sqlite3_stmt *statement;
        std::string stmnt = "SELECT id, m, s, min_interp_e, number_of_steps,";
        stmnt.append("max_checked_e, interp_accuracy FROM interpolations");
        const char *sql = stmnt.c_str();
        int error_flag = 0;

        __db_pms_id.resize(3 * (__max_s + 1), -1);
        __min_e.resize(3 * (__max_s + 1), Core::NaN);
        __num_steps.resize(3 * (__max_s + 1));
        __max_e.resize(3 * (__max_s + 1), Core::NaN);
        __interp_precision.resize(3 * (__max_s + 1), Core::NaN);

        if(sqlite3_prepare_v2(db, sql, -1, &statement, NULL) == SQLITE_OK)
        {
            int rc = sqlite3_step(statement);
            while(rc == SQLITE_ROW)
            {
                int m = sqlite3_column_int(statement, 1);
                int s = sqlite3_column_int(statement, 2);
                int ms_index = local_index(m, s);
                __db_pms_id[ms_index] = sqlite3_column_int(statement, 0);
                __min_e[ms_index] = sqlite3_column_double(statement, 3);
                __num_steps[ms_index] = sqlite3_column_int(statement, 4);
                __max_e[ms_index] = sqlite3_column_double(statement, 5);
                __interp_precision[ms_index] = sqlite3_column_double(
                    statement,
                    6
                );
                rc = sqlite3_step(statement);
            }
            if (rc != SQLITE_DONE)
                error_flag = -1;
        } else
            error_flag = -1;
        sqlite3_finalize(statement);

        if(error_flag == -1)
            throw Core::Error::IO(
                "Eccentricity expansion file could not be read in "
                "EccentricityExpansionCoefficients::load_metadata()!"
            );
    }

    void EccentricityExpansionCoefficients::load_max_ignore_eccentricity(
        sqlite3* db,
        double precision
    )
    {
        __max_ignore_eccentricity.resize(3 * (__max_s + 1), Core::NaN);

        for(int m=2; m >= -2; m -= 2)
        {
            for(int s=0; s <= __max_s; ++s)
            {
                int destination_i = local_index(m, s);
                if(s<=2) {
                    __max_ignore_eccentricity[destination_i] = -Core::Inf;
                } else {
                    sqlite3_stmt *statement;
                    std::ostringstream stmnt;
                    stmnt.setf(std::ios_base::scientific);
                    stmnt.precision(16);
                    stmnt << (
                        "SELECT MIN(b.step_number) FROM interpolations a "
                        "LEFT JOIN interpolation_data b ON a.id = b.p_id "
                        "WHERE "
                    )
                        << "a.m = " << m
                        << " AND a.s = " << s
                        << " AND ABS(b.y_value) >= " << precision/double(s);
                    int error_flag = 0;

                    if(
                        sqlite3_prepare_v2(db,
                                           stmnt.str().c_str(),
                                           -1,
                                           &statement,
                                           NULL)
                        ==
                        SQLITE_OK
                    ) {
                        if(sqlite3_step(statement) != SQLITE_ROW)
                            error_flag=-1;
                        else {
                            int max_ignore_step = (
                                sqlite3_column_int(statement, 0)
                                -
                                1
                            );
                            __max_ignore_eccentricity[destination_i] = (
                                max_ignore_step >= 0
                                ? step_to_e(m, s, max_ignore_step)
                                : -Core::Inf
                            );
                            if(sqlite3_step(statement) != SQLITE_DONE)
                                error_flag=-1;
                        }
                    } else
                        error_flag = -1;
                    sqlite3_finalize(statement);

                    if(error_flag == -1)
                        throw Core::Error::IO(
                            "Eccentricity expansion file could not be read in "
                            "EccentricityExpansionCoefficients::"
                            "load_max_ignore_eccentricity()!"
                        );
                }
            }
        }
    }

    double EccentricityExpansionCoefficients::load_specific_e(
        int m,
        int s,
        int e_step
    ) const
    {
        assert(__useable);

        bool error_flag = false;
        double result;
        int ms_index = local_index(m, s);
        sqlite3 *db;

        int rc = sqlite3_open(__file_name.c_str(),&db);
        if(rc!=SQLITE_OK) {
            sqlite3_close(db);
            throw Core::Error::IO(
                "Unable to open eccentricity expansion file: "
                +
                __file_name
                +
                "!"
            );
        }
        try {
            sqlite3_stmt *statement;
            std::string stmnt="SELECT y_value FROM interpolation_data";
            stmnt.append(" WHERE p_id = ");
            stmnt.append(std::to_string(__db_pms_id[ms_index]));
            stmnt.append(" AND step_number=");
            stmnt.append(std::to_string(e_step));
            const char *sql = stmnt.c_str();
            int theCode=sqlite3_prepare_v2(db,sql,-1,&statement,NULL);
            if(theCode==SQLITE_OK)
            {
                int rc = sqlite3_step(statement);
                while(rc==SQLITE_ROW)
                {
                    result=( sqlite3_column_double(statement,0) );
                    rc=sqlite3_step(statement);
                }

                if (rc!=SQLITE_DONE)
                {
                    error_flag=true;
                    std::cout<<"Loop finished without being done.\n";
                }
            }
            else
            {
                 error_flag=true;
                 std::cout<<"prepare was not ok. error code " << std::to_string(theCode) <<  "\n";
                 std::cout<<sqlite3_errmsg(db);
            }
            sqlite3_finalize(statement);

            if(error_flag)
                throw Core::Error::IO(
                    "Unable to retrieve expansion id "
                    +
                    std::to_string(__db_pms_id[ms_index])
                    +
                    " in eccentricity expansion file during load_specific_e!"
                );
        } catch(...) {
            sqlite3_close(db);
            throw;
        }
        sqlite3_close(db);
        return result;
    }

    void EccentricityExpansionCoefficients::get_interp_data(sqlite3* db)
    {
        __pms_interp_data.resize(3*(__max_s+1));
        for(int s=0; s <= __max_s; s++)
            for(int m=-2; m<=2; m+=2)
                __pms_interp_data[local_index(m, s)] = load_coefficient(db,
                                                                        m,
                                                                        s);
    }

    double EccentricityExpansionCoefficients::get_specific_e(
        int m,
        int s,
        int e_step
    ) const
    {
        return (
            __load_all
            ? __pms_interp_data[local_index(m, s)][e_step]
            : load_specific_e(m, s, e_step)
        );
    }

    std::vector<double> EccentricityExpansionCoefficients::find_pms_boundary_values(
        int m,
        int s,
        double e
    ) const
    {
        std::vector<double> results (4);
        int lo_i = e_to_nearest_step(m,s,e,true);
        int hi_i = lo_i+1;

        results[0] = step_to_e(m, s, lo_i);
        results[1] = step_to_e(m, s, hi_i);
        results[2] = get_specific_e(m, s, lo_i);
        results[3] = get_specific_e(m, s, hi_i);
        return results;
    }

    double EccentricityExpansionCoefficients::return_known_e(
        int m,
        int s,
        double e
    ) const
    {
        if(
            ( (s==0 || e==1.0) && m==0 )
            ||
            ( m==2 && s==2 && e==0.0 )
        )
            return 1.0;

        if(
            ( s==0 && std::abs(m)==2 )
            ||
            e==1.0
            ||
            e==0.0
        )
            return 0.0;

        if( e < __min_e[local_index(m,s)] )
            return 0.0;

        assert(false);
    }

    bool EccentricityExpansionCoefficients::check_known_e(
        int m,
        int s,
        double e
    ) const
    {
        if(
            s==0
            ||
            e==1.0
            ||
            e==0.0
            ||
            e<__min_e[local_index(m,s)]
        )
            return true;
        return false;
    }

    inline int EccentricityExpansionCoefficients::e_to_nearest_step(
        int m,
        int s,
        double e,
        bool flr
    ) const
    {
        int li = local_index(m,s);
        double square_step;
        square_step = ( e-double(__min_e[li]) )
                        *
                        (double(__num_steps[li]) - 1)
                        /
                        ( 1-double(__min_e[li]) );
        if(flr) return int( floor(square_step) );
        else return int( ceil(square_step) );
    }

    inline double EccentricityExpansionCoefficients::step_to_e(
        int m,
        int s,
        int step
    ) const
    {
        int li = local_index(m,s);
        double result = ( double(step) / (double(__num_steps[li]) - 1) )
                            *
                            ( 1-double(__min_e[li]) )
                            +
                            double(__min_e[li]);
        return result;
    }

    inline int EccentricityExpansionCoefficients::local_index(int m,
                                                              int s) const
    {
        return s * 3 + (m + 2) / 2;
    }

    void EccentricityExpansionCoefficients::prepare(
        const std::string &tabulated_pms_fname,
        double precision,
        bool pre_load,
        bool disable_precision_fail
    )
    {
        __expansion_precision = precision;
        sqlite3 *db;
        int rc;

        __load_all = pre_load;
        __allow_precision_fail = !disable_precision_fail;

        rc = sqlite3_open(tabulated_pms_fname.c_str(), &db);
        if(rc!=SQLITE_OK) {
            sqlite3_close(db);
            throw Core::Error::IO(
                "Unable to open eccentricity expansion file: "
                +
                tabulated_pms_fname
                +
                "!"
            );
        }
        try {
            set_max_s(db);
            load_metadata(db);
            load_max_ignore_eccentricity(db, precision);
            if(__load_all)
            {
                get_interp_data(db);
            }
            else
            {
                __file_name=tabulated_pms_fname;
            }
        } catch(...) {
            sqlite3_close(db);
            throw;
        }
        sqlite3_close(db);
        __useable = true;
    }

    double EccentricityExpansionCoefficients::interp_precision(int m, int s) const
    {
        if(!__useable)
            throw Core::Error::Runtime(
                "Asking for EccentricityExpansionCoefficients::"
                "interp_precision() before reading interpolation data"
            );
        if(m != 0l && std::abs(m) != 2)
            throw Core::Error::BadFunctionArguments(
                "Asking EccentricityExpansionCoefficients::interp_precision() "
                "for p_{m,s} with m other than +-2 and 0"
            );
        return __interp_precision[local_index(m, s)];
    }

    int EccentricityExpansionCoefficients::required_expansion_order(
        double e,
        int m
    ) const
    {
        if(!__useable)
            throw Core::Error::Runtime(
                "Asking EccentricityExpansionCoefficients::"
                "required_expansion_order() before reading interpolation data"
            );
        int s = __max_s;
        while(__max_ignore_eccentricity[local_index(m, s)] >= e)
            --s;
        if(s == __max_s && __allow_precision_fail) {
            std::ostringstream message;
            message << "Tabulated eccentricity interpolation is unsufficient "
                    << "to achieve the specified precision in tidal potential "
                    << "expansion for e = "
                    << e
                    << ", m = "
                    << m;
            throw Core::Error::Runtime(message.str());
        }
        return s;
    }

    std::pair<double,double>
        EccentricityExpansionCoefficients::get_expansion_range(
            int m,
            int max_s
        ) const
    {
        if(!__useable)
            throw Core::Error::Runtime(
                "Asking EccentricityExpansionCoefficients::"
                "get_expansion_range() before reading interpolation data"
            );
#ifdef VERBOSE_DEBUG
        std::cerr << "Max ignore e(m="
                  << m
                  << ", m'="
                  << max_s
                  << ") = "
                  << __max_ignore_eccentricity[local_index(m, max_s)]
                  << "Max ignore e(m="
                  << m
                  << ", m'="
                  << max_s + 1
                  << ") = "
                  << __max_ignore_eccentricity[local_index(m, max_s + 1)]
                  << std::endl;
#endif

        return std::make_pair(
            __max_ignore_eccentricity[local_index(m, max_s)],
            __max_ignore_eccentricity[local_index(m, max_s + 1)]
        );
    }

    double EccentricityExpansionCoefficients::operator()(
        int m,
        int s,
        double e,
        bool deriv
    ) const
    {
        if(s<0)
            return operator()(-m, -s, e, deriv);
        if(!__useable)
            throw Core::Error::Runtime(
                "Attempting to evaluate p_ms before reading interpolation data"
            );

        if(m != 0 && std::abs(m) != 2)
            throw Core::Error::BadFunctionArguments(
                "Asking for p_{m,s} with m other than +-2 and 0"
            );
        if(s > __max_s)
            throw Core::Error::BadFunctionArguments(
                "Attempting to evaluate larger s than is available!"
            );
        if(e > __max_e[local_index(m, s)])
            throw Core::Error::BadFunctionArguments(
                "Attempting to evaluate larger e than is accounted for!"
            );

        if(check_known_e(m, s, e)) {
            //TODO: handle deriv properly
            return (deriv ? Core::NaN : return_known_e(m, s, e));
        }

        std::vector<double> e_and_y_values(4);
        e_and_y_values = find_pms_boundary_values(m, s, e);

        double slope = (
            (e_and_y_values[3] - e_and_y_values[2])
            /
            (e_and_y_values[1] - e_and_y_values[0])
        );

        if(deriv)
            return slope;
        return slope * (e - e_and_y_values[0]) + e_and_y_values[2];
    }

} //End Evolve namespace.
