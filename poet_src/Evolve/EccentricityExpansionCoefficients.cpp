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
        
        std::vector<double> pms;
        
        sqlite3_stmt *statement;
        std::string stmnt="SELECT y_value,step_number";
        stmnt.append(" FROM interpolation_data WHERE p_id = ");
        stmnt.append(std::to_string(__db_index[local_index(m,s)]));
        stmnt.append(" ORDER BY step_number DESC");
        const char *sql = stmnt.c_str();
        if(sqlite3_prepare_v2(db,sql,-1,&statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(statement);
            int i_b = sqlite3_column_int(statement,1);
            assert (i_b==__step_num[local_index(m,s)]);
            pms.resize(i_b+1);
            while(rc==SQLITE_ROW)
            {
                pms[i_b]=( sqlite3_column_double(statement,0) );
                i_b--;
                rc=sqlite3_step(statement);
            }
            
            if (rc!=SQLITE_DONE) error_flag=true;
        }
        else error_flag=true;
        sqlite3_finalize(statement);
        
        if(error_flag)
            throw Core::Error::IO(
                "Unable to retrieve expansion id "
                +
                std::to_string(__db_index[local_index(m,s)])
                +
                " in eccentricity expansion file during load_coefficient!"
            );
        
        return pms;
    }
    
    void EccentricityExpansionCoefficients::get_max_s(sqlite3* db)
    {
        sqlite3_stmt *statement;
        std::string stmnt="SELECT MAX(s) FROM interpolations";
        const char *sql = stmnt.c_str();
        int result_s = 0;
        bool error_flag = false;
        
        int codeOne=sqlite3_prepare_v2(db,sql,-1,&statement,NULL);
        if(codeOne==SQLITE_OK)
        {
            int rc = sqlite3_step(statement);
            while(rc==SQLITE_ROW)
            {
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
            msg.append(" EccentricityExpansionCoefficients::get_max_s()!");
            throw Core::Error::IO(msg);
        }
        
        __max_s = result_s;
    }
    
    void EccentricityExpansionCoefficients::load_metadata(sqlite3* db)
    {
        sqlite3_stmt *statement;
        std::string stmnt="SELECT id,m,s,min_interp_e,number_of_steps,";
        stmnt.append("max_checked_e,interp_accuracy FROM interpolations");
        const char *sql = stmnt.c_str();
        int error_flag = 0;
        
        __pms_metadata.resize(3*(__max_s+1));
        __db_index.resize(3*(__max_s+1),-1);
        __min_e.resize(3*(__max_s+1),Core::NaN);
        __step_num.resize(3*(__max_s+1),Core::NaN);
        __max_e.resize(3*(__max_s+1),Core::NaN);
        __accur.resize(3*(__max_s+1),Core::NaN);
        
        if(sqlite3_prepare_v2(db,sql,-1,&statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(statement);
            while(rc==SQLITE_ROW)
            {
                int m=sqlite3_column_int(statement,1);
                int s=sqlite3_column_int(statement,2);
                __db_index[local_index(m,s)]=sqlite3_column_int(statement,0);
                __min_e[local_index(m,s)]=sqlite3_column_double(statement,3);
                __step_num[local_index(m,s)]=sqlite3_column_int(statement,4);
                __max_e[local_index(m,s)]=sqlite3_column_double(statement,5);
                __accur[local_index(m,s)]=sqlite3_column_double(statement,6);
                rc=sqlite3_step(statement);
            }
            if (rc!=SQLITE_DONE) error_flag=-1;
        }
        else error_flag=-1;
        sqlite3_finalize(statement);
        
        if(error_flag==-1)
        {
            std::string msg="Eccentricity expansion file could not be read in ";
            msg.append("EccentricityExpansionCoefficients::load_metadata()!");
            throw Core::Error::IO(msg);
        }
    }
    
    void EccentricityExpansionCoefficients::load_e_switches(
        sqlite3* db,
        double precision
    )
    {
        int em = 2;
        
        std::vector<int> mp2;
        std::vector<int> m0;
        std::vector<int> mm2;
        
        mp2.resize(__max_s+1),Core::NaN);
        m0.resize(__max_s+1),Core::NaN);
        mm2.resize(__max_s+1),Core::NaN);
        __order_switches.resize(__max_s+1),Core::NaN);
        
        while (em >= -2)
        {
            for(int es=0; es <= __max_s; es++)
            {
                sqlite3_stmt *statement;
                std::string stmnt="SELECT MIN(b.step_number) FROM";
                stmnt.append(" interpolations a LEFT JOIN interpolation_data ");
                stmnt.append("b ON a.id = b.p_id WHERE a.m = ");
                stmnt.append(std::to_string(em));
                stmnt.append(" AND a.s = ");
                stmnt.append(std::to_string(es));
                stmnt.append(" AND ABS(b.y_value) >= ");
                stmnt.append(std::to_string(precision/double(es)));
                const char *sql = stmnt.c_str();
                int error_flag = 0;
                
                if(sqlite3_prepare_v2(db,sql,-1,&statement,NULL)==SQLITE_OK)
                {
                    int rc = sqlite3_step(statement);
                    while(rc==SQLITE_ROW)
                    {
                        int li = local_index(em,es);
                        switch(em)
                        {
                            case 2: mp2[li]=sqlite3_column_int(statement,0);
                                    break;
                            case 0: m0[li]=sqlite3_column_int(statement,0);
                                    break;
                            case -2: mm2[li]=sqlite3_column_int(statement,0);
                                    break;
                        }
                        rc=sqlite3_step(statement);
                    }
                    if (rc!=SQLITE_DONE) error_flag=-1;
                }
                else error_flag=-1;
                sqlite3_finalize(statement);
                
                if(error_flag==-1)
                {
                    std::string msg="Eccentricity expansion file could not be ";
                    msg.append("read in EccentricityExpansionCoefficients::");
                    msg.append("load_e_switches()!");
                    throw Core::Error::IO(msg);
                }
            }
            
            em = em - 2;
        }
        
        for(int i = 0; i <= __max_s; i++(
        {
            if ( (mp2[i] <= m0[i]) && (mp2[i] <= mm2[i]) )
                __order_switches[i]=mp2[i]-1;
            else if ( (m0[i] <= mp2[i]) && (m0[i] <= mm2[i]) )
                __order_switches[i]=m0[i]-1;
            else if ( (mm2[i] <= mp2[i]) && (mm2[i] <= m0[i]) )
                __order_switches[i]=mm2[i]-1;
            else
                __order_switches[i]=0;
            
            if(__order_switches[i] < 0)
                __order_switches[i] = 0;
        }
    }
    
    double EccentricityExpansionCoefficients::load_specific_e(
        int m,
        int s,
        int e_step
    ) const
    {
        bool error_flag = false;
        double result;
        
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
            stmnt.append(std::to_string(__db_index[local_index(m,s)]));
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
                    std::to_string(__db_index[local_index(m,s)])
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
    
    void EccentricityExpansionCoefficients::get_expansions(sqlite3* db)
    {
        __pms_expansions.resize(3*(__max_s+1));
        
        for(int s=0; s <= __max_s; s++)
        {
            __pms_expansions[local_index(-2,s)]=load_coefficient(db,-2,s);
            __pms_expansions[local_index(0,s)]=load_coefficient(db,0,s);
            __pms_expansions[local_index(2,s)]=load_coefficient(db,2,s);
        }
    }
    
    double EccentricityExpansionCoefficients::get_specific_e(
        int m,
        int s,
        int e_step
    ) const
    {
        if(!__load_all) return load_specific_e(m,s,e_step);
        else return __pms_expansions[local_index(m,s)][e_step];
    }
    
    std::vector<double> EccentricityExpansionCoefficients::find_pms_boundary_values(
        int m,
        int s,
        double e
    ) const
    {
        std::vector<double> results (4);
        
        int lo_i=e_to_nearest_step(m,s,e,true);
        int hi_i=lo_i+1;
        
        results[0]=step_to_e(m,s,lo_i);
        results[1]=step_to_e(m,s,hi_i);
        results[2]=get_specific_e(m,s,lo_i);
        results[3]=get_specific_e(m,s,hi_i);
        
        return results;
    }
    
    double EccentricityExpansionCoefficients::return_known_e(
        int m,
        int s,
        double e
    ) const
    {
        if(
            ((s==0||e==1.0) && m==0)
            ||
            (m==2&&s==2&&e==0.0)
        ) return 1.0;
        
        if(
            (s==0&&std::abs(m)==2)
            ||
            e==1.0
            ||
            e==0.0
        ) return 0.0;
        
        if( e<__min_e[local_index(m,s)] ) return 0.0;
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
            ||
            (m==2&&s==2&&e==0.0)
        ) return true;
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
                        (double(__step_num[li]) - 1)
                        /
                        ( 1-double(__min_e[li]) );
        if(flr) return int(floor(square_step));
        else return int(ceil(square_step));
    }
    
    inline double EccentricityExpansionCoefficients::step_to_e(
        int m,
        int s,
        int step
    ) const
    {
        int li = local_index(m,s);
        double result = ( double(step) / (double(__step_num[li]) - 1) )
                            *
                            ( 1-double(__min_e[li]) )
                            +
                            double(__min_e[li]);
        return result;
    }
    
    int EccentricityExpansionCoefficients::current_largest_s(int m)
    {
        switch(m) {
            case -2 : return __max_s_for_m2;
            case 0  : return __max_s_for_0;
            case 2  : return __max_s_for_p2;
            default : throw Core::Error::BadFunctionArguments(
                        "Asking EccentricityExpansionCoefficients::current_largest_s() for p_{m,s} with m other than +-2 and 0"
                    );
        };
    }
    
    inline int EccentricityExpansionCoefficients::local_index(int m,int s) const
    {
        switch(m) {
            case -2 : return (s*3)+0;
            case 0  : return (s*3)+1;
            case 2  : return (s*3)+2;
        };
    }

    EccentricityExpansionCoefficients::EccentricityExpansionCoefficients(
        const std::string &tabulated_pms_fname,
        double precision,
        bool load_style
    )
    {
        __useable = false;
        
        sqlite3 *db;
        int rc;
        
        __load_all = load_style;
        
        rc = sqlite3_open(tabulated_pms_fname.c_str(),&db);
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
            get_max_s(db);
            load_metadata(db);
            load_e_switches(db,precision);
            if(__load_all)
            {
                get_expansions(db);
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

    void EccentricityExpansionCoefficients::read(const std::string &gary,int bob)
    {
        gary.c_str();
        bob++;
    }

    double EccentricityExpansionCoefficients::max_precision(int m, int s) const
    {
        if(m!=0&&std::abs(m)!=2)
            throw Core::Error::BadFunctionArguments(
                "Asking EccentricityExpansionCoefficients::max_precision() for p_{m,s} with m other than +-2 and 0"
            );
        return __accur[local_index(m,s)];
    }

    int EccentricityExpansionCoefficients::required_expansion_order(
        double e,
        int m
    ) const
    {
        int largest_s=__max_s;
        for(int s = 0; s <= __max_s; s++)
        {
            // We search the entire thing because it is possible that s
            // activates at >e but s+1 activates at <=e, in which case we would
            // want to include up to s+1
            if( step_to_e(m,s,__order_switches[s])<=e ) largest_s = s;
        }
        largest_s++;
        if(largest_s>__max_s) largest_s--;
        return largest_s;
    }

    double EccentricityExpansionCoefficients::operator()(
        int m,
        int s,
        double e,
        unsigned max_e_power,
        bool deriv
    ) const
    {
        
        if(!__useable)
            throw Core::Error::Runtime(
                "Attempting to evaluate Pms before reading in eccentricity expansion coefficients!"
            );
        
        if(m!=0&&std::abs(m)!=2)
            throw Core::Error::BadFunctionArguments(
                "Asking for p_{m,s} with m other than +-2 and 0"
            );
        if(s>__max_s)
            throw Core::Error::BadFunctionArguments(
                "Attempting to evaluate larger s than is available!"
            );
        if(e>__max_e[local_index(m,s)])
            throw Core::Error::BadFunctionArguments(
                "Attempting to evaluate larger e than is accounted for!"
            );
        
        if(deriv) return Core::NaN;

        if(check_known_e(m,s,e))
            return return_known_e(m,s,e);
        
        std::vector<double> e_and_y_values (4);
        e_and_y_values = find_pms_boundary_values(m,s,e);
        
        double slope = (e_and_y_values[3]-e_and_y_values[2])
                            / (e_and_y_values[1]-e_and_y_values[0]);
        return slope*(e-e_and_y_values[0])+e_and_y_values[2];
        
    }

} //End Evolve namespace.
