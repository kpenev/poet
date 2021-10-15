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
        
        sqlite3_stmt **statement;
        std::string instruc2="SELECT y_value,step_number"
                +
                " FROM interpolation_data WHERE p_id = "
                +
                std::to_string(__db_index[local_index(m,s)])
                +
                " ORDER BY step_number DESC";
        const char *sql = instruc2.c_str();
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            int i_b = sqlite3_column_int(*statement,1);
            assert (i_b==__step_num[local_index(m,s)]);
            pms.resize(i_b+1);
            while(rc==SQLITE_ROW)
            {
                pms[i_b]=( sqlite3_column_double(*statement,0) );
                i_b--;
                rc=sqlite3_step(*statement);
            }
            
            if (rc!=SQLITE_DONE) error_flag=true;
        }
        else error_flag=true;
        sqlite3_finalize(*statement);
        
        if(error_flag)
            throw Core::Error::IO(
                "Unable to retrieve expansion id "
                +
                std::to_string(__db_index[local_index(m,s)])
                +
                " in eccentricity expansion file!"
            );
        
        return pms;
    }
    
    void EccentricityExpansionCoefficients::get_max_s(sqlite3* db)
    {
        sqlite3_stmt **statement;
        std::string instruc2="SELECT MAX(s) FROM interpolations";
        const char *sql = instruc2.c_str();
        int result_s = 0;
        bool error_flag = false;
        
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            while(rc==SQLITE_ROW)
            {
                result_s = sqlite3_column_double(*statement,0)
                rc=sqlite3_step(*statement);
            }
            if (rc!=SQLITE_DONE) error_flag=true;
        }
        else error_flag=true;
        sqlite3_finalize(*statement);
        
        if(error_flag==true)
        {
            std::ostringstream msg;
            msg << "Eccentricity expansion file could not be read in"
                        +
                        " EccentricityExpansionCoefficients::get_max_s()!"
            throw Core::Error::IO(msg.str());
        }
        
        __max_s = result_s;
    }
    
    void EccentricityExpansionCoefficients::load_metadata(sqlite3* db)
    {
        sqlite3_stmt **statement;
        std::string instruc2="SELECT id,m,s,min_interp_e,number_of_steps,"
                +
                "max_checked_e,interp_accuracy FROM interpolations"
        const char *sql = instruc2.c_str();
        int i = 0;
        int error_flag = 0;
        
        __pms_metadata.resize(3*(__max_s+1));
        __db_index.resize(3*(__max_s+1),-1);
        __min_e.resize(3*(__max_s+1),Core::NaN);
        __step_num.resize(3*(__max_s+1),Core::NaN);
        __max_e.resize(3*(__max_s+1),Core::NaN);
        __accur.resize(3*(__max_s+1),Core::NaN);
        
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            while(rc==SQLITE_ROW)
            {
                int m=sqlite3_column_int(*statement,1);
                int s=sqlite3_column_int(*statement,2);
                __db_index[local_index(m,s)]=sqlite3_column_int(*statement,0);
                __min_e[local_index(m,s)]=sqlite3_column_double(*statement,3);
                __step_num[local_index(m,s)]=sqlite3_column_int(*statement,4);
                __max_e[local_index(m,s)]=sqlite3_column_double(*statement,5);
                __accur[local_index(m,s)]=sqlite3_column_double(*statement,6);
                rc=sqlite3_step(*statement);
            }
            if (rc!=SQLITE_DONE) error_flag=-1;
        }
        else error_flag=-1;
        sqlite3_finalize(*statement);
        
        if(error_flag==-1)
        {
            std::ostringstream msg;
            msg << "Eccentricity expansion file could not be read in "
                        +
                        "EccentricityExpansionCoefficients::load_metadata()!"
            throw Core::Error::IO(msg.str());
        }
    }
    
    double EccentricityExpansionCoefficients::load_specific_e(int m,int s,int e_step)
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
            sqlite3_stmt **statement;
            std::string instruc2="SELECT y_value FROM interpolation_data"
                    +
                    "WHERE p_id = "
                    +
                    std::to_string(__db_index[local_index(m,s)])
                    +
                    ",step_number="
                    +
                    std::to_string(e_step);
            const char *sql = instruc2.c_str();
            if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
            {
                int rc = sqlite3_step(*statement);
                while(rc==SQLITE_ROW)
                {
                    result=( sqlite3_column_double(*statement,0) );
                    rc=sqlite3_step(*statement);
                }
                
                if (rc!=SQLITE_DONE) error_flag=true;
            }
            else error_flag=true;
            sqlite3_finalize(*statement);
            
            if(error_flag)
                throw Core::Error::IO(
                    "Unable to retrieve expansion id "
                    +
                    std::to_string(__db_index[local_index(m,s)])
                    +
                    " in eccentricity expansion file!"
                );
        } catch(...) {
            sqlite3_close(db);
        }
        
        sqlite3_close(db);
        
        return result;
    }
    
    void EccentricityExpansionCoefficients::get_expansions(sqlite3* db)
    {
        __pms_expansions.resize(3*(__max_s+1));
        
        for(int s==0; s < __max_s; s++)
        {
            __pms_expansions[local_index(-2,s)]=load_coefficient(db,-2,s);
            __pms_expansions[local_index(0,s)]=load_coefficient(db,0,s);
            __pms_expansions[local_index(2,s)]=load_coefficient(db,2,s);
        }
    }
    
    double EccentricityExpansionCoefficients::get_specific_e(int m,int s,int e_step)
    {
        if(!__load_all) return load_specific_e(m,s,e_step);
        else return __pms_expansions[local_index(m,s)][e_step];
    }
    
    std::vector<double> EccentricityExpansionCoefficients::find_pms_boundary_values(
        int m,
        int s,
        double e
    )
    {
        std::vector<double> results (4);
        
        int lo_i=e_to_nearest_step(m,s,e,true);
        int hi_i=e_to_nearest_step(m,s,e,false);
        
        results[0]=step_to_e(m,s,lo_i);
        results[1]=step_to_e(m,s,hi_i);
        results[3]=get_specific_e(m,s,lo_i);
        results[4]=get_specific_e(m,s,hi_i);
        
        return results;
    }
    
    double EccentricityExpansionCoefficients::return_known_e(int m,int s,double e)
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
        if(e<__min_e[local_index(m,s)]) return 0.0;
        
        if(e_to_nearest_step(m,s,e,true)==e_to_nearest_step(m,s,e,false))
            return get_specific_e(m,s,e_to_nearest_step(m,s,e,true));
    }
    
    bool EccentricityExpansionCoefficients::check_known_e(int m,int s,double e)
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
            e_to_nearest_step(m,s,e,true)==e_to_nearest_step(m,s,e,false)
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
    )
    {
        double square_step;
        square_step = (e-__min_e[local_index(m,s)])
                        *
                        __step_num[local_index(m,s)]
                        /
                        (1-__min_e[local_index(m,s)]);
        if(flr)
            return int(floor(square_step));
        else
            return int(ceil(square_step));
    }
    
    inline double EccentricityExpansionCoefficients::step_to_e(int m,int s,int step)
    {
        return (step / __step_num[local_index(m,s)])
                    *
                    (1-__min_e[local_index(m,s)])
                    +
                    __min_e[local_index(m,s)];
    }
    
    int EccentricityExpansionCoefficients::current_largest_s(int m)
    {
        switch(m) {
            case -2 : return __max_s_for_m2;
            case 0  : return __max_s_for_0;
            case 2  : return __max_s_for_p2;
            default : throw Core::Error::BadFunctionArguments(
                        "Asking "
                        +
                        "EccentricityExpansionCoefficients::current_largest_s()"
                        +
                        " for p_{m,s} with m other than +-2 and 0"
                    );
        };
    }
    
    inline int EccentricityExpansionCoefficients::local_index(int m, int s)
    {
        switch(m) {
            case -2 : return (s*3)+0;
            case 0  : return (s*3)+1;
            case 2  : return (s*3)+2;
        };
    }

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        double precision,
        bool load_style
    )
    {
        sqlite3 *db;
        int rc;
        
        __load_all = load_style;
        
        rc = sqlite3_open(tabulated_pms_fname.c_str(),&db)
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
            if(__load_all)
            {
                __file_name=tabulated_pms_fname;
                get_expansions(db);
            }
        } catch(...) {
            sqlite3_close(db);
        }
        sqlite3_close(db);
        
        __useable = true;
    }

    double EccentricityExpansionCoefficients::max_precision(int m, int s)
    {
        if(m!=0&&std::abs(m)!=2)
            throw Core::Error::BadFunctionArguments(
                "Asking "
                +
                "EccentricityExpansionCoefficients::max_precision()"
                +
                " for p_{m,s} with m other than +-2 and 0"
            );
        return __accur[local_index(m,s)];
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
                "Attempting to evaluate Pms before reading in eccentricity "
                +
                "expansion coefficients!"
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
        
        if(check_known_e(int m,int s,double e))
            return return_known_e(int m,int s,double e);
        
        std::vector<double> e_and_y_values (4);
        e_and_y_values = find_pms_boundary_values(m,s,e);
        double slope = (e_and_y_values[3]-e_and_y_values[2])
                            / (e_and_y_values[1]-e_and_y_values[0]);
        return slope*(e-e_and_y_values[0])+e_and_y_values[2];
        
    }

} //End Evolve namespace.
