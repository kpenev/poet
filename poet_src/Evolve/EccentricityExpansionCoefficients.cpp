#define BUILDING_LIBRARY
#include "EccentricityExpansionCoefficients.h"

namespace Evolve {

    void EccentricityExpansionCoefficients::get_expansions(sqlite3* db,std::vector<int> id_list)             // bring back pair, store error for a series?
    { // sanity checking? get m and s and fill in proper thing. COUNT (new function, prefill __pms(etc))
        // different statement count first for knowing wherefore many of your friends are there romeo
        bool error_flag = false;
        
        for(int i_a==0; i_a < id_list.size() && !error_flag; i_a++)
        {
            sqlite3_stmt **statement;
            std::string instruc2="SELECT y_value,step_number FROM interpolation_data WHERE p_id = "+std::to_string(id_list[i_a])+" ORDER BY step_number DESC";
            const char *sql = instruc2.c_str();
            if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
            {
                std::vector<double> new_expansion;
                int rc = sqlite3_step(*statement);
                int i_b = sqlite3_column_int(*statement,1);
                new_expansion.resize(i_b+1);
                while(rc==SQLITE_ROW)
                {
                    new_expansion[i_b]=( sqlite3_column_double(*statement,0) );
                    rc=sqlite3_step(*statement);
                    i_b--;
                }
                
                if (rc==SQLITE_DONE) __pms_expansions[i_a] = new_expansion;
                else error_flag=true;
            }
            else error_flag=true;
            
            sqlite3_finalize(*statement);
        }
        
        if(error_flag) throw Core::Error::IO("Unable to retrieve expansion in eccentricity expansion file!");
    }
    
    std::vector<int> EccentricityExpansionCoefficients::identify_expansions(sqlite3* db,int max_s,double precision)
    { // but can keep but have more arguments?
        sqlite3_stmt **statement;
        std::string instruc2="SELECT id,m,s,interp_accuracy FROM interpolations WHERE interp_accuracy <= "+std::to_string(precision)+" ORDER BY s,m,accuracy DESC";
        const char *sql = instruc2.c_str();
        bool error_flag = false;
        
        int m = -2, s = 0, i = 0;
        std::vector<int> expansion_ids;
        
        expansion_ids.resize(3*(max_s+1));
        
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            while(rc==SQLITE_ROW && i+1<=3*(max_s+1))
            {
                if(
                    sqlite3_column_int(*statement,1) == m
                    &&
                    sqlite3_column_int(*statement,2) == s
                ) {
                    expansion_ids[i]=sqlite3_column_int(*statement,0);
                    __max_precision[i]=sqlite3_column_double(*statement,3);
                    i++;
                    if(m==2) {
                        m=-2;
                        s++;
                    }
                    else m+=2;
                }
                rc=sqlite3_step(*statement);
            }
            if (rc!=SQLITE_DONE) error_flag=true;
        }
        else error_flag=true;
        
        sqlite3_finalize(*statement);
        
        if(error_flag)
        {
            throw Core::Error::IO("Unable to search expansions in eccentricity expansion file!");
        }
        
        return expansion_ids;
    }
    ////////////////// if derivative, NaN and fix the code crash later. No pairs.
    int EccentricityExpansionCoefficients::get_max_s(sqlite3* db,double precision)
    { // hey there maybe do a minimum accuracy in python plz notice me when you're deleting comments
        sqlite3_stmt **statement;
        std::string instruc2="WITH meets_precision AS (SELECT m,s FROM m_and_s_to_accuracy WHERE accuracy <= "
                                +std::to_string(precision)+" GROUP BY s,m) "
                                +"SELECT s FROM meets_precision GROUP BY s HAVING COUNT(m)=3";  //aaaa   save actual precision?
        const char *sql = instruc2.c_str();
        bool keep_going=true;
        int compare_s = 0;
        int error_flag = 0;
        
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            while(rc==SQLITE_ROW && keep_going)
            {
                if(sqlite3_column_double(*statement,0)!=compare_s)
                {
                    keep_going=false;
                    if(compare_s==0) error_flag=1;
                }
                else
                {
                    compare_s++;
                    rc=sqlite3_step(*statement);
                }
            }
            if (rc!=SQLITE_DONE && keep_going==true) error_flag=-1;
            if (rc==SQLITE_DONE && keep_going==true && compare_s==0) error_flag=1;
        }
        else error_flag=-1;
        sqlite3_finalize(*statement);
        
        if(error_flag==1)
        {
            std::ostringstream msg;
            msg << "Eccentricity expansion file cannot accommodate requested precision ("
                << precision
                << ") in EccentricityExpansionCoefficients::get_precision()!";
            throw Core::Error::BadFunctionArguments(msg.str());
        }
        if(error_flag==-1)
        {
            std::ostringstream msg;
            msg << "Eccentricity expansion file could not be read in EccentricityExpansionCoefficients::get_precision()!"
            throw Core::Error::IO(msg.str());
        }
        
        return compare_s-1;
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

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        double precision
    )
    {
        sqlite3 *db; // check to make sure this is closed even when there's errors
        int rc;
        std::vector<int> expansion_ids;
        
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
            __max_s = get_max_s(db,precision);
            __pms_expansions.resize(3*(__max_s+1));
            __max_precision.resize(3*(__max_s+1));
            expansion_ids = identify_expansions(db,__max_s,precision);
            get_expansions(db,expansion_ids);
        } catch(...) {
            sqlite3_close(db);
        }
        
        sqlite3_close(db);
        
        __useable = true;
    }

    double EccentricityExpansionCoefficients::max_precision(int m, int s)
    {
        switch(m) {
            case -2 : return __max_precision[(s*3)+0];
            case 0  : return __max_precision[(s*3)+1];
            case 2  : return __max_precision[(s*3)+2];
            default : throw Core::Error::BadFunctionArguments(
                          "Asking EccentricityExpansionCoefficients::max_precision() for p_{m,s} with m other than +-2 and 0"
                      );
        };
    }

    double EccentricityExpansionCoefficients::operator()(
        int m,
        int s,
        double e,
        unsigned max_e_power,
        bool deriv // is this easy to do?     maybe not, make next bit NaN if it wants derivative
    ) const
    {
        if(!__useable)
            throw Core::Error::Runtime(
                "Attempting to evaluate Pms before reading in eccentricity "
                "expansion coefficients!"
            );
        
        if(s > __max_s) throw Core::Error::Runtime("Attempting to evaluate larger s than is available at requested precision!");
        
        if(deriv) return Core::NaN;
        
        std::vector<double> exp_coefs;

        switch(m) {
            case -2 : exp_coefs=__pms_expansions[(s*3)+0];
            case 0  : exp_coefs=__pms_expansions[(s*3)+1];
            case 2  : exp_coefs=__pms_expansions[(s*3)+2];
            default : throw Core::Error::BadFunctionArguments(
                          "Asking for p_{m,s} with m other than +-2 and 0"
                      );
        };
        // double check if I need to do 2*e+1 etc. sort of thing
        return chebyshev_clenshaw_recurrence(exp_coefs.data(),exp_coefs.size(),e); // using the right one? for t0 term? between psipy and boost
    }

} //End Evolve namespace.
