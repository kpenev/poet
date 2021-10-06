#define BUILDING_LIBRARY
#include "EccentricityExpansionCoefficients.h"

namespace Evolve {

    std::vector<double> EccentricityExpansionCoefficients::load_coefficient(sqlite3* db,int m,int s)
    {
        bool error_flag = false;
        int id;
        switch(m) {
            case -2 : id = (s*3)+3;
            case 0  : id = (s*3)+1;
            case 2  : id = (s*3)+2;
        };
        
        std::vector<double> pms;
        
        sqlite3_stmt **statement;
        std::string instruc2="SELECT y_value,step_number FROM interpolation_data WHERE p_id = "+std::to_string(id)+" ORDER BY step_number DESC";
        const char *sql = instruc2.c_str();
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            int i_b = sqlite3_column_int(*statement,1);
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
        
        if(error_flag) throw Core::Error::IO("Unable to retrieve expansion id " + std::to_string(id) + " in eccentricity expansion file!");
        
        return pms;
    }
    
    void EccentricityExpansionCoefficients::get_expansions(sqlite3* db)             // bring back pair, store error for a series?
    { // sanity checking? get m and s and fill in proper thing. COUNT (new function, prefill __pms(etc))
        // different statement count first for knowing wherefore many of your friends are there romeo
        __pms_expansions.resize(3*(__max_s+1));
        
        for(int s==0; s < __max_s; s++)
        {
            __pms_expansions[s*3]=load_coefficient(db,0,s)
            __pms_expansions[s*3+1]=load_coefficient(db,2,s)
            __pms_expansions[s*3+2]=load_coefficient(db,-2,s)
        }
    }
    
    ////////////////// if derivative, NaN and fix the code crash later. No pairs.
    void EccentricityExpansionCoefficients::get_max_s(sqlite3* db)
    { // hey there maybe do a minimum accuracy in python plz notice me when you're deleting comments
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
            msg << "Eccentricity expansion file could not be read in EccentricityExpansionCoefficients::get_max_s()!"
            throw Core::Error::IO(msg.str());
        }
        
        __max_s = result_s;
    }
    
    void EccentricityExpansionCoefficients::load_metadata(sqlite3* db)
    {
        sqlite3_stmt **statement;
        std::string instruc2="SELECT id,m,s,min_interp_e,number_of_steps,max_checked_e,interp_accuracy FROM interpolations"
        const char *sql = instruc2.c_str();
        int i = 0;
        int error_flag = 0;
        
        __pms_metadata.resize(3*(__max_s+1));
        
        if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
        {
            int rc = sqlite3_step(*statement);
            while(rc==SQLITE_ROW)
            {
                std::vector<double> data_row (7,0);
                for(int a=0;a<7;a++)
                    data_row[a] = sqlite3_column_double(*statement,a);
                __pms_metadata[i]=data_row;
                i++;
                rc=sqlite3_step(*statement);
            }
            if (rc!=SQLITE_DONE) error_flag=-1;
        }
        else error_flag=-1;
        sqlite3_finalize(*statement);
        
        if(error_flag==-1)
        {
            std::ostringstream msg;
            msg << "Eccentricity expansion file could not be read in EccentricityExpansionCoefficients::load_metadata()!"
            throw Core::Error::IO(msg.str());
        }
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
    
    std::vector<double> EccentricityExpansionCoefficients::grab_specific_pms(int m,int s)
    {
        if(__load_all)
        {
            switch(m) {
                case -2 : return __pms_expansions[(s*3)+2];
                case 0  : return __pms_expansions[(s*3)+0];
                case 2  : return __pms_expansions[(s*3)+1];
            };
        }
        else
        {
            std::vector<double> result;
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
                result = load_coefficient(db,m,s);
            } catch(...) {
                sqlite3_close(db);
            }
            
            sqlite3_close(db);
            
            return result;
        }
    }
    
    std::vector<double> EccentricityExpansionCoefficients::get_pms_boundary_values(int m,int s,double e)
    {
        std::vector<double> results (4);
        std::vector<double> pms = grab_specific_pms(m,s);
        
        return results;
    }

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        double precision,
        bool load_style
    )
    {
        sqlite3 *db; // check to make sure this is closed even when there's errors
        int rc;
        std::vector<int> expansion_ids;
        
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
        switch(m) {
            case -2 : return __pms_metadata[(s*3)+2][6];
            case 0  : return __pms_metadata[(s*3)+0][6];
            case 2  : return __pms_metadata[(s*3)+1][6];
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
        
        if(s > __max_s) throw Core::Error::BadFunctionArguments("Attempting to evaluate larger s than is available!");
        if(m!=0&&std::abs(m)!=2) throw Core::Error::BadFunctionArguments("Asking for p_{m,s} with m other than +-2 and 0");
        
        if(deriv) return Core::NaN;
        
        /* std::vector<double> exp_coefs;

        switch(m) {
            case -2 : exp_coefs=__pms_expansions[(s*3)+0];
            case 0  : exp_coefs=__pms_expansions[(s*3)+1];
            case 2  : exp_coefs=__pms_expansions[(s*3)+2];
            default : throw Core::Error::BadFunctionArguments(
                          "Asking for p_{m,s} with m other than +-2 and 0"
                      );
        };
        // double check if I need to do 2*e+1 etc. sort of thing
        return chebyshev_clenshaw_recurrence(exp_coefs.data(),exp_coefs.size(),e); */ // using the right one? for t0 term? between psipy and boost
        
        if(s==0)
        {
            switch(m) {
                case -2 : return 0.0;
                case 0  : return 1.0;
                case 2  : return 0.0;
            };
        }
        
        if(e==1.0)
        {
            if(m==0) return 1.0;
            else return 0.0;
        }
        
        std::vector<double> e_and_y_values (4);
        e_and_y_values = get_pms_boundary_values(m,s,e);
        double slope = (e_and_y_values[3]-e_and_y_values[2]) / (e_and_y_values[1]-e_and_y_values[0]);
        return slope*(e-e_and_y_values[0])+e_and_y_values[2];
        
    }

} //End Evolve namespace.
