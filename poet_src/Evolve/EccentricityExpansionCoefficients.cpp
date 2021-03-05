#define BUILDING_LIBRARY
#include "EccentricityExpansionCoefficients.h"

namespace Evolve {

    int EccentricityExpansionCoefficients::inner_index(int msign,
                                                       int s,
                                                       int epower) const
    {
        assert(std::abs(msign) < 2);

        return (epower - s + 2 * std::min(msign, s - msign)) / 2;
    }

    std::pair<double, double> EccentricityExpansionCoefficients::p_m2s( // boost
        double e,
        int s, 
        unsigned max_e_power,
        bool deriv
    ) const
    {
        std::pair<double, double> result(0.0, 0.0);

        double e2 = std::pow(e, 2);
        int min_n = std::max(1, -s - 1),
            gamma_ind1 = s + __max_e_power + 2,
            e_pow_ind = s + 2 * min_n - (deriv ? 1 : 0);
        assert(e_pow_ind + (deriv ? 1 : 0) >= 0);
        double e_pow = (e_pow_ind < 0 ? 0.0 : std::pow(e, e_pow_ind)),
               coef = (deriv ? s + 2 * min_n : 1);
        for(
            int gamma_ind2 = 0;
            gamma_ind2 <= inner_index(-1, s, max_e_power);
            ++gamma_ind2
        ) {
            result.second = (coef
                             *
                             __gamma_minus[gamma_ind1][gamma_ind2]
                             *
                             e_pow);
            result.first += result.second;
            e_pow *= e2;
            if(deriv) coef += 2;
        }
        return result;
    }

    std::pair<double, double> EccentricityExpansionCoefficients::p_0s(
        double e,
        int s, 
        unsigned max_e_power,
        bool deriv
    ) const
    {
        std::pair<double, double> result(0.0, 0.0);

        double e2 = std::pow(e, 2);
        int min_n = std::max(0, -s),
            alpha_ind1 = s + __max_e_power,
            e_pow_ind = s + 2 * min_n - (deriv ? 1 : 0);
        assert(e_pow_ind + (deriv ? 1 : 0) >= 0);
        double e_pow = (e_pow_ind < 0 ? 0.0 : std::pow(e, e_pow_ind)),
               coef = (deriv ? s + 2 * min_n : 1);
        for(
            int alpha_ind2 = 0;
            alpha_ind2 <= inner_index(0, s, max_e_power);
            ++alpha_ind2
        ) {
            result.second = coef * __alpha[alpha_ind1][alpha_ind2] * e_pow;
            result.first += result.second;
            e_pow *= e2;
            if(deriv) coef += 2;
        }
        return result;
    }

    std::pair<double, double> EccentricityExpansionCoefficients::p_p2s(
        double e,
        int s,
        unsigned max_e_power,
        bool deriv
    ) const
    {
        std::pair<double, double> result(0.0, 0.0);

        double e2 = std::pow(e, 2);
        int min_n = std::max(-1, -s + 1),
            gamma_ind1 = s + __max_e_power - 2,
            e_pow_ind = s + 2 * min_n - (deriv ? 1 : 0);
        assert(e_pow_ind + (deriv ? 1 : 0) >= 0);
        double e_pow = (e_pow_ind < 0 ? 0.0 : std::pow(e, e_pow_ind)),
               coef = (deriv ? s + 2 * min_n : 1);
        for(
            int gamma_ind2 = 0;
            gamma_ind2 <= inner_index(1, s, max_e_power);
            ++gamma_ind2
        ) {
            result.second = coef * __gamma_plus[gamma_ind1][gamma_ind2] * e_pow;
            result.first += result.second;
            e_pow *= e2;
            if(deriv) coef += 2;
        }
        return result;
    }
	
	void EccentricityExpansionCoefficients::get_expansion(sqlite3* db,std::vector<int> id)
	{ // sanity checking? get m and s and fill in proper thing. COUNT (new function, prefill __pms(etc))
		sqlite3_stmt **statement;
		std::string instruc2="SELECT coefficient_value FROM cheb_expansion_coeffs WHERE id = "+std::to_string(id)+" ORDER BY place_in_expansion";
		const char *sql = instruc2.c_str();
		// different statement count first for knowing wherefore many of your friends are there romeo
		bool error_flag = false;
		
		std::vector<double> new_expansion;
		// iterator bob = specific part of __pms for happy typing time (if I want, I have the power)
		if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
		{
			
			int rc = sqlite3_step(*statement);
			while(rc==SQLITE_ROW)
			{
				new_expansion[i]=( sqlite3_column_double(*statement,0) ); // just do __pms b/c we'll know how big it is
				rc=sqlite3_step(*statement);
			}
			
			// check error codes
			if (rc==SQLITE_DONE) __pms_expansions[3*sfs/2^17]=.push_back( new_expansion );
			else error_flag=true;
		}
		else error_flag=true;
		
		sqlite3_finalize(*statement);
		
		if(error_flag)
		{
			throw Core::Error::IO("Unable to retrieve expansion in eccentricity expansion file!");
			// do I need to also handle sqlite3?
		}
	}
	
	std::vector<int> EccentricityExpansionCoefficients::identify_expansions(sqlite3* db,int max_s,double precision)
	{ // but can keep but have more arguments?
		sqlite3_stmt **statement;
		//const char *sql = "SELECT id,m,s,accuracy FROM m_and_s_to_accuracy";
		std::string instruc2="SELECT id,m,s,accuracy FROM m_and_s_to_accuracy WHERE accuracy <= "+std::to_string(precision)+" ORDER BY s,m,accuracy DESC";
		const char *sql = instruc2.c_str();
		bool error_flag = false;
		
		int load_id, loaded_m, loaded_s;
		double loaded_precision;
		int m = -2, s = 0, i = 0;
		std::vector<int> expansion_ids;
		
		expansion_ids.resize(3*(max_s+1));
		
		if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
		{
			int rc = sqlite3_step(*statement);
			while(rc==SQLITE_ROW && i+1<=3*(max_s+1))
			{
				loaded_m=sqlite3_column_int(*statement,1);
				loaded_s=sqlite3_column_int(*statement,2);
				if(
					loaded_m == m
					&&
					loaded_s == s
				) {
					expansion_ids[i]=sqlite3_column_int(*statement,0);
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
	
// WITH meets_precision AS (
	// SELECT id,m,s,MAX(accuracy) as precision
	// FROM m_and_s_to_accuracy
	// GROUP BY m,s
	// HAVING accuracy <= goal_accuracy), 
// SELECT COUNT(m),s
// FROM meets_precision
// GROUP BY s
// HAVING COUNT(m)=3

//////////////////////////// eyyyy s has three ms, then c++
//// make sure it starts at zero

// WITH meets_precision AS (SELECT id,m,s,MAX(accuracy) AS precision FROM m_and_s_to_accuracy WHERE accuracy <= .0066 GROUP BY s,m)
// SELECT COUNT(m),s
// FROM meets_precision
// GROUP BY s
// HAVING COUNT(m)=3;
	
	int EccentricityExpansionCoefficients::get_max_s(sqlite3* db,double precision)
	{ // hey there maybe do a minimum accuracy in python plz notice me when you're deleting comments
		sqlite3_stmt **statement;
		std::string instruc2="WITH meets_precision AS (SELECT m,s FROM m_and_s_to_accuracy WHERE accuracy <= "+std::to_string(precision)+" GROUP BY s,m) ";
		instruc2 += "SELECT s FROM meets_precision GROUP BY s HAVING COUNT(m)=3";
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

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        double precision
    )
    {
        sqlite3 *db; // check to make sure this is closed even when there's errors
		int rc;
		int max_s;
		
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
			max_s = get_max_s(db,precision); // maximum s given precision reqs
		} catch(...) {
			sqlite3_close(db);
		}
		// allocate __pms space b/c we know how many exist now
		__pms_expansions.resize(3*(s+1));
		
		identify_expansion(db,precision);
		//errors
		__useable = true;
		
		sqlite3_close(db);
		
		if(max_precision.second) throw Core::Error::IO("Unable to verify precision in eccentricity expansion file!");
		else {
			if (error finding an expansion) //throw error
			else if (error getting an expansion) //throw error
			else // useable now
		}
    }

    std::pair<double, double> EccentricityExpansionCoefficients::operator()(
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
                "expansion coefficients!"
            );

        std::pair<double, double> zero(0.0, 0.0);

        if( // throw exception here if m or s not available b/c of us doing clever SQL stuff and stopping before a set of -2,0,2 that has a gap due to precision not being enough etc.
            s < -static_cast<int>(max_e_power) + m
            ||
            s > static_cast<int>(max_e_power) + m
        )
            return zero;

        switch(m) {
            case -2 : return (s == 0 ? zero : p_m2s(e, s, max_e_power, deriv));
            case 0  : return p_0s(e, s, max_e_power, deriv);
            case 2  : return (s == 0 ? zero : p_p2s(e, s, max_e_power, deriv));
            default : throw Core::Error::BadFunctionArguments(
                          "Asking for p_{m,s} with m other than +-2 and 0"
                      );
        };
    }

} //End Evolve namespace.
