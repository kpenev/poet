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
	
	void EccentricityExpansionCoefficients::get_expansion(sqlite3* db,int id)
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
	
	void EccentricityExpansionCoefficients::identify_expansion(sqlite3* db,double precision)
	{ // but can keep but have more arguments?
		sqlite3_stmt **statement;
		const char *sql = "SELECT id,m,s,accuracy FROM m_and_s_to_accuracy";
		bool error_flag = false;
		
		int load_id, loaded_m, loaded_s;
		double loaded_precision;
		int m = 0, s = 0;
		
		if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
		{
			int rc = sqlite3_step(*statement);
			while(rc==SQLITE_ROW)
			{
				load_id=sqlite3_column_int(*statement,0);
				loaded_m=sqlite3_column_int(*statement,1);
				loaded_s=sqlite3_column_int(*statement,2);
				loaded_precision=sqlite3_column_double(*statement,3);
				if(
					loaded_m == m
					&&
					loaded_s == s
					&&
					loaded_precision <= precision
				) {
					get_expansion(db,load_id);
					switch(m) { // Should remember to change order for real table
						case 0: m=2;
						case 2: m=-2;
						case -2: {
							m=0;
							s++;
							break; }
					};
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
			// do I need to also handle sqlite3?
		}
	}
	
// WITH meets_precision AS (
	// SELECT m,s,accuracy
	// FROM m_and_s_to_accuracy
	// WHERE accuracy <= goal_accuracy
	// ORDER BY id),
  // one_of_each AS (
	// SELECT DISTINCT m,s
	// FROM meets_precision
	// GROUP BY s),
  // all_three AS (
	// SELECT COUNT(m),s
	// FROM one_of_each
	// GROUP BY s
	// HAVING COUNT(m)=3),//// make sure it starts at zero
  // do a join AS (
	// SELECT A.s
	// FROM all_three A
	// INNER JOIN all_three B
	  // ON 
	// WHERE )
  // contiguous AS (
	// values where 0 included up to first null
	// FROM all_three )
// SELECT COUNT(m)
// FROM contiguous
	
	std::pair<double, bool> EccentricityExpansionCoefficients::get_precision(sqlite3* db)
	{
		std::pair<double, bool> result(0.0,false);
		
		sqlite3_stmt **statement;  //// nevermind
		const char *sql = "SELECT accuracy FROM m_and_s_to_accuracy WHERE m=0,s=0 ORDER BY id";
		int m = 0, s = 0;
		if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
		{
			int rc = sqlite3_step(*statement);
			while(rc==SQLITE_ROW)
			{
				result.first=sqlite3_column_double(*statement,0);
				rc=sqlite3_step(*statement);
			}
			if (rc!=SQLITE_DONE) result.second=true;
		}
		else result.second=true;
		sqlite3_finalize(*statement);
		return result;
	}

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        double precision
    )
    {
        sqlite3 *db; // check to make sure this is closed even when there's errors
		int rc;
		std::pair<double, bool> max_precision;
		
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
		// try catch
		max_precision = get_precision(db); // maximum s given precision reqs
		if(!max_precision.second) {
			if(max_precision.first != precision) {
				std::ostringstream msg;
				msg << "Eccentricity expansion file '"
					<< tabulated_pms_fname
					<< "' stops at less precision ("
					<< max_precision.first
					<<") than requested ("
					<< precision 
					<< ") in EccentricityExpansionCoefficients::read()!";
				throw Core::Error::BadFunctionArguments(msg.str());
				
				precision = max_precision.first; // why is this here????
			}
			// allocate __pms space b/c we know how many exist now
			identify_expansion(db,precision);
			//errors
			__useable = true;
		}
		
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
