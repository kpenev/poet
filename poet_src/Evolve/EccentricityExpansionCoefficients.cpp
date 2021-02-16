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

    std::pair<double, double> EccentricityExpansionCoefficients::p_m2s(
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
	// while learning_func() etc.
	// or not, might as well put logic in the function itself
	bool learning_func(sqlite3* db)
	{
		//
		sqlite3_stmt **statement;
		const char *sql = "specific to the function so why make it an argument?";
		int value;
		
		if(sqlite3_prepare_v2(db,sql,-1,statement,NULL)==SQLITE_OK)
		{
			
			int result = sqlite3_step(*statement);
			while(result==SQLITE_ROW)
			{
				//column
				value = sqlite3_column_(int or double)(*statement,the column we want (int))
				result=sqlite3_step(*statement);
			}
			
			// check error codes
			// I could possibly make another function that just handles all of these
			// and closes out any memory handling stuff if need be
			if (result==SQLITE_DONE)
				return false;
			if(sqlite3_finalize(*statement)!=SQLITE_OK)
				// do an error
			
			return value;
		}
	}

	void EccentricityExpansionCoefficients::get_line(void *data,int numberOfColumns,char **fieldsInRow,char **columnNames)
	{
		// I assume I successfully managed to make this go in order of increasing i, so
		__last_line = std::strtod(*fieldsInRow); // double to int though
	}
	
	void EccentricityExpansionCoefficients::identify_expansion(void *data,int numberOfColumns,char **fieldsInRow,char **columnNames)
	{
		//
		__loaded_m; //1
		__loaded_s;//2
		__loaded_precision=std::strtod(*fieldsInRow[3]); //3
	}
	
	void EccentricityExpansionCoefficients::get_expansion(void *data,int numberOfColumns,char **fieldsInRow,char **columnNames)
	{
		//
	}

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        double precision
    )
    {
        sqlite3 *db;
		int rc;
		
		__last_line = 0; // Maybe this should go in the constructor? Or maybe it shouldn't be available to the entire class in the first place
		
		rc = sqlite3_open(tabulated_pms_fname.c_str(),&db)
		
		/*std::ifstream tabulated_coef(tabulated_pms_fname.c_str());*/
        if(rc) throw Core::Error::IO(
            "Unable to open eccentricity expansion file: "
            +
            tabulated_pms_fname
            +
            "!"
        );
        /* tabulated_coef >> __max_e_power;
        if(max_e_power >= 0) {
            if(__max_e_power<static_cast<unsigned>(max_e_power)) {
                std::ostringstream msg;
                msg << "Eccentricity expansion file '"
                    << tabulated_pms_fname
                    << "' stops at lower eccentricity power ("
                    << __max_e_power
                    <<") than requested ("
                    << max_e_power 
                    << ") in EccentricityExpansionCoefficients::read()!";
                throw Core::Error::BadFunctionArguments(msg.str());
            }
            __max_e_power = max_e_power;
        } */

        __alpha.resize(2 * __max_e_power + 1);
        __gamma_plus.resize(2 * __max_e_power + 1);
        __gamma_minus.resize(2 * __max_e_power + 1);

		int m = 0;
		int s = 0;
		std::string poll_last_line ("SELECT id FROM m_and_s_to_accuracy ORDER BY id");
		std::string poll_tab1 ("SELECT m,s,accuracy FROM m_and_s_to_accuracy WHERE id = ");
		std::string poll_tab2_a ("SELECT coefficient_value FROM cheb_expansion_coeffs WHERE id = ");
		std::string poll_tab2_b (" ORDER BY place_in_expansion");
		rc = sqlite3_exec(db,poll_last_line.c_str(),get_line,data,error_message);
		if(rc!=SQLITE_OK) {
			throw Core::Error::IO("Unable to find number of lines in eccentricity expansion file: "+tabulated_pms_fname+"!");
			// do I need to also handle sqlite3?
		}
		for(
			// all lines in database 1
			i = 0;
			i <= __last_line;
			i++
		) {
			// Read a line
			std::string instruc1 (poll_tab1);
			instruc1+=std::to_string(i);
			rc = sqlite3_exec(db,instruc1.c_str(),identify_expansion,data,error_message);
			// check rc for issues, call errors if needed
			if(rc!=SQLITE_OK){
				throw Core::Error::IO("Unable to search expansions in eccentricity expansion file: "+tabulated_pms_fname+"!");
				// do I need to also handle sqlite3?
			}
			// If m is correct and s is correct
			if(__loaded_m == m && __loaded_s == s) {
				// If precision meets or exceeds requirement
				if(__loaded_precision<=precision) {
					// Grab data from database 2
					std::string instruc2 (poll_tab2_a);
					instruc2+=std::to_string(i);
					instruc2+=poll_tab2_b;
					rc = sqlite3_exec(db,instruc2.c_str(),get_expansion,data,error_message);
					// check rc for issues, call errors if needed
					if(rc!=SQLITE_OK){
						throw Core::Error::IO("Unable to find number of lines in eccentricity expansion file: "+tabulated_pms_fname+"!");
						// do I need to also handle sqlite3?
					}
					__pms_expansions.resize(by one);
					__pms_expansions[whichever]=the_new_data;
					// Depending on value of m
						// Advance m or (s and m)
					switch(m) {
						case 0:
							m=2;
							break;
						case 2:
							m=-2;
							break;
						case -2:
							m=0;
							s++;
							break;
					}
				}
			}
		/*for(
            int epower = 0;
            epower <= static_cast<int>(__max_e_power);
            ++epower
        ) {
            for(int s = -epower - 2; s <= epower-2; s += 2) {
                if(s) {
                    assert(s + static_cast<int>(__max_e_power) + 2 >= 0);
                    assert(s + __max_e_power + 2 < __gamma_minus.size());
                    std::vector<double> 
                        &destination = __gamma_minus[s + __max_e_power + 2];
                    if(destination.size() == 0)
                        destination.resize(inner_index(-1, s, __max_e_power)
                                           +
                                           1);
                    assert(inner_index(-1, s, epower) >= 0);
                    assert(inner_index(-1, s, epower)
                            <static_cast<int>(destination.size()));
                    tabulated_coef >> destination[inner_index(-1,
                                                              s,
                                                              epower)];
                }
            }
            for(int s = -epower; s <= epower; s += 2) {
                assert(s + static_cast<int>(__max_e_power) >= 0);
                assert(s + __max_e_power < __alpha.size());
                std::vector<double> &destination=__alpha[s + __max_e_power];
                if(destination.size() == 0) 
                    destination.resize(inner_index(0, s, __max_e_power) + 1);
                assert(inner_index(0, s, epower) >= 0);
                assert(inner_index(0, s, epower)
                        <static_cast<int>(destination.size()));
                tabulated_coef >> destination[inner_index(0, s, epower)];
            }
            for(int s = -epower + 2; s <= epower + 2; s += 2) {
                if(s) {
                    assert(s + static_cast<int>(__max_e_power) - 2 >= 0);
                    assert(s + __max_e_power - 2 < __gamma_plus.size());
                    std::vector<double> 
                        &destination=__gamma_plus[s + __max_e_power - 2];
                    if(destination.size() == 0)
                        destination.resize(inner_index(1, s, __max_e_power)
                                           +
                                           1);
                    assert(inner_index(1, s, epower) >= 0);
                    assert(inner_index(1, s, epower)
                            <static_cast<int>(destination.size()));
                    tabulated_coef >> destination[inner_index(1, s, epower)];
                }
            }*/
        }
		sqlite3_close(db);
        __useable = true;
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

        if(
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
