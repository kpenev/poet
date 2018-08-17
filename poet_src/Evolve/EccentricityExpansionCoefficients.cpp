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
            gamma_ind1 = s + __max_e_power + 2;
        double e_pow = std::pow(e, s + 2 * min_n - (deriv ? 1 : 0)),
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
            alpha_ind1 = s + __max_e_power;
        double e_pow = std::pow(e, s + 2 * min_n - (deriv ? 1 : 0)),
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
        int min_n = std::max(-1, -s + 1), gamma_ind1 = s + __max_e_power - 2;
        double e_pow = std::pow(e, s + 2 * min_n - (deriv ? 1 : 0)),
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

    void EccentricityExpansionCoefficients::read(
        const std::string &tabulated_pms_fname,
        int max_e_power
    )
    {
        std::ifstream tabulated_coef(tabulated_pms_fname.c_str());
        if(!tabulated_coef) throw Core::Error::IO(
            "Unable to open eccentricity expansion file: "
            +
            tabulated_pms_fname
            +
            "!"
        );
        tabulated_coef >> __max_e_power;
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
        }

        __alpha.resize(2 * __max_e_power + 1);
        __gamma_plus.resize(2 * __max_e_power + 1);
        __gamma_minus.resize(2 * __max_e_power + 1);

        for(
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
            }
        }
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
