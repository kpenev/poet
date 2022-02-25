#define BUILDING_LIBRARY
#include "EvolvingStellarCore.h"

namespace Star {

    double EvolvingStellarCore::positive_definite_quantity(
        size_t quantity,
        double age,
        unsigned deriv_order,
        double min_value
    ) const
    {
        if(deriv_order == 0) {
            return std::max(any_age_quantity(quantity, age, 0), min_value);
        } else {
            if(any_age_quantity(quantity, age, 0) < min_value)
                return 0.0;
            else
                return any_age_quantity(quantity, age, deriv_order);
        }
    }

    double EvolvingStellarCore::current_positive_definite_quantity(
        size_t quantity,
        unsigned deriv_order,
        double min_value
    ) const
    {
        if(deriv_order == 0) {
            return std::max(current_age_quantity(quantity, 0), min_value);
        } else {
            if(current_age_quantity(quantity, 0) < min_value)
                return 0.0;
            else
                return current_age_quantity(quantity, deriv_order);
        }

    }
}//End Star namespace.
