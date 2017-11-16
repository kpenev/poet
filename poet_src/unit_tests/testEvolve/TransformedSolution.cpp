#include "TransformedSolution.h"

namespace Evolve {

    void TransformedSolution::add_transformation(
        const std::vector<const Core::OneArgumentDiffFunction *> &transforms,
        double change_age
    )
    {
        for(unsigned q = 0; q < NUM_REAL_QUANTITIES - 1; ++q)
            __transforms[q].push_back(transforms[q]);

        __change_ages.push_back(change_age);
    }

    const std::vector< const std::list<double> * > &
        TransformedSolution::operator()(
            const std::vector< const std::list<double> * > &solution
        )
    {
        if(__transformed_orbit[AGE]) delete __transformed_orbit[AGE];
        __transformed_orbit[AGE] = solution[AGE];

        for(unsigned q = 0; q < NUM_REAL_QUANTITIES - 1; ++q) {
            std::list<double>::const_iterator
                change_i = __change_ages.begin();
            change_i++;
            std::list<const Core::OneArgumentDiffFunction *>::const_iterator
                transform_i = __transforms[q].begin();
            std::list<double> *transformed_quantity = new std::list<double>();
            for(
                std::list<double>::const_iterator
                    var_i = solution[q]->begin(),
//                    deriv_i = deriv->begin(),
                    t_i = solution[AGE]->begin();
                var_i != solution[q]->end();
                ++var_i,
//                ++deriv_i,
                ++t_i
            ) {
                if(change_i != __change_ages.end() && *t_i >= *change_i) {
                    ++transform_i;
                    ++change_i;
                }
                transformed_quantity->push_back((**transform_i)(*var_i));
/*                const FunctionDerivatives *dvar = (*transform_i)->deriv(*var_i);
                __transformed_orbit[var_type].push_back(dvar->order(0));
                __transformed_deriv[var_type].push_back(dvar->order(1)
                                                        *
                                                        (*deriv_i));*/
            }
            if(__transformed_orbit[q]) delete __transformed_orbit[q]; 
            __transformed_orbit[q] = transformed_quantity;
        }
        return __transformed_orbit;
    }

    TransformedSolution::~TransformedSolution()
    {
        for(unsigned q = 0; q < NUM_REAL_QUANTITIES - 1; ++q)
            if(__transformed_orbit[q]) delete __transformed_orbit[q];
    }

}//End Evolve namespace.
