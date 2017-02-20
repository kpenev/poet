#include "CombinedStoppingCondition.h"

namespace Evolve {

    void CombinedStoppingCondition::update_meta_information(
            const StoppingCondition *rhs)
    {
        size_t num_new_subcond=rhs->num_subconditions();
        __num_subconditions+=num_new_subcond;
        __types.reserve(__types.size()+num_new_subcond);
        for(size_t i=0; i<num_new_subcond; i++) __types.push_back(rhs->type(i));
    }

    void CombinedStoppingCondition::add_subcondition_values(
            const StoppingCondition *cond, 
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            size_t &first_index,
            std::valarray<double> &values,
            std::valarray<double> &derivs) const
    {
        std::valarray<double> sub_stop_deriv(cond->num_subconditions()),
            temp_array=(*cond)(evol_mode, orbit, derivatives, sub_stop_deriv);
        for(size_t subcond_ind=0; subcond_ind<temp_array.size();
                subcond_ind++) {
            values[first_index]=temp_array[subcond_ind];
            derivs[first_index]=sub_stop_deriv[subcond_ind];
            ++first_index;
        }
    }

    std::vector<StoppingCondition *>::iterator 
    CombinedStoppingCondition::find_condition(unsigned &index)
    {
        std::vector<StoppingCondition *>::iterator sc_iter=
            __sub_conditions.begin();
        while(index>=(*sc_iter)->num_subconditions()) {
#ifdef DEBUG
            assert(sc_iter!=__sub_conditions.end());
#endif
            index-=(*sc_iter)->num_subconditions();
            ++sc_iter;
        }
        return sc_iter;
    }

    CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
            CombinedStoppingCondition &rhs)
    {
        update_meta_information(&rhs);
        __sub_conditions.insert(__sub_conditions.end(),
                rhs.__sub_conditions.begin(), rhs.__sub_conditions.end());
        if(__delete_subcond) {
#ifdef DEBUG
            assert(rhs.__delete_subcond);
#endif
            const_cast<CombinedStoppingCondition &>(rhs).__delete_subcond=false;
        }
        return *this;
    }

    CombinedStoppingCondition &CombinedStoppingCondition::operator|=(
            StoppingCondition *rhs)
    {
        update_meta_information(rhs);
        __sub_conditions.push_back(rhs);
        return *this;
    }

    std::valarray<double> CombinedStoppingCondition::operator()(
            Core::EvolModeType evol_mode,
            const std::valarray<double> &orbit,
            const std::valarray<double> &derivatives,
            std::valarray<double> &stop_deriv) const
    {
        std::valarray<double> result(__num_subconditions);
        stop_deriv.resize(__num_subconditions, Core::NaN);
        size_t i=0;
        for(std::vector<StoppingCondition *>::const_iterator
                cond=__sub_conditions.begin(); cond!=__sub_conditions.end();
                ++cond)
            add_subcondition_values(*cond, evol_mode, orbit, derivatives, i,
                                    result, stop_deriv);
        return result;
    }

    void CombinedStoppingCondition::reached(short deriv_sign, unsigned index)
    {
        
        std::vector<StoppingCondition *>::iterator sc_iter=find_condition(index);
        (*sc_iter)->reached(deriv_sign, index);
    }

    short CombinedStoppingCondition::expected_crossing_deriv_sign(unsigned index)
    {
        std::vector<StoppingCondition *>::iterator sc_iter=find_condition(index);
        return (*sc_iter)->expected_crossing_deriv_sign(index);
    }

    CombinedStoppingCondition::~CombinedStoppingCondition()
    {
        if(!__delete_subcond) return;
        for(std::vector<StoppingCondition *>::const_iterator
                i=__sub_conditions.begin(); i!=__sub_conditions.end(); ++i)
            delete *i;
    }

} //End Evolve namespace.

