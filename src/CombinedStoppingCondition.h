#ifndef __COMBINED_STOPPING_CONDITION_H
#define __COMBINED_STOPPING_CONDITION_H

#include "StoppingCondition.h"
#include <vector>

///\brief A class combining the the outputs of multiple stopping conditions.
///
///\ingroup OrbitSolver_group
class CombinedStoppingCondition : public StoppingCondition {
private:
	///The conditions that are combined.
	std::vector<StoppingCondition *> __sub_conditions;

	///\brief The number of subconditinos included, allowing for
	///subconditions with multiple entries.
	unsigned __num_subconditions;

	///Whether to delete the sub-conditions then *this is destroyed.
	bool __delete_subcond;

	///\brief The types of the subconditinos, including subconditions of
	///subconditions.
	std::vector<StoppingConditionType> __types;

	///Updates all private methods except the sub_condition lists for adding
	///rhs.
	void update_meta_information(const StoppingCondition *rhs);

	///Adds the values of the given sub-condition to the given array.
	void add_subcondition_values(
			///The stopping condition to process.
			const StoppingCondition *cond,

			///See operator()
			EvolModeType evol_mode,

			///See operator()
			const std::valarray<double> &orbit,
			
			///See operator()
			const std::valarray<double> &derivatives,

			///The index in the values and derivs arrays where to start
			///inserting. On output it is updated to the position immediately
			///after the last insert.
			size_t &first_index,

			///An array to overwrite with the sub-condition values.
			std::valarray<double> &values,

			///An array to overwrite with the sub-condition derivatives.
			std::valarray<double> &derivs) const;
public:
	///Create an empty stopping condition (identical to NoStopCondition).
	CombinedStoppingCondition() :
		__sub_conditions(), __num_subconditions(0), __delete_subcond(true) {}

	///Adds the conditions in RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(CombinedStoppingCondition &rhs);

	///Adds RHS to the conditions of *this.
	CombinedStoppingCondition &operator|=(StoppingCondition *rhs);

	///\brief Returns the values of all stopping sub_conditions, not
	///necessarily in the order added.
	///
	///See StoppingCondition::operator()() for a description of the
	///arguments.
	std::valarray<double> operator()(
			EvolModeType evol_mode,
			const std::valarray<double> &orbit,
			const std::valarray<double> &derivatives,
			std::valarray<double> &stop_deriv) const;

	///Disables the destruction of the subconditions when *this is destroyed.
	void no_delete_subcond() {__delete_subcond=false;}

	virtual size_t num_subconditions() const {return __num_subconditions;}

	StoppingConditionType type(unsigned index=0) const
	{return __types[index];}

	///See StoppingCondition::reached().
	virtual void reached(short deriv_sign, unsigned index=0);

	///See StoppingCondition::expected_crossing_deriv_sign().
	short expected_crossing_deriv_sign(unsigned index=0);

	///\brief Deletes all subconditions, unless no_delete_subcond has been
	///previously called.
	~CombinedStoppingCondition();
};

#endif
