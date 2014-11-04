/**\file
 * 
 * \brief The implementations of the various stopping condition methods.
 *
 * \ingroup OrbitSolver_group
 */

#include "StoppingConditions.h"

std::ostream &operator<<(std::ostream &os,
		const StoppingConditionType &stop_cond_type)
{
	switch(stop_cond_type) {
		case NO_STOP: os << "NO_STOP"; break;
		case SYNCHRONIZED: os << "SYNCHRONIZED"; break;
		case BREAK_LOCK: os << "BREAK_LOCK"; break;
		case PLANET_DEATH: os << "PLANET_DEATH"; break;
		case WIND_SATURATION: os << "WIND_SATURATION"; break;
		case EXTERNAL: os << "EXTERNAL";
	}
	return os;
}
