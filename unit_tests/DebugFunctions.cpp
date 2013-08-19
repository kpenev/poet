#include "DebugFunctions.h"

std::string stop_info_to_str(const StopInformation &stop)
{
	std::ostringstream result;
	result << stop << std::endl;
	return result.str();
}
