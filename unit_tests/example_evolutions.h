#ifndef __EXAMPLE_EVOLUTIONS_H
#define __EXAMPLE_EVOLUTIONS_H

#include "Common.h"
#include "../OrbitSolver.h"
#include "../YRECIO.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

void output_solution(const OrbitSolver &solver, const StellarSystem &system,
		const std::string &filename);
void calculate_test();

#endif
