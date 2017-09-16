/**\file
 *
 * \brief Outputs the interpolated stellar quantities derived from the YREC
 * tracks.
 *
 * \ingroup UnitTests_group
 */

#ifndef __OUTPUT_YREC_EVOLUTION_H
#define __OUTPUT_YREC_EVOLUTION_H

#include "../StellarEvolution.h"
#include "../Functions.h"
#include "../StellarZone.h"
#include "../YRECIO.h"
#include "../MESAIO.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <assert.h>

///\brief Outputs the given evolution properly formatted with all quantities
///in columns with comments labeling them etc.
void output_evolution(std::ostream &os, const StellarEvolution &evolution,
		double mass);

#endif
