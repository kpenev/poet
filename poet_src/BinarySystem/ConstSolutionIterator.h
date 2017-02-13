/**\file
 *
 * \brief Defines the SolutionIterator class.
 *
 */

#ifndef __SOLUTION_ITERATOR_H
#define __SOLUTION_ITERATOR_H

#include "BinarySystem.h"
#include "../IO/IOColumns.h"
#include "EvolvingStar.h"
#include "OrbitSolver.h"
#include "DissipatingZone.h"

///\brief Iterates over the tabulated solution after an evolution
///calculation.
class ConstSolutionIterator {
private:
	///\brief Iterators to the tabulated real values quantities.
	///
	///The order is defined by OutCol::OutputColumns:
	std::vector< std::list<double>::const_iterator > __real_iterators;

	///Location for undefined iterators to point to.
	static std::list<double> __placeholder_list;

	///All undefined entries in __real_iterators are set to this.
	static std::list<double>::const_iterator __placeholder_iterator;

	///Iterator over the tabulated evolution mode.
	std::list<Core::EvolModeType>::const_iterator __mode;

	///Iterator over the tabulated wind saturation state.
	std::list<bool>::const_iterator __wind_saturation;

	///One past the last tabulated age.
	std::list<double>::const_iterator __last_age;

	///The mass of the star.
	double __mstar,

		   ///The mass of the planet.
		   __mplanet;

	///\brief Total angular momentum vector of the star in the reference
	///frame of the surface zone.
	Eigen::Vector3d __stellar_angmom;

	///The star in the system.
	const InterpolatedEvolutionStar &__star;

	///\brief Creates lists for non-orbital quantities and sets the
	///corresponding iterators.
	void create_missing_lists(
			///The ages at which to tabulate quantities.
			const std::list<double> &tabulation_ages
	);

	///\brief Handles the case when no evolution was actually calculated.
	void fix_no_evolution(
			///The starting age if no orbit was calculated. Ignored if solver
			///contains an orbit.
			double start_age,

			///The starting age if no orbit was calculated. Ignored if solver
			///contains an orbit.
			double end_age,

			///The time step if no orbit was calculated. Ignored if solver 
			///contains an orbit.
			double timestep,

			///A list of ages for which an output line must be written.
			///Ignored if solver contains an orbit
			const std::list<double> &required_ages=std::list<double>()
	);
public:
	///Start iterating over a solution.
	ConstSolutionIterator(
			///The solver which calculated the evolution.
			const OrbitSolver &solver,

			///The binary system that was evolved.
			const BinarySystem &system,

			///The star that was evolved need modification privileges only if
			///no evolution was calculated.
			const InterpolatedEvolutionStar &star,

			///The starting age if no orbit was calculated. Ignored if solver
			///contains an orbit.
			double start_age,

			///The starting age if no orbit was calculated. Ignored if solver
			///contains an orbit.
			double end_age,

			///The time step if no orbit was calculated. Ignored if solver 
			///contains an orbit.
			double timestep,

			///A list of ages for which an output line must be written.
			///Ignored if solver contains an orbit
			const std::list<double> &required_ages=std::list<double>()
	);

	///Move to the next tabulated point.
	const ConstSolutionIterator &operator++();

	///Did we move past the last tabulated evolution point.
	operator bool() {return __real_iterators[OutCol::AGE]!=__last_age;}

	///Get a real value quantity for the current step.
	double real_quantity(OutCol::OutputColumns quantity);

	///Get the evolution mode for the current step.
    Core::EvolModeType evolution_mode() {return *__mode;}

	///Was the wind saturated for the current step?
	bool wind_saturation() {return *__wind_saturation;}
};

#endif
