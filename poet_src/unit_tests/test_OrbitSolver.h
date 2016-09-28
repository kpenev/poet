/**\file
 * 
 * \brief Declares the test suite that exercises the OrbitSolver class and
 * some other clasess necessary to accomplish this.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_ORBIT_SOLVER_H
#define __TEST_ORBIT_SOLVER_H

#include "Common.h"
#include "../OrbitSolver.h"
#include "DebugFunctions.h"
#include <iostream>
#include <vector>

/**\brief A class that can be passed to the solution testing function instead
 * of the solver that transforms the solutions before testing.
 *
 * The tranformed orbit is obtained by applying a list of
 * OneArgumentDiffFunction objects to each of the tabulated variables. Each
 * transformation function is used only in a pre-specified age range. 
 *
 * \ingroup UnitTests_group
 */
class TransformedSolution {
private:
	///The transformed orbit values.
	std::vector< std::list<double> > __transformed_orbit,

		///The transformed derivatives.
		__transformed_deriv;

	///The functions to use to transform the semimajor axis.
	std::list<const OneArgumentDiffFunction *> __a_transforms,

		///\brief The function to use  to transfom the convective zone
		///angular momentum.
		__Lconv_transforms,
		
		///\brief The function to use  to transfom the radiative zone
		///angular momentum.
		__Lrad_transforms;

	///The ages at which the orbit was tabulated.
	const std::list<double> *__ages;

	///The boundaries between consecutive trasformation functions.
	std::list<double> __change_ages;

	///The evolution mode at each tabulated age
	const std::list<EvolModeType> *__evolution_mode;

	///Applies a list of transformations to a variable.
	void transform(
			///The variable to transform.
			const std::list<double>* var,

			///The derivative of the variable to transform.
			const std::list<double>* deriv,
			
			///The list of trasformation functions.
			const std::list<const OneArgumentDiffFunction*> &transforms,
			
			///Which variable is being transformed.
			EvolVarType var_type);
public:
	///\brief Default constructor, use add_transformation repeatedly to build
	///the final solution.
	TransformedSolution() : __transformed_orbit(3), __transformed_deriv(3) {};

	///Create a single piece transformed solution.
	TransformedSolution(
			///The function to use to transform the semimajor axis.
			const OneArgumentDiffFunction &a_transform,

			///The function to use to transform the convective zone angular
			///momentum.
			const OneArgumentDiffFunction &Lconv_transform,

			///The function to use to transform the radiative zone angular
			///momentum.
			const OneArgumentDiffFunction &Lrad_transform,

			///The age at which this transformation starts to apply
			double start_age) : __transformed_orbit(3),__transformed_deriv(3)
	{add_transformation(a_transform, Lconv_transform, Lrad_transform,
			start_age);}

	///Add more pieces to the transformation.
	void add_transformation(
			///The function to use to transform the semimajor axis after the
			///last change age up to the given one.
			const OneArgumentDiffFunction &a_transform,

			///The function to use to transform the convective zone angular
			///momentum after the last change age up to the given one.
			const OneArgumentDiffFunction &Lconv_transform,

			///The function to use to transform the radiative zone angular
			///momentum after the last change age up to the given one.
			const OneArgumentDiffFunction &Lrad_transform,

			///The age up to which this transformation applies.
			double change_age);

	///\brief The solver to apply this transformation to.
	///
	///It should already have a solution stored.
	void operator()(const OrbitSolver &solver);

	///The value of a transformed variable at the tabulated ages.
	const std::list<double>
		*get_tabulated_var(
				///Which variable to return.
				EvolVarType var_type) const;

	///The derivative of a transformed variable at the tabulated ages.
	const std::list<double> *get_tabulated_var_deriv(
			EvolVarType var_type) const;

	///Returns a list of the evolution modes.
	const std::list<EvolModeType> *get_tabulated_evolution_mode() const
	{return __evolution_mode;} 
};

/**\brief Evolution mode that changes at specified ages.
 *
 * \ingroup UnitTests_group
 */
class ExpectedEvolutionMode {
private:
	///The ages at which evolution mode changes occur.
	std::list<double> __age_breaks;

	///The evolution modes that apply between consecutive __age_breaks.
	std::list<EvolModeType> __expected_mode;
public:
	///Create.
	ExpectedEvolutionMode() {}

	///Add an evolution mode that applies up to the given age.
	void add_break(double age, EvolModeType mode)
	{__age_breaks.push_back(age); __expected_mode.push_back(mode);}
	

	///The evolution mode that corresponds to the given age.
	EvolModeType operator()(double age) const;
};

/**\brief The test suite that exercises the OrbitSolver class.
 *
 * \ingroup UnitTests_group
 */
class test_OrbitSolver : public Test::Suite {
private:
	///\brief Tests the latest evolution calculated by the solver against the
	///given tracks.
	template<class SOLVER_LIKE>
	void test_solution(const SOLVER_LIKE &solver,
			const OneArgumentDiffFunction &expected_a,
			const OneArgumentDiffFunction &expected_Lconv,
			const OneArgumentDiffFunction &expected_Lrad, double min_age,
			double max_age, const ExpectedEvolutionMode &expected_mode,
			bool debug_mode=false);
protected:
	///No fixtures at this time
	void setup() {};

	///No fixtures at this time
	void tear_down() {};
public:
	///Create the test suite.
	test_OrbitSolver();

	///\brief Tests the evolution of stellar systems with disk locking times
	///longer than the stellar lifetimes for a non-evolving star.
	void test_disk_locked_no_stellar_evolution();

	///\brief Tests the evolution of stellar systems with disk locking times
	///longer than the stellar lifetimes, for evolving stars.
	void test_disk_locked_with_stellar_evolution();

	///\brief Tests the rotational evolution of a star without a planet and
	///without a disk locking phase.
	///
	///Calculates each evolution with NO_PLANET evolution mode, and with
	///FAST_PLANET/SLOW_PLANET evolution mode, but zero planet mass.
	void test_no_planet_evolution();

	///\brief Tests the evolution of the orbit plus stellar rotation,
	///starting with the planet already present and ensuring it does not die
	///or lock to the star.
	void test_unlocked_evolution();

	///\brief Tests the evolution of the orbit plus stellar rotation for the
	///case where the star is locked in synchronous rotation with the orbit.
	void test_locked_evolution();

	///\brief Tests an evolution that starts locked to a disk then
	///immediately is locked to the planet, which then falls below the Roche
	///radius.
	void test_disklocked_to_locked_to_noplanet();

	///\brief Tests an evolution that starts locked to a disk then
	///the orbital period is shorter than the stellar spin period and
	///eventually the planet falls below the Roche radius.
	void test_disklocked_to_fast_to_noplanet();

	///\brief Tests an evolution that starts locked to a disk then
	///the orbital period is shorter than the stellar spin period and
	///then is locked.
	void test_disklocked_to_fast_to_locked();

	///\brief Tests an evolution that starts locked to a disk then
	///immediately is locked to the planet and then the lock is broken
	///leaving the orbital period shorter than the stellar spin period.
	void test_disklocked_to_locked_to_fast();
};

#endif
