#ifndef __TEST_ORBIT_SOLVER_H
#define __TEST_ORBIT_SOLVER_H

#include "Common.h"
#include "../OrbitSolver.h"
#include <iostream>
#include <vector>

///A class that can be passed to the solution testing function instead of the
///solver that transforms the solutions before testing them
class TransformedSolution {
private:
	std::vector< std::list<double> > __transformed_orbit,
		__transformed_deriv;

	std::list<const OneArgumentDiffFunction *> __a_transforms,
		__Lconv_transforms, __Lrad_transforms;

	const std::list<double> *__ages;
	std::list<double> __change_ages;

	const std::list<EvolModeType> *__evolution_mode;

	void transform(const std::list<double>* var,
			const std::list<double>* deriv,
			const std::list<const OneArgumentDiffFunction*> &transforms,
			EvolVarType var_type);
public:
	TransformedSolution() : __transformed_orbit(3), __transformed_deriv(3) {};

	TransformedSolution(const OneArgumentDiffFunction &a_transform,
			const OneArgumentDiffFunction &Lconv_transform,
			const OneArgumentDiffFunction &Lrad_transform,
			double start_age) : __transformed_orbit(3),__transformed_deriv(3)
	{add_transformation(a_transform, Lconv_transform, Lrad_transform,
			start_age);}

	void add_transformation(const OneArgumentDiffFunction &a_transform,
			const OneArgumentDiffFunction &Lconv_transform,
			const OneArgumentDiffFunction &Lrad_transform,
			double change_age);

	void operator()(const OrbitSolver &solver);

	///Returns the values of the variables at the tabulated ages transformed
	///through the functions specified on construction.
	const std::list<double>
		*get_tabulated_var(EvolVarType var_type) const;

	///Returns the derivative of the independent variable at the points
	///specified by get_tabulated_indep_var transformed
	///through the functions specified on construction.
	const std::list<double> *get_tabulated_var_deriv(
			EvolVarType var_type) const;

	///Returns a list of the evolution modes.
	const std::list<EvolModeType> *get_tabulated_evolution_mode() const
	{return __evolution_mode;} 
};

class ExpectedEvolutionMode {
private:
	std::list<double> __age_breaks;
	std::list<EvolModeType> __expected_mode;
public:
	ExpectedEvolutionMode() {}
	void add_break(double age, EvolModeType mode)
	{__age_breaks.push_back(age); __expected_mode.push_back(mode);}
	EvolModeType operator()(double age) const;
};

class test_OrbitSolver : public Test::Suite {
private:
	///Tests the latest evolution calculated by the solver against the given
	///tracks.
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
	test_OrbitSolver();

	///Generates and computes the evolution of stellar systems with disk
	///locking times larger than the stellar lifetimes, creating a star
	///which does not evolve.
	void test_disk_locked_no_stellar_evolution();

	///Generates and computes the evolution of stellar systems with disk
	///locking times larger than the stellar lifetimes, creating a star
	///that evolves over time.
	void test_disk_locked_with_stellar_evolution();

	///Tests the rotational evolution of a star without a planet and without
	///a disk locking phase. Calculates each evolution with NO_PLANET 
	///evolution mode, and with FAST_PLANET/SLOW_PLANET evolution mode, but
	///zero planet mass.
	void test_no_planet_evolution();

	///Tests the evolution of the orbit plus stellar rotation, starting with
	///the planet already present and ensuring it does not die or lock to
	///the star.
	void test_unlocked_evolution();

	///Tests the evolution of the orbit plus stellar rotation for the case
	///where the star is locked in synchronous rotation with the orbit.
	void test_locked_evolution();

	void test_disklocked_to_locked_to_noplanet();

	void test_disklocked_to_fast_to_noplanet();

	void test_disklocked_to_fast_to_locked();

	void test_disklocked_to_locked_to_fast();
};

#endif
