/**\file
 *
 * \brief Defines some of the methods of the test suite that exercises the
 * OrbitSolver class and the other clasess necessary to accomplish this.
 *
 * \ingroup UnitTests_group
 */
#include "test_OrbitSolver.h"

const double tstart=2.0*min_age;
PolynomialEvolutionTrack nan_func(std::valarray<double>(NaN, 2), tstart,1.0);
PolynomialEvolutionTrack one_func(std::valarray<double>(1.0, 1), tstart,1.0);
PolynomialEvolutionTrack two_func(std::valarray<double>(2.0, 1), tstart,1.0);
PolynomialEvolutionTrack two_hundred_func(std::valarray<double>(200.0, 1), tstart,1.0);
PolynomialEvolutionTrack zero_func(std::valarray<double>(), tstart, 1.0);
const double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;

void TransformedSolution::transform(const std::list<double>* var,
		const std::list<double>* deriv,
		const std::list<const OneArgumentDiffFunction*> &transforms,
		EvolVarType var_type)
{
	std::list<double>::const_iterator change_i=__change_ages.begin();
	change_i++;
	std::list<const OneArgumentDiffFunction *>::const_iterator
		transform_i=transforms.begin();
	for(std::list<double>::const_iterator
			var_i=var->begin(), deriv_i=deriv->begin(), t_i=__ages->begin();
			var_i!=var->end(); var_i++, deriv_i++, t_i++) {
		if(change_i!=__change_ages.end() && *t_i>*change_i) {
			transform_i++; change_i++;
		}
		const FunctionDerivatives *dvar=(*transform_i)->deriv(*var_i);
		__transformed_orbit[var_type].push_back(dvar->order(0));
		__transformed_deriv[var_type].push_back(dvar->order(1)*(*deriv_i));
	}
}

void TransformedSolution::add_transformation(
		const OneArgumentDiffFunction &a_transform,
		const OneArgumentDiffFunction &Lconv_transform,
		const OneArgumentDiffFunction &Lrad_transform,
		double change_age)
{
	__a_transforms.push_back(&a_transform);
	__Lconv_transforms.push_back(&Lconv_transform);
	__Lrad_transforms.push_back(&Lrad_transform);
	__change_ages.push_back(change_age);
}

void TransformedSolution::operator()(const OrbitSolver &solver)
{
	__evolution_mode=solver.get_tabulated_evolution_mode();
	__ages=solver.get_tabulated_var(AGE);

	transform(solver.get_tabulated_var(SEMIMAJOR),
			solver.get_tabulated_var_deriv(SEMIMAJOR),
			__a_transforms, SEMIMAJOR);

	transform(solver.get_tabulated_var(LCONV),
			solver.get_tabulated_var_deriv(LCONV),
			__Lconv_transforms, LCONV);

	transform(solver.get_tabulated_var(LRAD),
			solver.get_tabulated_var_deriv(LRAD),
			__Lrad_transforms, LRAD);

}

const std::list<double> *TransformedSolution::get_tabulated_var(
		EvolVarType var_type) const
{
	if(var_type==AGE) return __ages;
	return &(__transformed_orbit[var_type]);
}

const std::list<double> *TransformedSolution::get_tabulated_var_deriv(
		EvolVarType var_type) const
{
	assert(var_type!=AGE);
	return &(__transformed_deriv[var_type]);
}

EvolModeType ExpectedEvolutionMode::operator()(double age) const
{
	std::list<EvolModeType>::const_iterator mode_i=__expected_mode.begin();
	for(std::list<double>::const_iterator age_i=__age_breaks.begin();
			age_i!=__age_breaks.end(); age_i++) {
		std::list<double>::const_iterator next_age_i=age_i;
		next_age_i++;
		if(*age_i<=age && (next_age_i==__age_breaks.end() ||
					*next_age_i>=age)) return *mode_i;
		mode_i++;
	}
	assert(false);
}

template<class SOLVER_LIKE>
void test_OrbitSolver::test_solution(const SOLVER_LIKE &solver,
		const OneArgumentDiffFunction &expected_a,
		const OneArgumentDiffFunction &expected_Lconv,
		const OneArgumentDiffFunction &expected_Lrad, double min_age,
		double max_age, const ExpectedEvolutionMode &expected_mode,
		bool debug_mode)
{
	const std::list<double>
		*tabulated_ages=solver.get_tabulated_var(AGE),
		*tabulated_a=solver.get_tabulated_var(SEMIMAJOR),
		*tabulated_Lconv=solver.get_tabulated_var(LCONV),
		*tabulated_Lrad=solver.get_tabulated_var(LRAD),
		*tabulated_da_dt=solver.get_tabulated_var_deriv(SEMIMAJOR),
		*tabulated_dLconv_dt=solver.get_tabulated_var_deriv(LCONV),
		*tabulated_dLrad_dt=solver.get_tabulated_var_deriv(LRAD);

	const std::list<EvolModeType> *tabulated_modes=
		solver.get_tabulated_evolution_mode();

	std::ostringstream msg_start;
	std::ostringstream msg;
	msg.precision(16);
	msg_start.precision(16);
	msg << msg_start.str() << tabulated_ages->size() << " tabulated ages, ";
	bool all_same_size=true;

	msg << tabulated_a->size() << " tabulated a, " << tabulated_da_dt->size()
		<< " tabulated da_dt, " << tabulated_Lconv->size()
		<< " tabulated Lconv, " << tabulated_dLconv_dt->size()
		<< " tabulated dLconv_dt, " << tabulated_Lrad->size()
		<< " tabulated Lrad, " << tabulated_dLrad_dt->size()
		<< " tabulated dLrad_dt, " << tabulated_modes->size()
		<< " tabulated modes";
	unsigned num_ages=tabulated_ages->size();
	all_same_size=tabulated_a->size()==num_ages &&
		tabulated_da_dt->size()==num_ages &&
		tabulated_Lconv->size()==num_ages &&
		tabulated_dLconv_dt->size()==num_ages &&
		tabulated_Lrad->size()==num_ages &&
		tabulated_dLrad_dt->size()==num_ages &&
		tabulated_modes->size()==num_ages;
	if(debug_mode) std::cout << msg.str() << std::endl;
	TEST_ASSERT_MSG(all_same_size, msg.str().c_str());

	msg.str("");
	msg << "Expected age range: " << min_age << " to " << max_age
		<< ", actual age range: " << tabulated_ages->front()
		<< " to " << tabulated_ages->back();
	if(debug_mode) std::cout << msg.str() << std::endl;
	TEST_ASSERT_MSG(tabulated_ages->front()==min_age &&
			tabulated_ages->back()==max_age, msg.str().c_str());

	std::list<double>::const_iterator a_iter=tabulated_a->begin(),
		Lconv_iter=tabulated_Lconv->begin(), Lrad_iter=tabulated_Lrad->begin(),
		da_dt_iter=tabulated_da_dt->begin(),
		dLconv_dt_iter=tabulated_dLconv_dt->begin(),
		dLrad_dt_iter=tabulated_dLrad_dt->begin();

	std::list<EvolModeType>::const_iterator mode_iter=
		tabulated_modes->begin();

	for(std::list<double>::const_iterator age_i=tabulated_ages->begin();
			age_i!=tabulated_ages->end(); age_i++) {
		const FunctionDerivatives *answer_a=expected_a.deriv(*age_i),
			  *answer_Lconv=expected_Lconv.deriv(*age_i),
			  *answer_Lrad=expected_Lrad.deriv(*age_i);
		double a_value=answer_a->order(0), a_deriv=answer_a->order(1),
			   Lconv_value=answer_Lconv->order(0),
			   Lconv_deriv=answer_Lconv->order(1),
			   Lrad_value=answer_Lrad->order(0),
			   Lrad_deriv=answer_Lrad->order(1);
		EvolModeType mode_value=expected_mode(*age_i);
		delete answer_a; delete answer_Lconv; delete answer_Lrad;

		double a=*a_iter, Lconv=*Lconv_iter, Lrad=*Lrad_iter,
			   da_dt=*da_dt_iter, dLconv_dt=*dLconv_dt_iter,
			   dLrad_dt=*dLrad_dt_iter;
		EvolModeType mode=*mode_iter;

		for(std::list<double>::const_iterator next_age_i=age_i;
				++next_age_i!=tabulated_ages->end() &&
				(*next_age_i)==(*age_i);) {
			mode_iter++; a_iter++; Lconv_iter++; Lrad_iter++;
			da_dt_iter++; dLconv_dt_iter++; dLrad_dt_iter++; age_i++;
			if(mode_value==*mode_iter) {
				mode=*mode_iter;
				a=*a_iter;
				Lconv=*Lconv_iter;
				Lrad=*Lrad_iter;
				da_dt=*da_dt_iter;
				dLconv_dt=*dLconv_dt_iter;
				dLrad_dt=*dLrad_dt_iter;
			}
		}

		std::ostringstream age_msg_start;
		age_msg_start.precision(16);
		age_msg_start << msg_start.str() << "age=" << *age_i
			<< ", mode=" << mode << ", a=" << a << ", da_dt=" << da_dt
			<< ", Lconv=" << Lconv << ", dLconv_dt=" << dLconv_dt
			<< ", Lrad=" << Lrad << ", dLrad_dt=" << dLrad_dt;
		msg.str("");
		msg << age_msg_start.str() << " age is out of range.";
		TEST_ASSERT_MSG(*age_i>=min_age && *age_i<=max_age,
				msg.str().c_str());


		msg.str("");
		msg << age_msg_start.str() << ": mode is not "
			<< expected_mode(*age_i) << ", but " << mode;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(mode_value==mode, msg.str().c_str());


		msg.str("");
		msg << age_msg_start.str() << ": a is not "
			<< a_value << ", difference=" << a-a_value;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(check_diff(a, a_value, 1e-5, 0.0),
				msg.str().c_str());

		msg.str("");
		msg << age_msg_start.str() << ": Lconv is not "
			<< Lconv_value << ", difference=" << Lconv-Lconv_value;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(check_diff(Lconv, Lconv_value, 1e-5, 0.0),
				msg.str().c_str());

		msg.str("");
		msg << age_msg_start.str() << ": Lrad is not "
			<< Lrad_value << ", difference=" << Lrad-Lrad_value;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(check_diff(Lrad, Lrad_value, 1e-5, 0.0),
				msg.str().c_str());

		msg.str("");
		msg << age_msg_start.str() << ": da_dt is not "
			<< a_deriv << ", difference=" << da_dt-a_deriv;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(check_diff(da_dt, a_deriv, 1e-5, 0.0),
				msg.str().c_str());

		msg.str("");
		msg << age_msg_start.str() << ": dLconv_dt is not "
			<< Lconv_deriv << ", difference="
			<< dLconv_dt-Lconv_deriv;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(check_diff(dLconv_dt, Lconv_deriv, 1e-5, 0.0),
				msg.str().c_str());

		msg.str("");
		msg << age_msg_start.str() << ": dLrad_dt is not "
			<< Lrad_deriv << ", difference="
			<< dLrad_dt-Lrad_deriv;
		if(debug_mode) std::cout << msg.str() << std::endl;
		TEST_ASSERT_MSG(check_diff(dLrad_dt, Lrad_deriv, 1e-5, 0.0),
				msg.str().c_str());

		mode_iter++; a_iter++; Lconv_iter++; Lrad_iter++;
		da_dt_iter++; dLconv_dt_iter++; dLrad_dt_iter++;
	}
}

void test_OrbitSolver::test_disk_locked_no_stellar_evolution()
{
	try {
		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_DISK);
		MockStellarEvolution no_evol(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		Star star(1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, Inf, no_evol, 1.0, 1.0,
				1.0);
		Planet planet(&star, 1.0, 1.0, 1.0);
		StellarSystem system(&star, &planet);
		OrbitSolver solver(tstart, 1.0, 5e-5);
		solver(system, Inf, Inf, 1.0, tstart, LOCKED_TO_DISK,
				std::valarray<double>(1.0, 1));
		test_solution(solver, nan_func,
				PolynomialEvolutionTrack(std::valarray<double>(1.0, 1),
					tstart, 1.0),
				PolynomialEvolutionTrack(std::valarray<double>(1.0, 1),
					tstart, 1.0),
				tstart, 1.0, expected_mode);

		solver(system, Inf, Inf, 1.0, tstart, LOCKED_TO_DISK,
				std::valarray<double>(0.0, 1));
		PolynomialEvolutionTrack *one_poly=new PolynomialEvolutionTrack(
				std::valarray<double>(1.0, 1), tstart, 1.0);
		test_solution(solver, nan_func,
				PolynomialEvolutionTrack(std::valarray<double>(1.0, 1),
					tstart, 1.0),
				ExponentialPlusFunc(one_poly, -std::exp(tstart/2.0), -0.5),
				tstart, 1.0, expected_mode);
		delete one_poly;
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
					ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
					ex.what()).c_str());
	}
}

void test_OrbitSolver::test_disk_locked_with_stellar_evolution()
{
	try {
		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_DISK);
		MockStellarEvolution evol1(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 2),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 2),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		Star star1(1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, Inf, evol1, 1.0, 1.0,
				1.0);
		Planet planet1(&star1, 1.0, 1.0, 1.0);
		StellarSystem system1(&star1, &planet1);
		OrbitSolver solver(tstart, 1.0, 5e-5);

		solver(system1, Inf, Inf, 1.0, tstart, LOCKED_TO_DISK,
				std::valarray<double>(tstart-1+std::exp(-tstart/2.0), 1));
		std::valarray<double> temp_array=std::valarray<double>(1.0, 2);
		temp_array[0]=-1;
		PolynomialEvolutionTrack *temp_poly=
			new PolynomialEvolutionTrack(temp_array, tstart, 1.0);
		test_solution(solver, nan_func,
				PolynomialEvolutionTrack(std::valarray<double>(1.0, 2),
					tstart, 1.0),
				ExponentialPlusFunc(temp_poly, 1, -0.5),
				tstart, 1.0, expected_mode);

		solver(system1, Inf, Inf, 1.0, tstart, LOCKED_TO_DISK,
				std::valarray<double>(tstart-1+2.0*std::exp(-tstart/2.0), 1));
		test_solution(solver, nan_func,
				PolynomialEvolutionTrack(std::valarray<double>(1.0, 2),
					tstart, 1.0),
				ExponentialPlusFunc(temp_poly, 2, -0.5),
				tstart, 1.0, expected_mode);
		delete temp_poly;

		std::valarray< std::valarray<double> >
			Ic_arr(std::valarray<double>(1.0, 1), 2);
		Ic_arr[1][0]=-1.0/6.0;
		MockStellarEvolution evol2(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				Ic_arr,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0/6.0, 1), 2),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(0.5, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 2));
		Star star2(1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, Inf, evol2, 1.0, 1.0,
				1.0);
		Planet planet2(&star2, 1.0, 1.0, 1.0);
		StellarSystem system2(&star2, &planet2);
		solver(system2, Inf, Inf, 1.0, tstart, LOCKED_TO_DISK,
				std::valarray<double>((tstart/2.0+1.0)/6.0, 1));
		std::valarray<double> Lc_coef(1.0, 2); Lc_coef[1]=-1.0/6.0;
		test_solution(solver, nan_func,
				PolynomialEvolutionTrack(Lc_coef, tstart, 1.0),
				PolynomialEvolutionTrack(
					std::valarray<double>(1.0/6.0, 2), tstart, 1.0),
				tstart, 1.0, expected_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
					ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
					ex.what()).c_str());
	}
}

void test_OrbitSolver::test_no_planet_evolution()
{
	try {
		ExpectedEvolutionMode no_planet_mode, slow_planet_mode,
							  fast_planet_mode;
		no_planet_mode.add_break(tstart, NO_PLANET);
		slow_planet_mode.add_break(tstart, SLOW_PLANET);
		fast_planet_mode.add_break(tstart, FAST_PLANET);
		const double rt2=std::sqrt(2.0);
		MockStellarEvolution evol1(-1,std::valarray< std::valarray<double> >(
										std::valarray<double>(1.0, 1), 1),
									std::valarray< std::valarray<double> >(
										std::valarray<double>(1.0, 1), 1),
									std::valarray< std::valarray<double> >(
										std::valarray<double>(1.0, 1), 1),
									std::valarray< std::valarray<double> >(
										std::valarray<double>(1.0, 1), 1),
									std::valarray< std::valarray<double> >(
										std::valarray<double>(1.0, 1), 1)),
							 evol2(-1,std::valarray< std::valarray<double> >(
										 std::valarray<double>(1.0, 1), 1),
									 std::valarray< std::valarray<double> >(
										 std::valarray<double>(1.0, 1), 2),
									 std::valarray< std::valarray<double> >(
										 std::valarray<double>(1.0, 1), 2),
									 std::valarray< std::valarray<double> >(
										 std::valarray<double>(1.0, 1), 1),
									 std::valarray< std::valarray<double> >(
										 std::valarray<double>(1.0, 1), 1));
		Star star_no_wind_no_coupling(1.0, 1.0, 0.0, 2.0, Inf, 0.0, 0.0,
				0.0, evol1, 1.0, 1.0, 1.0);
		Planet planet0(&star_no_wind_no_coupling, 0.0, 0.0, 1.0);
		Planet planet1(&star_no_wind_no_coupling, 1.0, 1.0, 1.0);
		StellarSystem system0(&star_no_wind_no_coupling, &planet0);
		StellarSystem system1(&star_no_wind_no_coupling, &planet1);
		OrbitSolver solver(tstart, 1.0, 5e-7);


		std::valarray<double> initial_orbit(0.0, 3);
		solver(system1, Inf, Inf, 2.0, tstart, NO_PLANET,
				std::valarray<double>(0.0, 2));
		test_solution(solver, nan_func, zero_func, zero_func, tstart, 1.0,
				no_planet_mode);
		initial_orbit[0]=std::pow(2.0, 6.5);
		solver(system0, Inf, 0.0, 2.0, tstart, FAST_PLANET, initial_orbit);
		test_solution(solver, two_func, zero_func, zero_func, tstart, 1.0,
				fast_planet_mode);

		solver(system1, Inf, Inf, 2.0, tstart, NO_PLANET,
				std::valarray<double>(1.0, 2));
		test_solution(solver, nan_func, one_func, one_func, tstart, 1.0,
				no_planet_mode);
		initial_orbit.resize(3, 1.0); initial_orbit[0]=std::pow(2.0, 6.5);
		solver(system0, Inf, 0.0, 2.0, tstart, FAST_PLANET, initial_orbit);
		test_solution(solver, two_func, one_func, one_func, tstart, 1.0,
				fast_planet_mode);
		initial_orbit[0]=std::pow(200.0, 6.5);
		solver(system0, Inf, 0.0, 200.0, tstart, SLOW_PLANET, initial_orbit);
		test_solution(solver, two_hundred_func, one_func, one_func, tstart, 1.0,
				slow_planet_mode);

		initial_orbit.resize(2, 0.0); initial_orbit[0]=1;
		solver(system1, Inf, Inf, 2.0, tstart, NO_PLANET,
				initial_orbit);
		test_solution(solver, nan_func, one_func, zero_func, tstart, 1.0,
				no_planet_mode);
		initial_orbit.resize(3, 1.0);
		initial_orbit[0]=std::pow(2.0, 6.5);
		initial_orbit[2]=0.0;
		solver(system0, Inf, 0.0, 2.0, tstart, FAST_PLANET,
				initial_orbit);
		test_solution(solver, two_func, one_func, zero_func, tstart, 1.0,
				fast_planet_mode);
		initial_orbit.resize(3, 1.0);
		initial_orbit[0]=std::pow(200.0, 6.5);
		initial_orbit[2]=0.0;
		solver(system0, Inf, 0.0, 200.0, tstart, SLOW_PLANET,
				initial_orbit);
		test_solution(solver, two_hundred_func, one_func, zero_func, tstart,
				1.0, slow_planet_mode);

		Star star_no_wind(1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0,
				0.0, evol1, 1.0, 1.0, 1.0);
		initial_orbit.resize(2);
		initial_orbit[0]=0.5*(1.0+std::exp(-tstart));
		initial_orbit[1]=0.5*(1.0-std::exp(-tstart));
		StellarSystem system2(&star_no_wind, &planet1);
		StellarSystem system2_0(&star_no_wind, &planet0);
		PolynomialEvolutionTrack half_func(std::valarray<double>(0.5, 1),
				tstart, 1.0);
		solver(system2, Inf, Inf, 2.0, tstart, NO_PLANET, initial_orbit);
		test_solution(solver, nan_func,
				ExponentialPlusFunc(&half_func, 0.5, -1.0),
				ExponentialPlusFunc(&half_func, -0.5, -1.0), tstart, 1.0,
				no_planet_mode);
		initial_orbit.resize(3);
		initial_orbit[0]=std::pow(2.0, 6.5);
		initial_orbit[1]=0.5*(1.0+std::exp(-tstart));
		initial_orbit[2]=0.5*(1.0-std::exp(-tstart));
		solver(system2_0, Inf, 0.0, 2.0, tstart, FAST_PLANET, initial_orbit);
		test_solution(solver, two_func,
				ExponentialPlusFunc(&half_func, 0.5, -1.0),
				ExponentialPlusFunc(&half_func, -0.5, -1.0), tstart, 1.0,
				fast_planet_mode);
		initial_orbit[0]=std::pow(200.0, 6.5);
		solver(system2_0, Inf, 0.0, 200.0, tstart, SLOW_PLANET, initial_orbit);
		test_solution(solver, two_hundred_func,
				ExponentialPlusFunc(&half_func, 0.5, -1.0),
				ExponentialPlusFunc(&half_func, -0.5, -1.0), tstart, 1.0,
				slow_planet_mode);

		Star star_no_coupling(1.0, 1.0, 2.0, std::sqrt(0.5), Inf, 0.0, 1.0,
				0.0, evol2, 1.0, 1.0, 1.0);
		initial_orbit.resize(2);
		initial_orbit[0]=1.0/(1.0+tstart);
		initial_orbit[1]=1.0;
		StellarSystem system3(&star_no_coupling, &planet1);
		StellarSystem system3_0(&star_no_coupling, &planet0);
		std::valarray<double> late_denom_coef(3);
		late_denom_coef[0]=2.0*rt2-2.0;
		late_denom_coef[1]=4.0*rt2;
		late_denom_coef[2]=2.0*rt2;
		PolynomialEvolutionTrack one_func_early(std::valarray<double>(1.0,1),
				tstart, std::pow(2.0,0.25)-1),
				one_plus_t(std::valarray<double>(1.0, 2), tstart, 1.0),
				late_denom2(late_denom_coef, std::pow(2.0, 0.25)-1, 1.0);
		FunctionToPower late_denom(&late_denom2, 0.5);
		FunctionRatio early_solution(&one_func_early, &one_plus_t),
					  late_solution(&one_plus_t, &late_denom);
		PiecewiseFunction full_solution;
		full_solution.add_piece(&early_solution);
		full_solution.add_piece(&late_solution);
		solver(system3, Inf, Inf, 2.0, tstart, NO_PLANET, initial_orbit);
		test_solution(solver, nan_func, full_solution, one_func,
				tstart, 1.0, no_planet_mode);
		initial_orbit.resize(3, 1.0);
		initial_orbit[0]=std::pow(2.0, 6.5);
		initial_orbit[1]=1.0/(1.0+tstart);
		solver(system3_0, Inf, 0.0, 2.0, tstart, FAST_PLANET,
				initial_orbit);
		test_solution(solver, two_func, full_solution, one_func,
				tstart, 1.0, fast_planet_mode);
		initial_orbit[0]=std::pow(200.0, 6.5);
		solver(system3_0, Inf, 0.0, 200.0, tstart, SLOW_PLANET,
				initial_orbit);
		test_solution(solver, two_hundred_func, full_solution, one_func,
				tstart, 1.0, slow_planet_mode);

		Star star(1.0, 1.0, 100.0, 0.1, 1.0, 0.0, 0.0, 0.0, evol1, 1.0, 1.0,
				1.0);
		double b1=(std::sqrt(2)-1)/(2.0*rt2),
			   b2=(std::sqrt(2)+1)/(2.0*rt2);
		initial_orbit.resize(2);
		initial_orbit[0]=b1*std::exp(-tstart/(2.0+rt2)) +
			b2*std::exp(-tstart/(2.0-rt2));
		initial_orbit[1]=0.5/rt2*std::exp(-1.0/(2.0+rt2)*(tstart))-
				0.5/rt2*std::exp(-1.0/(2.0-rt2)*(tstart));
		StellarSystem system4(&star, &planet1);
		StellarSystem system4_0(&star, &planet0);

		ExponentialPlusFunc 
			Lc1(&zero_func, b1, -1.0/(2.0+rt2)),
			Lc2(&zero_func, b2, -1.0/(2.0-rt2)),
			Lr1(&zero_func, 0.5/rt2, -1.0/(2.0+rt2)),
			Lr2(&zero_func, -0.5/rt2, -1.0/(2.0-rt2));

		solver(system4, Inf, Inf, 2.0, tstart, NO_PLANET,
				initial_orbit);
		test_solution(solver, nan_func, FuncPlusFunc(&Lc1, &Lc2),
				FuncPlusFunc(&Lr1, &Lr2), tstart, 1.0, no_planet_mode);
		initial_orbit.resize(3);
		initial_orbit[0]=std::pow(2.0, 6.5);
		initial_orbit[1]=b1*std::exp(-tstart/(2.0+rt2)) +
			b2*std::exp(-tstart/(2.0-rt2));
		initial_orbit[2]=0.5/rt2*std::exp(-1.0/(2.0+rt2)*(tstart))-
				0.5/rt2*std::exp(-1.0/(2.0-rt2)*(tstart));
		solver(system4_0, Inf, 0.0, 2.0, tstart, FAST_PLANET,
				initial_orbit);
		test_solution(solver, two_func, FuncPlusFunc(&Lc1, &Lc2),
				FuncPlusFunc(&Lr1, &Lr2), tstart, 1.0, fast_planet_mode);
		initial_orbit[0]=std::pow(200.0, 6.5);
		solver(system4_0, Inf, 0.0, 200.0, tstart, SLOW_PLANET,
				initial_orbit);
		test_solution(solver, two_hundred_func, FuncPlusFunc(&Lc1, &Lc2),
				FuncPlusFunc(&Lr1, &Lr2), tstart, 1.0, slow_planet_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}
}

void test_OrbitSolver::test_unlocked_evolution()
{
	try {
		ExpectedEvolutionMode slow_planet_mode,
							  fast_planet_mode;
		slow_planet_mode.add_break(tstart, SLOW_PLANET);
		fast_planet_mode.add_break(tstart, FAST_PLANET);
		MockStellarEvolution no_evol_highI(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1)),
							 no_evol_lowI(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0e-9, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		const double mplanet=100;
		double Q=1e8*mplanet, alpha=-4.5*std::sqrt(AstroConst::G/
				(AstroConst::solar_radius*AstroConst::solar_mass))*
			mplanet*AstroConst::jupiter_mass/Q*AstroConst::Gyr/
			AstroConst::solar_radius,
			a6p5_offset=std::pow(2.0, 6.5)-6.5*alpha,
			a0=std::pow(a6p5_offset, 1.0/6.5),
			Lscale=-AstroConst::jupiter_mass*mplanet/
				std::pow(AstroConst::solar_radius, 1.5)*
				std::sqrt(AstroConst::G/
						(mplanet*AstroConst::jupiter_mass+
						 AstroConst::solar_mass))*
				AstroConst::day;
		Star star_highI_no_wind_no_coupling(1.0, Q, 0.0, 1.0, Inf, 0.0, 0.0,
				0.0, no_evol_highI, 1.0, 1.0, 1.0);
		Star star_lowI_no_wind_no_coupling(1.0, Q, 0.0, 1.0, Inf, 0.0, 0.0,
				0.0, no_evol_lowI, 1.0, 1.0, 1.0);
		Planet planet_highI(&star_highI_no_wind_no_coupling, mplanet,
				0.0001, 1.0),
			   planet_lowI(&star_lowI_no_wind_no_coupling, mplanet,
					   0.0001, 1.0);
		StellarSystem system_highI(&star_highI_no_wind_no_coupling,
								&planet_highI),
					  system_lowI(&star_lowI_no_wind_no_coupling,
							  &planet_lowI);
		OrbitSolver solver(tstart, 1.0, 1e-7);
		std::valarray<double> initial_orbit(3);
		initial_orbit[0]=a6p5_offset+6.5*tstart*alpha;
		initial_orbit[1]=Lscale*(std::pow(initial_orbit[0], 1.0/13.0) -
				std::sqrt(a0));
		initial_orbit[2]=0.0;
		std::valarray<double> a6p5_poly_coef(2);
		a6p5_poly_coef[0]=a6p5_offset; a6p5_poly_coef[1]=6.5*alpha;
		PolynomialEvolutionTrack a6p5_evol(a6p5_poly_coef, tstart, 1.0);
		FunctionToPower a_evol(&a6p5_evol, 1.0/6.5),
						sqrta_evol(&a6p5_evol, 1.0/13.0);
		ExponentialPlusFunc Lconv_unscaled(&sqrta_evol, -std::sqrt(a0), 0);
		solver(system_highI, Inf, tstart, a0, tstart, FAST_PLANET, initial_orbit);
		test_solution(solver, a_evol,
				ScaledFunction(&Lconv_unscaled, Lscale), zero_func,
				tstart, 1.0, fast_planet_mode);

		a0=2.6;
		a6p5_offset=std::pow(a0, 6.5);
		initial_orbit[0]=a6p5_offset-6.5*tstart*alpha;
		initial_orbit[1]=Lscale*(std::pow(initial_orbit[0], 1.0/13.0) -
				std::sqrt(a0));
		initial_orbit[2]=0.0;
		a6p5_poly_coef[0]=a6p5_offset; a6p5_poly_coef[1]=-6.5*alpha;
		PolynomialEvolutionTrack a6p5_evol_slow(a6p5_poly_coef, tstart, 1.0);
		FunctionToPower a_evol_slow(&a6p5_evol_slow, 1.0/6.5),
						sqrta_evol_slow(&a6p5_evol_slow, 1.0/13.0);
		ExponentialPlusFunc Lconv_unscaled_slow(&sqrta_evol_slow,
				-std::sqrt(a0), 0);
		solver(system_highI, Inf, tstart, a0, tstart, SLOW_PLANET,
				initial_orbit, true);
		test_solution(solver, a_evol_slow,
				ScaledFunction(&Lconv_unscaled_slow, Lscale), zero_func,
				tstart, 1.0, slow_planet_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}
}

///\brief The equation that should be solved in order to get the semimamjor
///axis at a gien time.
///
///The paramaters should be: alpha, beta, kappa, C, t
double locked_unsat_eq(double a, void *params)
{
	double *dbl_par=static_cast<double *>(params);
	double alpha=dbl_par[0], beta=dbl_par[1], kappa=dbl_par[2],
		   c=dbl_par[3], t=dbl_par[4];
	return 0.1*alpha*std::pow(a, 5) - 0.5*beta*std::pow(a, 3) + kappa*t - c;
}

///\mbrief The equation that should be solved in order to get the semimamjor
///axis at a given time for the locked evolution when the wind is saturated.
///
///The paramaters should be: alpha, beta, kappa, C, t
double locked_sat_eq(double a, void *params)
{
	double *dbl_par=static_cast<double *>(params);
	double alpha=dbl_par[0], beta=dbl_par[1], kappa=dbl_par[2],
		   c=dbl_par[3], t=dbl_par[4];
	return alpha*a*a/4.0 - 1.5*beta*std::log(a) + kappa*t - c;
}

///\brief The derivative of the equation that should be solved in order to
///get the semimamjor axis at a given time for the locked evolution when the
///wind is not saturated.
///
///The paramaters should be: alpha, beta.
double locked_unsat_deriv(double a, void *params)
{
	double *dbl_par=static_cast<double *>(params);
	double alpha=dbl_par[0], beta=dbl_par[1];
	return 0.5*alpha*std::pow(a, 4) - 1.5*beta*std::pow(a, 2);
}

///\brief The derivative of the equation that should be solved in order to
///get the semimamjor axis at a gien time for the locked evolution when the
///wind is saturated.
///
///The paramaters should be: alpha, beta.
double locked_sat_deriv(double a, void *params)
{
	double *dbl_par=static_cast<double *>(params);
	double alpha=dbl_par[0], beta=dbl_par[1];
	return alpha*a/2.0 - 1.5*beta/a;
}

///\brief The equation and its derivative that should be solved in order to
///get the semimamjor axis at a given time for the locked evolution when the
///wind is not saturated.
///
///The paramaters should be: alpha, beta, kappa, C, t
void locked_unsat_eq_deriv(double a, void *params, double *f, double *df)
{
	double *dbl_par=static_cast<double *>(params);
	double alpha=dbl_par[0], beta=dbl_par[1], kappa=dbl_par[2],
		   c=dbl_par[3], t=dbl_par[4];
	*f=0.1*alpha*std::pow(a, 5) - 0.5*beta*std::pow(a, 3) + kappa*t - c;
	*df=0.5*alpha*std::pow(a, 4) - 1.5*beta*std::pow(a, 2);
}

///\brief The equation and its derivative that should be solved in order to
///get the semimamjor axis at a given time for the locked evolution when the
///wind is saturated.
///
///The paramaters should be: alpha, beta, kappa, C, t
void locked_sat_eq_deriv(double a, void *params, double *f, double *df)
{
	double *dbl_par=static_cast<double *>(params);
	double alpha=dbl_par[0], beta=dbl_par[1], kappa=dbl_par[2],
		   c=dbl_par[3], t=dbl_par[4];
	*f=alpha*a*a/4.0 - 1.5*beta*std::log(a) + kappa*t - c;
	*df=alpha*a/2.0 - 1.5*beta/a;
}

void test_OrbitSolver::test_locked_evolution()
{
	try {
		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_PLANET);
		const double Ic=0.001, Kwind=1e-3, Kwind_s=1;
		MockStellarEvolution no_evol(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(Ic, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		double a1=3, Lscale=AstroConst::jupiter_mass/
			std::pow(AstroConst::solar_radius, 1.5)*
			std::sqrt(AstroConst::G/
					(AstroConst::jupiter_mass+AstroConst::solar_mass))*
			AstroConst::day,
			beta=Ic*std::sqrt(AstroConst::G*(AstroConst::solar_mass+
						AstroConst::jupiter_mass))*AstroConst::day/
				std::pow(AstroConst::solar_radius, 1.5),
			wsat_s=0.1, //Must be adjusted if a1 is adjusted
			wsat=Inf,
			kappa=Kwind*std::pow(AstroConst::G*
					(AstroConst::solar_mass+AstroConst::jupiter_mass)/
					std::pow(AstroConst::solar_radius, 3), 1.5)*
				std::pow(AstroConst::day, 3),
			kappa_s=Kwind_s*wsat_s*wsat_s*std::sqrt(AstroConst::G*
					(AstroConst::solar_mass+AstroConst::jupiter_mass)/
					std::pow(AstroConst::solar_radius, 3))*AstroConst::day,
			int_const=Lscale/10.0*std::pow(a1, 5) - beta/2.0*std::pow(a1, 3)+
				kappa,
			int_const_s=Lscale/4.0*std::pow(a1, 2) - 3.0*beta/2.0*std::log(a1)+
				kappa_s;
		double solver_params[]={Lscale, beta, kappa, int_const, 0.0};
		double a0=solve(a1, 0.0, 1e-9, &locked_unsat_eq, &locked_unsat_deriv,
				&locked_unsat_eq_deriv, static_cast<void*>(solver_params));
		solver_params[4]=tstart;
		double astart=solve(a0, 0.0, 1e-9, &locked_unsat_eq,
				&locked_unsat_deriv, &locked_unsat_eq_deriv,
				static_cast<void*>(solver_params));
		solver_params[2]=kappa_s; solver_params[3]=int_const_s;
		solver_params[4]=0;
		double a0_s=solve(a1, 0.0, 1e-9, &locked_sat_eq,
				&locked_sat_deriv, &locked_sat_eq_deriv,
				static_cast<void*>(solver_params));
		solver_params[4]=tstart;
		double astart_s=solve(a1, 0.0, 1e-9, &locked_sat_eq,
				&locked_sat_deriv, &locked_sat_eq_deriv,
				static_cast<void*>(solver_params));
		solver_params[4]=1.0;
/*		double w0=std::sqrt(
					   AstroConst::G*
					   (AstroConst::solar_mass+AstroConst::jupiter_mass)/
					   std::pow(a0*AstroConst::solar_radius, 3))*
				   AstroConst::day,
			   w0_s=std::sqrt(
					   AstroConst::G*
					   (AstroConst::solar_mass+AstroConst::jupiter_mass)/
					   std::pow(a0_s*AstroConst::solar_radius, 3))*
				   AstroConst::day,
			   w1=std::sqrt(
					   AstroConst::G*
					   (AstroConst::solar_mass+AstroConst::jupiter_mass)/
					   std::pow(a1_check*AstroConst::solar_radius, 3))*
				   AstroConst::day,
			   w1_s=std::sqrt(
					   AstroConst::G*
					   (AstroConst::solar_mass+AstroConst::jupiter_mass)/
					   std::pow(a1_check_s*AstroConst::solar_radius, 3))*
				   AstroConst::day;*/
		std::valarray<double> a_transform_coef(0.0, 6), t_coef(2), t_coef_s(2);
		a_transform_coef[5]=Lscale/10.0; a_transform_coef[3]=-beta/2.0;
		t_coef[0]=int_const; t_coef[1]=-kappa;
		t_coef_s[0]=int_const_s; t_coef_s[1]=-kappa_s;
		std::valarray<double> Lconv_term1_coef(0.0, 2), Lconv_term2_coef(0.0, 3),
			identity_coef(0.0, 2), a_term1_coef_s(0.0, 3),
			Lconv_term1_coef_s(0.0, 2), Lc_beta_coef_s(0.0, 2);
		Lconv_term1_coef[1]=1.0/beta*std::pow(10.0/Lscale, 3.0/10.0);
		Lconv_term2_coef[2]=-2.0/std::pow(beta, 3);
		identity_coef[1]=1;
		a_term1_coef_s[2]=Lscale/4.0;
		Lconv_term1_coef_s[1]=1.0/beta*std::pow(4.0/Lscale, 3.0/4.0);
		Lc_beta_coef_s[1]=1.0/beta;
		PolynomialEvolutionTrack a_transform(a_transform_coef, 0.0, Inf),
								 a_transform1_s(a_term1_coef_s, 0.0, Inf),
								 transformed_a_evol(t_coef, tstart, 1.0),
								 transformed_a_evol_s(t_coef_s, tstart, 1.0),
								 Lconv_term1_poly(Lconv_term1_coef, 0.0,Inf),
								 Lconv_term2_poly(Lconv_term2_coef, 0.0,Inf),
								 Lconv_term1_poly_s(Lconv_term1_coef_s, 0.0,
										 Inf),
								 Lc_beta_s(Lc_beta_coef_s, 0.0, Inf),
								 identity(identity_coef, 0.0, Inf);
		FunctionToPower L_transform1(&Lconv_term1_poly, -10.0/3.0),
						L_transform2(&Lconv_term2_poly, -1.0),
						L_transform1_s(&Lconv_term1_poly_s, -4.0/3.0);
		LogFunction log_a(&identity), log_Lc_beta_s(&Lc_beta_s);
		ScaledFunction a_transform2_s(&log_a, -1.5*beta),
					   L_transform2_s(&log_Lc_beta_s, beta);
		FuncPlusFunc L_transform(&L_transform1, &L_transform2),
					 a_transform_s(&a_transform1_s, &a_transform2_s),
					 L_transform_s(&L_transform1_s, &L_transform2_s);
		std::valarray<double> initial_orbit(2);

		initial_orbit[0]=astart;
		initial_orbit[1]=0.0;
		std::cout.precision(16);
		std::cout.setf(std::ios_base::scientific);
		Star star_not_saturated_wind_no_coupling(1.0, 0.0, Kwind, wsat, Inf,
				0.0, 0.0, 0.0, no_evol);
		Planet planet1(&star_not_saturated_wind_no_coupling, 1.0, 1.0, 1.0);
		StellarSystem system1(&star_not_saturated_wind_no_coupling, &planet1);
		OrbitSolver solver(tstart, 1.0, 1e-9);
		solver(system1, Inf, 0.0, a0, tstart, LOCKED_TO_PLANET, initial_orbit, true);
		TransformedSolution to_check(a_transform, L_transform, identity, -Inf);
		to_check(solver);
		test_solution(to_check, transformed_a_evol, transformed_a_evol,
				zero_func, tstart, 1.0, expected_mode);

		initial_orbit[0]=astart_s;
		initial_orbit[1]=0.0;
		Star star_saturated_wind_no_coupling(1.0, 0.0, Kwind_s, wsat_s, Inf,
				0.0, 0.0, 0.0, no_evol);
		Planet planet2(&star_saturated_wind_no_coupling, 1.0, 1.0, 1.0);
		StellarSystem system2(&star_saturated_wind_no_coupling, &planet2);
		solver(system2, Inf, 0.0, a0_s, tstart, LOCKED_TO_PLANET,
				initial_orbit, true);
		TransformedSolution to_check_s(a_transform_s, L_transform_s,
				identity, -Inf);
		to_check_s(solver);
		test_solution(to_check_s, transformed_a_evol_s, transformed_a_evol_s,
				zero_func, tstart, 1.0, expected_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}
}

void test_OrbitSolver::test_disklocked_to_locked_to_noplanet()
{
	try {
		const double Ic=0.001, Kwind=8e-6, tdisk=1, tfinal=2;
		MockStellarEvolution no_evol(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(Ic, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		double afinal=1.75, Lscale=AstroConst::jupiter_mass/
			std::pow(AstroConst::solar_radius, 1.5)*
			std::sqrt(AstroConst::G/
					(AstroConst::jupiter_mass+AstroConst::solar_mass))*
			AstroConst::day,
			beta=Ic*std::sqrt(AstroConst::G*(AstroConst::solar_mass+
						AstroConst::jupiter_mass))*AstroConst::day/
				std::pow(AstroConst::solar_radius, 1.5),
			wsat=Inf,
			kappa=Kwind*std::pow(AstroConst::G*
					(AstroConst::solar_mass+AstroConst::jupiter_mass)/
					std::pow(AstroConst::solar_radius, 3), 1.5)*
				std::pow(AstroConst::day, 3),
			int_const=Lscale/10.0*std::pow(afinal, 5) -
				beta/2.0*std::pow(afinal, 3) + kappa*tfinal;
		double solver_params[]={Lscale, beta, kappa, int_const, tdisk};
		double ainitial=solve(afinal, 0.0, 1e-9, &locked_unsat_eq,
				&locked_unsat_deriv, &locked_unsat_eq_deriv,
				static_cast<void*>(solver_params)),
			   wdisk=std::sqrt(
					   AstroConst::G*
					   (AstroConst::solar_mass+AstroConst::jupiter_mass)/
					   std::pow(ainitial*AstroConst::solar_radius, 3))*
				   AstroConst::day;
		Star star_not_saturated_wind_no_coupling(1.0, 1.0e-10, Kwind, wsat, Inf,
				0.0, wdisk, tdisk, no_evol);
		Planet planet1(&star_not_saturated_wind_no_coupling, 1.0, 1.0, 1.0);
		StellarSystem system1(&star_not_saturated_wind_no_coupling, &planet1);
		double adestr=planet1.minimum_semimajor(tdisk)*AU_Rsun,
			   tdestr=(beta/2.0*std::pow(adestr, 3) + int_const -
					   Lscale/10.0*std::pow(adestr, 5))/kappa,
			   Lc_at_destr=beta/std::pow(adestr, 1.5) +
				   Lscale*std::sqrt(adestr);
/*		std::cout << std::endl << "alpha=" << Lscale
			<< std::endl << "beta=" << beta
			<< std::endl << "kappa=" << kappa 
			<< std::endl << "wdisk=" << wdisk << std::endl;
		std::cout << std::endl << "ainitial=" << ainitial << std::endl;
		std::cout << std::endl << "Destruction a=" << adestr
			<< ", t=" << tdestr << ", " << ", Lc=" << Lc_at_destr
			<< ", no planet time="
			<< (Lscale*(std::pow(adestr, 5)-std::pow(afinal, 5))/10.0 -
					beta*(std::pow(adestr, 3)-std::pow(afinal, 3))/2.0)/kappa
			<< std::endl;*/

		std::valarray<double> a_transform_coef(0.0, 6), t_coef(2);
		a_transform_coef[5]=Lscale/10.0; a_transform_coef[3]=-beta/2.0;
		t_coef[0]=int_const; t_coef[1]=-kappa;
		std::valarray<double> Lconv_term1_coef(0.0, 2), Lconv_term2_coef(0.0, 3),
			identity_coef(0.0, 2), noplanet_Lconv_m2_coef(2);
		Lconv_term1_coef[1]=1.0/beta*std::pow(10.0/Lscale, 3.0/10.0);
		Lconv_term2_coef[2]=-2.0/std::pow(beta, 3);
		noplanet_Lconv_m2_coef[0]=std::pow(Lc_at_destr, -2) -
			2.0*Kwind*tdestr/std::pow(Ic, 3);
		noplanet_Lconv_m2_coef[1]=2.0*Kwind/std::pow(Ic, 3);
		identity_coef[1]=1;
		PolynomialEvolutionTrack
			identity(identity_coef, -Inf, Inf),
			disk_a_evol(std::valarray<double>(NaN, 2), tstart, tdisk),
			locked_a_transform(a_transform_coef, -Inf, Inf),
			locked_aLconv_evol(t_coef, tdisk, tdestr),
			noplanet_a_evol(std::valarray<double>(NaN, 2), tdestr,
					tfinal),

			disk_Lconv_evol(std::valarray<double>(wdisk*Ic, 1), tstart,
					tdisk),
			noplanet_Lconv_m2_evol(noplanet_Lconv_m2_coef, tdestr, tfinal),
			Lconv_term1_poly(Lconv_term1_coef, -Inf, Inf),
			Lconv_term2_poly(Lconv_term2_coef, -Inf, Inf),
			Lrad_evol(std::valarray<double>(), tstart, tfinal);
		FunctionToPower L_transform1(&Lconv_term1_poly, -10.0/3.0),
						L_transform2(&Lconv_term2_poly, -1.0),
						noplanet_Lconv_evol(&noplanet_Lconv_m2_evol, -0.5);
		LogFunction log_a(&identity);
		FuncPlusFunc locked_Lconv_transform(&L_transform1, &L_transform2);
		PiecewiseFunction a_evol, Lconv_evol;
		a_evol.add_piece(&disk_a_evol);
		a_evol.add_piece(&locked_aLconv_evol);
		a_evol.add_piece(&noplanet_a_evol);
		Lconv_evol.add_piece(&disk_Lconv_evol);
		Lconv_evol.add_piece(&locked_aLconv_evol);
		Lconv_evol.add_piece(&noplanet_Lconv_evol);

		OrbitSolver solver(tstart, tfinal, 1e-8);
		solver(system1, Inf, 0.0, ainitial/AU_Rsun, tstart);
		TransformedSolution to_check(identity, identity, identity, tstart);
		to_check.add_transformation(locked_a_transform,
				locked_Lconv_transform, identity, tdisk);
		to_check.add_transformation(identity, identity, identity, tdestr);
		to_check(solver);

		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_DISK);
		expected_mode.add_break(tdisk, LOCKED_TO_PLANET);
		expected_mode.add_break(tdestr, NO_PLANET);

		test_solution(to_check, a_evol, Lconv_evol, Lrad_evol, tstart, tfinal,
				expected_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}

}

void test_OrbitSolver::test_disklocked_to_fast_to_noplanet()
{
	try {
		const double tdisk=1, tdestr=2, tend=3, wdisk=0.0, Ic=1,
			  Lscale=-AstroConst::jupiter_mass/
				  std::pow(AstroConst::solar_radius, 1.5)*
				  std::sqrt(AstroConst::G/
						  (AstroConst::jupiter_mass+AstroConst::solar_mass))*
				  AstroConst::day;
		MockStellarEvolution no_evol(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(Ic, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		double Q=1e8, alpha=-4.5*std::sqrt(AstroConst::G/
				(AstroConst::solar_radius*AstroConst::solar_mass))*
			AstroConst::jupiter_mass/Q*AstroConst::Gyr/
			AstroConst::solar_radius;
		Star star_no_wind_no_coupling(1.0, Q, 0.0, 1.0, Inf, 0.0, wdisk,
				tdisk, no_evol, 1.0, 0.0, 0.0);
		Planet planet1(&star_no_wind_no_coupling, 1.0, 1.0, 1.0);
		StellarSystem system1(&star_no_wind_no_coupling, &planet1);
		double adestr=planet1.minimum_semimajor(tdisk)*AU_Rsun,
			   a6p5_offset=std::pow(adestr, 6.5)-6.5*alpha*tdestr,
			   a_formation=std::pow(a6p5_offset + 6.5*alpha*tdisk, 1.0/6.5);
		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_DISK);
		expected_mode.add_break(tdisk, FAST_PLANET);
		expected_mode.add_break(tdestr, NO_PLANET);

		std::valarray<double> a6p5_poly_coef(2);
		a6p5_poly_coef[0]=a6p5_offset; a6p5_poly_coef[1]=6.5*alpha;
		PolynomialEvolutionTrack
			a6p5_evol(a6p5_poly_coef, tdisk, tdestr),
			Lconv_disk(std::valarray<double>(Ic*wdisk, 1), tstart, tdisk),
			Lconv_noplanet(std::valarray<double>(
						-Lscale*std::sqrt(a_formation)+Ic*wdisk, 1),
					tdestr, tend),
			a_disk(std::valarray<double>(NaN, 2), tstart, tdisk),
			a_noplanet(std::valarray<double>(NaN, 2), tdestr, tend),
			Lrad_evol(std::valarray<double>(), tstart, tend);
		FunctionToPower a_fast(&a6p5_evol, 1.0/6.5),
						sqrta_evol(&a6p5_evol, 1.0/13.0);
		ExponentialPlusFunc Lconv_unscaled(&sqrta_evol,
				Ic*wdisk/Lscale-std::sqrt(a_formation), 0);
		ScaledFunction Lconv_fast(&Lconv_unscaled, Lscale);

		PiecewiseFunction a_evol, Lconv_evol;
		a_evol.add_piece(&a_disk);
		a_evol.add_piece(&a_fast);
		a_evol.add_piece(&a_noplanet);
		Lconv_evol.add_piece(&Lconv_disk);
		Lconv_evol.add_piece(&Lconv_fast);
		Lconv_evol.add_piece(&Lconv_noplanet);

		OrbitSolver solver(tstart, tend, 5e-8);
		solver(system1, Inf, 0.0, a_formation/AU_Rsun, tstart);
		test_solution(solver, a_evol, Lconv_evol, Lrad_evol, tstart, tend,
				expected_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}
}

void test_OrbitSolver::test_disklocked_to_fast_to_locked()
{
	try {
		const double Q=1e8, alpha=(-4.5*std::sqrt(AstroConst::G/
					(AstroConst::solar_radius*AstroConst::solar_mass))*
				AstroConst::jupiter_mass/Q*AstroConst::Gyr/
				AstroConst::solar_radius),
			  Lscale=AstroConst::jupiter_mass/
				  std::pow(AstroConst::solar_radius, 1.5)*
				  std::sqrt(AstroConst::G/
						  (AstroConst::jupiter_mass+AstroConst::solar_mass))*
				  AstroConst::day,
			  beta=std::sqrt(AstroConst::G*(AstroConst::solar_mass+
						  AstroConst::jupiter_mass))*AstroConst::day/
				  std::pow(AstroConst::solar_radius, 1.5),
			  tdisk=1, async=2.5, tsync=2.0, tend=3,
			  a6p5_offset=std::pow(async, 6.5)-6.5*alpha*tsync,
			  a_formation=std::pow(a6p5_offset + 6.5*alpha*tdisk, 1.0/6.5),
			  Ic=Lscale*(std::sqrt(a_formation)-std::sqrt(async))/
				  (beta*(std::pow(async, -1.5)-0.5*std::pow(a_formation, -1.5))),
			  wdisk=0.5*beta/std::pow(a_formation, 1.5),
			  wlocked=beta/std::pow(async, 1.5);
/*		std::cout << std::endl << "arate=" << -alpha
			<< std::endl << "alpha=" << Lscale
			<< std::endl << "beta=" << beta
			<< std::endl << "Ic=" << Ic
			<< std::endl << "wdisk=" << wdisk
			<< std::endl << "wlocked=" << wlocked
			<< std::endl << "a0=" << a_formation
			<< std::endl;*/
		MockStellarEvolution no_evol(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(Ic, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		Star star_no_wind_no_coupling(1.0, Q, 0.0, 1.0, Inf, 0.0, wdisk,
				tdisk, no_evol, 1.0, 0.0, 0.0);
		Planet planet1(&star_no_wind_no_coupling, 1.0, 0.0, 1.0);
		StellarSystem system1(&star_no_wind_no_coupling, &planet1);
		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_DISK);
		expected_mode.add_break(tdisk, FAST_PLANET);
		expected_mode.add_break(tsync, LOCKED_TO_PLANET);

		std::valarray<double> a6p5_poly_coef(2);
		a6p5_poly_coef[0]=a6p5_offset; a6p5_poly_coef[1]=6.5*alpha;
		PolynomialEvolutionTrack
			a6p5_evol(a6p5_poly_coef, tdisk, tsync),
			Lconv_disk(std::valarray<double>(Ic*wdisk, 1), tstart, tdisk),
			Lconv_locked(std::valarray<double>(Ic*wlocked, 1), tsync, tend),
			a_disk(std::valarray<double>(NaN, 2), tstart, tdisk),
			a_locked(std::valarray<double>(async, 1), tsync, tend),
			Lrad_evol(std::valarray<double>(), tstart, tend);
		FunctionToPower a_fast(&a6p5_evol, 1.0/6.5),
						sqrta_evol(&a6p5_evol, 1.0/13.0);
		ExponentialPlusFunc Lconv_unscaled(&sqrta_evol,
				-Ic*wdisk/Lscale-std::sqrt(a_formation), 0);
		ScaledFunction Lconv_fast(&Lconv_unscaled, -Lscale);

		PiecewiseFunction a_evol, Lconv_evol;
		a_evol.add_piece(&a_disk);
		a_evol.add_piece(&a_fast);
		a_evol.add_piece(&a_locked);
		Lconv_evol.add_piece(&Lconv_disk);
		Lconv_evol.add_piece(&Lconv_fast);
		Lconv_evol.add_piece(&Lconv_locked);

		OrbitSolver solver(tstart, tend, 5e-8);
		solver(system1, Inf, 0.0, a_formation/AU_Rsun, tstart);
		test_solution(solver, a_evol, Lconv_evol, Lrad_evol, tstart, tend,
				expected_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}

}

void test_OrbitSolver::test_disklocked_to_locked_to_fast()
{
	try {
		const double a0=3.2, abreak=3.0, adeath=2.0, tdisk=1, tbreak=2, tdeath=3,
			  tend=4,
			  Rp=adeath*AstroConst::solar_radius/AstroConst::jupiter_radius/
				  (2.44*std::pow(AstroConst::solar_mass/
								 AstroConst::jupiter_mass, 1.0/3.0)),
			  beta=std::sqrt(AstroConst::G*(AstroConst::solar_mass+
						  AstroConst::jupiter_mass))*AstroConst::day/
				  std::pow(AstroConst::solar_radius, 1.5),
			  wdisk=beta/std::pow(a0, 1.5),
			  gamma=(std::pow(abreak, 6.5)-std::pow(adeath, 6.5))/
				  (tdeath-tbreak),
			  Q=9.0*13.0/4.0*std::sqrt(AstroConst::G/AstroConst::solar_mass)*
					  AstroConst::jupiter_mass*AstroConst::Gyr/
					  (gamma*std::pow(AstroConst::solar_radius, 1.5)),
			  alphaL=AstroConst::jupiter_mass/
				  std::pow(AstroConst::solar_radius, 1.5)*
				  std::sqrt(AstroConst::G/
						  (AstroConst::jupiter_mass+AstroConst::solar_mass))*
				  AstroConst::day,		
			  Ic=((std::pow(a0, 5)-std::pow(abreak, 5))*std::pow(abreak, 3.5)*6.5
					  - 5.0*gamma*abreak*abreak)/
				  ((std::pow(a0, 3)-std::pow(abreak, 3))*std::pow(abreak,3.5)*6.5
				   - 3.0*gamma)*
				  alphaL/(5.0*beta),
			  kappa=(std::pow(a0, 5) - std::pow(abreak, 5))*alphaL/10.0 -
				  (std::pow(a0, 3) - std::pow(abreak, 3))*beta*Ic/2.0,
			  Kwind=kappa/(std::pow(AstroConst::G*
						  (AstroConst::solar_mass+AstroConst::jupiter_mass)/
						  std::pow(AstroConst::solar_radius, 3), 1.5)*
					  std::pow(AstroConst::day, 3)),
			  locked_a_int_const=alphaL*std::pow(a0, 5)/10.0 -
				  beta*Ic*std::pow(a0, 3)/2.0 + kappa*tdisk;
/*		std::cout << std::endl << "a0=" << a0
			<< std::endl << "alphaL=" << alphaL
			<< std::endl << "beta=" << beta
			<< std::endl << "gamma=" << gamma
			<< std::endl << "Ic=" << Ic
			<< std::endl << "wdisk=" << wdisk
			<< std::endl << "kappa=" << kappa
			<< std::endl << "Kwind=" << Kwind
			<< std::endl << "Q=" << Q
			<< std::endl << "Rp=" << Rp
			<< std::endl;*/
		MockStellarEvolution no_evol(-1,
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(Ic, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1),
				std::valarray< std::valarray<double> >(
					std::valarray<double>(1.0, 1), 1));
		Star star_no_coupling(1.0, Q, Kwind, Inf, Inf, 0.0, wdisk,
				tdisk, no_evol);
		Planet planet1(&star_no_coupling, 1.0, Rp, 1.0);
		StellarSystem system1(&star_no_coupling, &planet1);
		ExpectedEvolutionMode expected_mode;
		expected_mode.add_break(tstart, LOCKED_TO_DISK);
		expected_mode.add_break(tdisk, LOCKED_TO_PLANET);
		expected_mode.add_break(tbreak, FAST_PLANET);
		expected_mode.add_break(tdeath, NO_PLANET);

		std::valarray<double> a_locked_transform_coef(0.0, 6),
			a_locked_evol_coef(2), identity_coef(0.0, 2),
			a6p5_fast_evol_coef(2), Lconv_locked_term1_coef(0.0, 2),
			Lconv_locked_term2_coef(0.0, 3);
		a_locked_transform_coef[5]=alphaL/10.0;
		a_locked_transform_coef[3]=-beta*Ic/2.0;
		a_locked_evol_coef[0]=locked_a_int_const;
		a_locked_evol_coef[1]=-kappa;
		identity_coef[1]=1.0;
		a6p5_fast_evol_coef[0]=gamma*tbreak+std::pow(abreak, 6.5);
		a6p5_fast_evol_coef[1]=-gamma;
		Lconv_locked_term1_coef[1]=1.0/(beta*Ic)*std::pow(10.0/alphaL, 3.0/10.0);
		Lconv_locked_term2_coef[2]=-2.0/std::pow(beta*Ic, 3);

		PolynomialEvolutionTrack
			identity(identity_coef, -Inf, Inf),
			a_disk_transform=identity,
			a_disk_evol(std::valarray<double>(NaN, 2), tstart, tdisk),
			a_locked_transform(a_locked_transform_coef, -Inf, Inf),
			a_locked_evol(a_locked_evol_coef, tdisk, tbreak),
			a_fast_transform=identity,
			a6p5_fast_evol(a6p5_fast_evol_coef, tbreak, tdeath),
			a_noplanet_transform=identity,
			a_noplanet_evol(std::valarray<double>(NaN, 2), tdeath, tend),

			Lconv_disk_transform=identity,
			Lconv_disk_evol(std::valarray<double>(wdisk*Ic, 1), tstart,
					tdisk),
			Lconv_locked_term1_poly(Lconv_locked_term1_coef, -Inf, Inf),
			Lconv_locked_term2_poly(Lconv_locked_term2_coef, -Inf, Inf),
			Lconv_locked_evol=a_locked_evol,
			Lconv_fast_transform(std::valarray<double>(NaN, 2), -Inf, Inf),
			Lconv_fast_evol(std::valarray<double>(NaN, 2), tbreak, tdeath),
			Lconv_noplanet_transform(std::valarray<double>(NaN,2), -Inf,Inf),
			Lconv_noplanet_evol(std::valarray<double>(NaN,2), tdeath, tend),

			Lrad_transform=identity,
			Lrad_evol(std::valarray<double>(), tstart, tend);

		FunctionToPower
			a_fast_evol(&a6p5_fast_evol, 1.0/6.5),
			Lconv_locked_transform1(&Lconv_locked_term1_poly, -10.0/3.0),
			Lconv_locked_transform2(&Lconv_locked_term2_poly, -1.0);

		FuncPlusFunc Lconv_locked_transform(&Lconv_locked_transform1,
				&Lconv_locked_transform2);

		PiecewiseFunction a_evol, Lconv_evol;
		a_evol.add_piece(&a_disk_evol);
		a_evol.add_piece(&a_locked_evol);
		a_evol.add_piece(&a_fast_evol);
		a_evol.add_piece(&a_noplanet_evol);
		Lconv_evol.add_piece(&Lconv_disk_evol);
		Lconv_evol.add_piece(&Lconv_locked_evol);
		Lconv_evol.add_piece(&Lconv_fast_evol);
		Lconv_evol.add_piece(&Lconv_noplanet_evol);

		OrbitSolver solver(tstart, tend, 1e-9);
		solver(system1, Inf, 0.0, a0/AU_Rsun, tstart);

		TransformedSolution to_check(a_disk_transform, Lconv_disk_transform,
				Lrad_transform, tstart);
		to_check.add_transformation(a_locked_transform,
				Lconv_locked_transform, Lrad_transform, tdisk);
		to_check.add_transformation(a_fast_transform, Lconv_fast_transform,
				Lrad_transform, tbreak);
		to_check.add_transformation(a_noplanet_transform,
				Lconv_noplanet_transform, Lrad_transform, tdeath);
		to_check(solver);

		test_solution(to_check, a_evol, Lconv_evol, Lrad_evol, tstart, tend,
				expected_mode);
	} catch (Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()+": "+ex.get_message()).c_str());
	} catch (std::exception &ex) {
		TEST_ASSERT_MSG(false, (std::string("Unexpected exception thrown: ")+
				ex.what()).c_str());
	}

}


test_OrbitSolver::test_OrbitSolver()
{
	TEST_ADD(test_OrbitSolver::test_disk_locked_no_stellar_evolution);
	TEST_ADD(test_OrbitSolver::test_disk_locked_with_stellar_evolution);
	TEST_ADD(test_OrbitSolver::test_no_planet_evolution);
	TEST_ADD(test_OrbitSolver::test_unlocked_evolution); 
	TEST_ADD(test_OrbitSolver::test_locked_evolution);
	TEST_ADD(test_OrbitSolver::test_disklocked_to_locked_to_noplanet);
	TEST_ADD(test_OrbitSolver::test_disklocked_to_fast_to_noplanet);
	TEST_ADD(test_OrbitSolver::test_disklocked_to_fast_to_locked); 
	TEST_ADD(test_OrbitSolver::test_disklocked_to_locked_to_fast);
}

#ifdef STANDALONE
int main()
{
	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_OrbitSolver tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
