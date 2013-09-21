/**\file
 *
 * \brief Produces an executable that calculates and outputs a relatively
 * complicated evolution of a planet around a star.
 *
 * \ingroup UnitTests_group
 */

#include "example_evolutions.h"

///AU/\f$\mathrm{R}_\odot\f$.
const double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;

///Outputs the solution calculated by the given solver.
void output_solution(const OrbitSolver &solver, const StellarSystem &system,
		const std::string &filename)
{
	std::list<double>::const_iterator
		age_i=solver.get_tabulated_var(AGE)->begin(),
		a_i=solver.get_tabulated_var(SEMIMAJOR)->begin(),
		Lconv_i=solver.get_tabulated_var(LCONV)->begin(),
		Lrad_i=solver.get_tabulated_var(LRAD)->begin();

	std::list<double>::const_iterator
		last_age=solver.get_tabulated_var(AGE)->end();

	std::list<EvolModeType>::const_iterator
		mode_i=solver.get_tabulated_evolution_mode()->begin();
	std::ofstream outf(filename.c_str());
	outf.precision(16);
	outf << '#' << std::setw(24) << "age"
		<< std::setw(25) << "semimajor"
		<< std::setw(25) << "Lconv"
		<< std::setw(25) << "Lrad"
		<< std::setw(25) << "Iconv"
		<< std::setw(25) << "Irad"
		<< std::setw(25) << "mode"
		<< std::endl;
	while(age_i!=last_age) {
		outf << std::setw(25) << *age_i
			<< std::setw(25) << *a_i
			<< std::setw(25) << *Lconv_i
			<< std::setw(25) << *Lrad_i
			<< std::setw(25)
			<< system.get_star().moment_of_inertia(*age_i, convective)
			<< std::setw(25)
			<< system.get_star().moment_of_inertia(*age_i, radiative)
			<< std::setw(25) << *mode_i
			<< std::endl;
		age_i++; a_i++; Lconv_i++; Lrad_i++; mode_i++;
	}
	outf.close();
}

///\brief Sets up a test case for which an exact analytical solution is known
///and calculates it.
void calculate_test()
{
	const double tstart=2.0*min_age, Q=1e8,
		  alpha=(-4.5*std::sqrt(AstroConst::G/
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
		  wdisk=0.5*beta/std::pow(a_formation, 1.5);

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
	Planet planet(&star_no_wind_no_coupling, 1.0, 0.0, 1.0);
	StellarSystem system(&star_no_wind_no_coupling, &planet);
	OrbitSolver solver(tstart, tend, 5e-8);
	solver(system, Inf, 0.0, a_formation/AU_Rsun, tstart);
	output_solution(solver, system, "DiskFastLocked_test.txt");
}

///Calculates a realistic evolution chosen to be comlicated (no analytical
///solution is available).
void calculate_full()
{
	YRECEvolution stellar_evolution;
	stellar_evolution.load_state("../interp_state_data_phs4");
	const double Mstar=1, Q=3e6, Kwind=0.17, wsat=2.2,
		  coupling_timescale=0.03, wdisk=0.9, tdisk=5e-3, Mplanet=10,
		  Rplanet=1, a_formation=10.0, tstart=1e-3, tend=10.0;
	Star star(Mstar, Q, Kwind, wsat, coupling_timescale, 0.0, wdisk, tdisk,
			stellar_evolution);
	Planet planet(&star, Mplanet, Rplanet, a_formation/AU_Rsun);
	StellarSystem system(&star, &planet);
	OrbitSolver solver(tstart, tend, 5e-6);
	solver(system, Inf, 0.0, a_formation/AU_Rsun, tstart);
	output_solution(solver, system, "FullEvolution.txt");
//	solver(system, 1e-3, 0.0, a_formation/AU_Rsun, tstart);
//	output_solution(solver, system, "FullEvolutionHD.txt");
}

void calculate_slow()
{
	YRECEvolution stellar_evolution;
	stellar_evolution.load_state("../interp_state_data_phs4");
	const double Mstar=0.90000000000000013323,
		  Q=1e6,
		  Kwind=0.155,
		  wsat=2.454,
		  coupling_timescale=0.012,
		  wdisk=2*M_PI/1.4,
		  tdisk=2.5e-3,
		  Mplanet=25,
		  Rplanet=0.714,
		  P0=5.9000000000000003553,
		  a_formation=AstroConst::G*Mstar*AstroConst::solar_mass*
			  std::pow(P0*AstroConst::day, 2)/4/M_PI/M_PI,
		  tstart=MIN_AGE;
	Star star(Mstar, Q, Kwind, wsat, coupling_timescale, 0.0, wdisk, tdisk,
			stellar_evolution);
	double tend=std::min((const double)MAX_END_AGE,
				star.get_lifetime());/**/
	Planet planet(&star, Mplanet, Rplanet, a_formation/AU_Rsun);
	StellarSystem system(&star, &planet);
	OrbitSolver solver(tstart, tend, 1e-5, SPIN_THRES, MAIN_SEQ_START);
	solver(system, Inf, PLANET_FORM_AGE, a_formation/AstroConst::AU, tstart);
	output_solution(solver, system, "SlowEvolution.txt");
}

///Calculates a realistic evolution chosen to be comlicated.
int main()
{
//	try {
//		calculate_test();
//		calculate_full();
		calculate_slow();
/*	} catch(Error::General &ex) {
		std::cerr << "Unexpected exception thrown: " << ex.what() << ":	"
			<< ex.get_message() << std::endl;
	}*/
}
