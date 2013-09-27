/**\file
 *
 * \brief Produces an executable that calculates orbital evolutions of
 * planet-star systems.
 *
 */

#include "poet.h"

void CommandLineOptions::define_options()
{
	__low_mass_windK=arg_dbl0("K", "low-mass-K", "<double>",
			"The wind strength for low mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). Default: 0.35");

	__high_mass_windK=arg_dbl0(NULL, "high-mass-K", "<double>",
			"The wind strength for low mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). Default: 0");

	__low_mass_wind_saturation=arg_dbl0(NULL, "low-mass-wind-sat",
			"<double>", "The frequency at which the wind saturates for low "
			"mass stars in units of rad/day. Default: 1.84.");

	__high_mass_wind_saturation=arg_dbl0(NULL, "high-mass-wind-sat",
			"<double>", "The frequency at which the wind saturates for high "
			"mass stars in units of rad/day. Default: 0.");

	__core_env_coupling_timescale=arg_dbl0(NULL,
			"core-env-coupling-timescale", "<double>", "The timescale on "
			"which the core end envelope are coupled in Myr. Default: 5.");

	__lgQ=arg_dbl0(NULL, "lgQ", "<double>", "Log base 10 of the tidal "
			"quality factor of the star. Default: 6.");

	__star_mass=arg_dbl0("M", "Mstar", "<double>", "Mass of the star in "
			"solar masses. Default: 1");

	__planet_mass=arg_dbl0("m", "Mplanet", "<double>", "Mass of the planet "
			"in jupiter masses. Default: 1");

	__planet_radius=arg_dbl0("r", "Rplanet", "<double>", "Radius of the "
			"planet in jupiter radii. Default: 1");

	__planet_formation_age=arg_dbl0(NULL, "planet-formation-age", "<double>",
			"The age at which the planet forms. If it is smaller than the "
			"disk dissipation age, the planet forms at the disk dissipation "
			"age. Default 0.");

	__disk_lock_frequency=arg_dbl0(NULL, "w-disk", "<double>", "The spin at "
			"which the star is locked while the disk is present in rad/day. "
			"Default: 0.68");

	__disk_dissipation_age=arg_dbl0(NULL, "t-disk", "<double>", "The age at "
			"which the disk dissipates in Myr. Default: 5.");

	__planet_formation_semimajor=arg_dbl0("a", "init-semimajor", "<double>",
			"The semimajor axis at which the planet first appears in AU. "
			"Default: 0.05");

	__start_age=arg_dbl0(NULL, "t0", "<double>", "The starting age for the "
			"evolution in Gyr. If this argument is used, --w-disk, "
			"--w-disk-col, --t-disk and --t-disk-col are ignored and "
			"--init-semimajor or --init-semimajor-col is the semimajor axis "
			"that the planet has at the age specified by this argument. "
			"Further, if this argument is  used, --init-wsurf or "
			"--init-wsurf-col must also be specified and if the evolution "
			"of low mass stars is to be computed --init-wrad or "
			"--init-wrad-col is also required.");

	__end_age=arg_dbl0(NULL, "tmax", "<double>", "The maximum end age for "
			"the evolution in Gyr. The evolution could stop earlier if the "
			"star moves off the main sequence before this age. Default: "
			"Inf.");

	__start_wrad=arg_dbl0(NULL, "init-wrad", "<double>", "Initial spin of "
			"the stellar core in rad/day for low mass stars. This argument "
			"is ignored, unless --t0 or --t0-col is also specified.");

	__start_wsurf=arg_dbl0(NULL, "init-wsurf", "<double>", "Initial spin of "
			"the stellar surface in rad/day. For low mass stars this is the "
			"rotation of the convective zone, and for high mass stars it is "
			"the unique rotation of the star. This argument is ignored, "
			"unless --t0 or --t0-col is also specified.");

	__max_timestep=arg_dbl0(NULL, "max-step", "<double>", "An upper limit to"
			" impose on the sover timestep.");

	__precision=arg_dbl0(NULL, "precision", "<double>", "The number of "
			"digits of precision to require of the solution. Need not be an "
			"integer. Default: 6.");

	__low_mass_windK_column=arg_int0(NULL, "low-mass-K-col", "<index>",
			"Index in the input file of the column which contains the wind "
			"strength for low mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). If this parameter is not "
			"specified the const value specified by --low-mass-K or its "
			"default is used.");

	__high_mass_windK_column=arg_int0(NULL, "high-mass-K-col", "<index>",
			"Index in the input file of the column which contains the wind "
			"strength for high mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). If this parameter is not "
			"specified the const value specified by --high-mass-K or its "
			"default is used.");

	__low_mass_wind_saturation_column=arg_int0(NULL, "low-mass-wind-sat-col",
			"<index>", "Index in the input file of the column which "
			"contains the frequency at which the wind saturates for low "
			"mass stars in units of rad/day. If this parameter is not "
			"specified the const value specified by --low-mass-wind-sat or "
			"its default is used.");

	__high_mass_wind_saturation_column=arg_int0(NULL,
			"high-mass-wind-sat-col", "<index>",
			"Index in the input file of the column which contains the "
			"frequency at which the wind saturates for high mass stars in "
			"units of rad/day. If this parameter is not specified the const "
			"value specified by --high-mass-wind-sat or its default is "
			"used.");

	__core_env_coupling_timescale_column=arg_int0(NULL,
			"core-env-coupling-timescale-col", "<index>", "Index in the "
			"input file of the column which contains the timescale on "
			"which the core end envelope are coupled in Myr. If this "
			"parameter is not specified the const value specified by "
			"--core-env-coupling-timescale or its default is used.");

	__lgQ_column=arg_int0(NULL, "lgQ-col", "<index>", "Index in the "
			"input file of the column which contains log base 10 of the "
			"tidal quality factor of the star. If this parameter is not "
			"specified the const value specified by --lgQ or its default "
			"is used.");

	__star_mass_column=arg_int0(NULL, "Mstar-col", "<index>", "Index in the "
			"input file of the column which contains the mass of the star in "
			"solar masses. If this parameter is not specified the const "
			"value specified by --Mstar or its default is used.");

	__planet_mass_column=arg_int0(NULL, "Mplanet-col", "<index>", "Index in "
			"the input file of the column which contains the mass of the "
			"planet in jovian masses. If this parameter is not specified "
			"the const value specified by --Mplanet or its default is "
			"used.");

	__planet_radius_column=arg_int0(NULL, "Rplanet-col", "<index>", "Index "
			"in the input file of the column which contains the radius of "
			"the planet in jovian radii. If this parameter is not specified "
			"the const value specified by --Rplanet or its default is "
			"used.");

	__planet_formation_age_column=arg_int0(NULL, "planet-formation-age-col",
			"<index>", "Index in the input file of the column which "
			"contains the age at which the planet forms. If that age is "
			" smaller than the disk dissipation age, the planet forms at the"
			" disk dissipation age.");

	__disk_lock_frequency_column=arg_int0(NULL, "w-disk-col", "<index>",
			"Index in the input file of the column which contains the spin "
			"at which the star is locked while the disk is present in "
			"rad/day. If this parameter is not specified the const "
			"value specified by --w-disk or its default is used.");

	__disk_dissipation_age_column=arg_int0(NULL, "t-disk-col", "<index>",
			"Index in the input file of the column which contains the age "
			"at which the disk dissipates in Myr. If this parameter is not "
			"specified the const value specified by --t-disk or its default "
			"is used.");

	__planet_formation_semimajor_column=arg_int0(NULL, "init-semimajor-col",
			"<index>", "Index in the input file of the column which contains"
			" the semimajor axis at which the planet first appears in AU. If"
			" this parameter is not specified the const value specified by "
			"--init-semimajor or its default is used.");

	__start_age_column=arg_int0(NULL, "t0-col", "<index>", "Index in the "
			"input file of the column which contains the starting age for "
			"the evolution in Gyr. If this argument is used, --w-disk, "
			"--w-disk-col, --t-disk and --t-disk-col are ignored and "
			"--init-semimajor or --init-semimajor-col is the semimajor axis "
			"that the planet has at the age specified by this argument. "
			"Further, if this argument is  used, --init-wsurf or "
			"--init-wsurf-col must also be specified and if the evolution "
			"of low mass stars is to be computed --init-wrad or "
			"--init-wrad-col is also required.");

	__end_age_column=arg_int0(NULL, "tmax-col", "<index>", "Index in the "
			"input file of the column which contains the maximum end age for"
			" the evolution in Gyr. If this parameter is not specified the "
			"const value specified by --tmax or its default is used.");

	__start_wrad_column=arg_int0(NULL, "init-wrad-col", "<index>", "Index "
			"in the input file of the column which contains the initial spin"
			" of the stellar core in rad/day for low mass stars. This "
			"argument is ignored, unless --t0 or --t0-col is also "
			"specified.");

	__start_wsurf_column=arg_int0(NULL, "init-wsurf-col", "<index>",
			"Index in the input file of the column which contains the "
			"initial spin of the stellar surface in rad/day. For low mass "
			"stars this is the rotation of the convective zone, and for "
			"high mass stars it is the unique rotation of the star. This "
			"argument is ignored, unless --t0 or --t0-col is also "
			"specified.");

	__input_fname=arg_file0("i", "input", "<file>", "The file to read the "
			"parameters of the planet-star systems to calculate evolutions "
			"for. If omitted, standard input is used instead. If no --*-col "
			"options are specified, this argument is ignored.");

	__output_fname=arg_file0("o", "output", "<file|index>", "If nothing is "
			"read from an input file, then this argument should be the name "
			"of the file to use to output the evolution. If at least one "
			"quantity is read from a list file, that file should also "
			"contain a column specifying the input file name for each "
			"parameter set, and this argument should be the index of that "
			"column. Default: poet.evol");

	__serialized_stellar_evolution=arg_file0(NULL, "serialized-stellar-evol",
			"<file>", "The file to read previously serialized stellar "
			"evolution from. Default: 'interp_state_data'.");

	__start_locked=arg_lit0(NULL, "start-locked", "Whether the planet "
			"should start locked to the star. If true, causes the value of "
			"--init-wsurf-col to be ignored.");
}

void CommandLineOptions::set_defaults()
{
	__low_mass_windK->dval[0]=0.35;
	__high_mass_windK->dval[0]=0;
	__low_mass_wind_saturation->dval[0]=1.84;
	__high_mass_wind_saturation->dval[0]=0;
	__core_env_coupling_timescale->dval[0]=5;
	__lgQ->dval[0]=6;
	__star_mass->dval[0]=1;
	__planet_mass->dval[0]=1;
	__planet_radius->dval[0]=1;
	__planet_formation_age->dval[0]=0;
	__disk_lock_frequency->dval[0]=0.68;
	__disk_dissipation_age->dval[0]=5;
	__planet_formation_semimajor->dval[0]=0.05;
	__start_age->dval[0]=NaN; 
	__end_age->dval[0]=Inf; 
	__start_wrad->dval[0]=NaN;
	__start_wsurf->dval[0]=NaN;
	__max_timestep->dval[0]=Inf;
	__precision->dval[0]=6;
	__serialized_stellar_evolution->filename[0]="interp_state_data";

	__low_mass_windK_column->ival[0]=0;
	__high_mass_windK_column->ival[0]=0;
	__low_mass_wind_saturation_column->ival[0]=0;
	__high_mass_wind_saturation_column->ival[0]=0;
	__core_env_coupling_timescale_column->ival[0]=0;
	__lgQ_column->ival[0]=0;
	__star_mass_column->ival[0]=0;
	__planet_mass_column->ival[0]=0;
	__planet_radius_column->ival[0]=0;
	__planet_formation_age_column->ival[0]=0;
	__disk_lock_frequency_column->ival[0]=0;
	__disk_dissipation_age_column->ival[0]=0;
	__planet_formation_semimajor_column->ival[0]=0;
	__start_age_column->ival[0]=0; 
	__start_wrad_column->ival[0]=0;
	__start_wsurf_column->ival[0]=0;

	__output_fname->filename[0]="poet.evol";
}

void CommandLineOptions::postprocess()
{
	if(--__low_mass_windK_column->ival[0]>=0 || 
			--__high_mass_windK_column->ival[0]>=0 ||
			--__low_mass_wind_saturation_column->ival[0]>=0 ||
			--__high_mass_wind_saturation_column->ival[0]>=0 ||
			--__core_env_coupling_timescale_column->ival[0]>=0 ||
			--__lgQ_column->ival[0]>=0 ||
			--__star_mass_column->ival[0]>=0 ||
			--__planet_mass_column->ival[0]>=0 ||
			--__planet_radius_column->ival[0]>=0 ||
			--__planet_formation_age_column->ival[0]>=0 ||
			--__disk_lock_frequency_column->ival[0]>=0 ||
			--__disk_dissipation_age_column->ival[0]>=0 ||
			--__planet_formation_semimajor_column->ival[0]>=0 ||
			--__start_age_column->ival[0]>=0 ||
			--__start_wrad_column->ival[0]>=0 ||
			--__start_wsurf_column->ival[0]>=0) {
		__input_from_list=true;
		std::istringstream outfname_col_str(__output_fname->filename[0]);
		outfname_col_str >> __output_fname_column;
	} else __input_from_list=false;
	if(__input_fname->count) __input_stream.open(__input_fname->filename[0]);
	__precision->dval[0]=std::pow(10.0, -__precision->dval[0]);
	__core_env_coupling_timescale->dval[0]*=1e-3;
	__disk_dissipation_age->dval[0]*=1e-3;
}

void CommandLineOptions::cleanup()
{
	if(__input_fname->count) __input_stream.close();
	free(__low_mass_windK);
	free(__low_mass_windK_column);
	free(__high_mass_windK);
	free(__high_mass_windK_column);
	free(__low_mass_wind_saturation);
	free(__low_mass_wind_saturation_column);
	free(__high_mass_wind_saturation);
	free(__high_mass_wind_saturation_column);
	free(__core_env_coupling_timescale);
	free(__core_env_coupling_timescale_column);
	free(__lgQ);
	free(__lgQ_column);
	free(__star_mass);
	free(__star_mass_column);
	free(__planet_mass);
	free(__planet_mass_column);
	free(__planet_radius);
	free(__planet_radius_column);
	free(__planet_formation_age);
	free(__planet_formation_age_column);
	free(__disk_lock_frequency);
	free(__disk_lock_frequency_column);
	free(__disk_dissipation_age);
	free(__disk_dissipation_age_column);
	free(__planet_formation_semimajor);
	free(__planet_formation_semimajor_column);
	free(__start_age);
	free(__start_age_column);
	free(__end_age);
	free(__end_age_column);
	free(__start_wrad);
	free(__start_wrad_column);
	free(__start_wsurf);
	free(__start_wsurf_column);
	free(__max_timestep);
	free(__precision);
	free(__serialized_stellar_evolution);
	free(__input_fname);
	free(__output_fname);
	free(__start_locked);
}

CommandLineOptions::CommandLineOptions(int argc, char **argv)
{
	define_options();
	arg_lit *help_option=arg_lit0("h", "help", "Print this help and exit.");
	struct arg_end *end = arg_end(100);
	void *argtable[] = {
		__low_mass_windK, 				__low_mass_windK_column,
		__high_mass_windK, 				__high_mass_windK_column,
		__low_mass_wind_saturation, 	__low_mass_wind_saturation_column,
		__high_mass_wind_saturation, 	__high_mass_wind_saturation_column,
		__core_env_coupling_timescale, 	__core_env_coupling_timescale_column,
		__lgQ,							__lgQ_column,
		__star_mass,					__star_mass_column,
		__planet_mass,					__planet_mass_column,
		__planet_radius,				__planet_radius_column,
		__planet_formation_age,			__planet_formation_age_column,
		__disk_lock_frequency,			__disk_lock_frequency_column,
		__disk_dissipation_age,			__disk_dissipation_age_column,
		__planet_formation_semimajor,	__planet_formation_semimajor_column,
		__start_age,					__start_age_column,
		__end_age,						__end_age_column,
		__start_wrad,					__start_wrad_column,
		__start_wsurf,					__start_wsurf_column,
		__max_timestep,
		__precision,
		__serialized_stellar_evolution,
		__start_locked,
		__input_fname,
		__output_fname,
		help_option, end};
	if(arg_nullcheck(argtable) != 0) {
		cleanup();
		throw Error::CommandLine("Failed to allocate argument table.");
	}
	set_defaults();
	int nerrors=arg_parse(argc, argv, argtable);
	if(help_option->count>0 || nerrors>0) {
        printf("Usage: %s", "SubPixPhot");
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		if(help_option->count==0)
			arg_print_errors(stdout, end, "poet");
		free(help_option);
		free(end);
		return;
	}
	postprocess();
	free(help_option);
	free(end);
}

///AU/\f$\mathrm{R}_\odot\f$.
const double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;

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

void calculate_evolution(double Mstar, double Q, double Kwind, double wsat, 
		double coupling_timescale, double wdisk, double tdisk, 
		double Mplanet, double Rplanet, double planet_formation_age,
		double a_formation, double start_wrad, double start_wsurf,
		double tstart, double tend, double max_time_step, bool start_locked,
		const StellarEvolution &stellar_evolution, OrbitSolver &solver,
		const std::string &outfname)
{
	Star star(Mstar, Q, Kwind, wsat, coupling_timescale, 0.0, wdisk, tdisk,
			stellar_evolution);
	if(star.is_low_mass() && 
			star.core_formation_age()>star.get_disk_dissipation_age() &&
			tstart<star.core_formation_age())
		throw Error::Runtime("At present the case when the disk dissipates "
				"before the stellar core starts to form is not supported.");
	Planet planet(&star, Mplanet, Rplanet, a_formation);
	StellarSystem system(&star, &planet);
	std::valarray<double> start_orbit(0.0, 1);
	EvolModeType start_evol_mode;
	if(std::isnan(tstart) || tstart<star.get_disk_dissipation_age())
		start_evol_mode=LOCKED_TO_DISK;
	if(std::isnan(tstart) ||
			(star.is_low_mass() && tstart<star.core_formation_age()))
		tstart=NaN;
	else if(tstart<star.get_disk_dissipation_age()) {
		start_orbit[0]=(star.is_low_mass() ? start_wrad : start_wsurf)*
			star.moment_of_inertia(tstart, radiative);
	} else if(start_locked) {
		start_evol_mode=LOCKED_TO_PLANET;
		start_orbit.resize(2);
		start_orbit[0]=a_formation*AU_Rsun;
		start_orbit[1]=(star.is_low_mass() ? start_wrad*
				star.moment_of_inertia(tstart, radiative) : 0);
	} else {
		if(planet.orbital_angular_velocity_semimajor(a_formation)<
				start_wsurf) start_evol_mode=SLOW_PLANET;
		else start_evol_mode=FAST_PLANET;
		start_orbit.resize(3);
		start_orbit[0]=std::pow(a_formation*AU_Rsun, 6.5);
		start_orbit[1]=start_wsurf*
			star.moment_of_inertia(tstart, convective);
		start_orbit[2]=(star.is_low_mass() ? start_wrad*
				star.moment_of_inertia(tstart, radiative) : 0);
	}
	solver(system, max_time_step, planet_formation_age, a_formation, tstart,
			start_evol_mode, start_orbit);
	output_solution(solver, system, outfname);
}

///Calculates a realistic evolution chosen to be comlicated.
int main(int argc, char **argv)
{
	try {
		CommandLineOptions options(argc, argv);
		YRECEvolution stellar_evolution;
		stellar_evolution.load_state(
				options.serialized_stellar_evolution());
		OrbitSolver solver(options.start_age(), options.end_age(),
				options.precision());
		double Kwind, wsat;
		if(options.star_mass()>stellar_evolution.get_mass_break()) {
			Kwind=options.high_mass_windK();
			wsat=options.high_mass_wind_saturation();
		} else {
			Kwind=options.low_mass_windK();
			wsat=options.low_mass_wind_saturation();
		}
		calculate_evolution(options.star_mass(),
				std::pow(10.0, options.lgQ()),
				Kwind, wsat, options.core_env_coupling_timescale(),
				options.disk_lock_frequency(),
				options.disk_dissipation_age(), options.planet_mass(),
				options.planet_radius(), options.planet_formation_age(),
				options.planet_formation_semimajor(), options.start_wrad(),
				options.start_wsurf(), options.start_age(),
				options.end_age(), options.max_timestep(),
				options.start_locked(), stellar_evolution, solver,
				options.output_filename());
	} catch(Error::General &err) {
		std::cerr << err.what() << ": " << err.get_message() << std::endl;
		return 1;
	}
}
