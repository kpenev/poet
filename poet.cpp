/**\file
 *
 * \brief Produces an executable that calculates orbital evolutions of
 * planet-star systems.
 *
 */

#include "poet.h"

const std::string CommandLineOptions::__input_column_names[]={
	"K", "K", "wsat", "wsat", "psat", "psat", "tcoup", "lgQ", "M", "m", "r",
	"tform", "wdisk", "pdisk", "tdisk", "aform", "pform", "t0", "tmax",
	"wrad0", "wsurf0", "locked", "maxdt", "prec", "outf"};

const std::string CommandLineOptions::__output_column_names[]={
	"t", "a", "Lconv", "Lrad", "L", "Iconv", "Irad", "I", "Wsurf", "Wrad",
	"Psurf", "Prad", "mode", "R", "Lum", "Rrad", "Mrad"};

const std::string CommandLineOptions::__output_column_descr[]={
	"Age in Gyr",
	"Semimajor axis in AU",
	"Convective zone angular momentum (low mass stars only else NaN).",
	"Radiative zone angular momentum (low mass stars only else NaN).",
	"Total angular momentum of the star.",
	"Convective zone moment of inertia (low mass stars only else NaN).",
	"Radiative zone moment of inertia (low mass stars only else NaN).",
	"Total moment of inertia of the star.",
	"Stellar surface angular velocity of the star in rad/day.",
	"Angular velocity of the radiative zone in rad/day (low mass stars only,"
		", else NaN)",
	"Stellar surface spin period in days.",
	"Spin period of the radiative zone in days (low mass stars only, else "
		"NaN)",
	"Evolution mode of the step starting at this age.",
	"Stellar radius",
	"Stellar luminosity",
	"Radius of the stellar radiative zone (low mass stars only, else NaN).",
	"Mass of the stellar radiative zone (low mass stars only, else NaN)."};

const double CommandLineOptions::__defaults[]={
		0.35,			//LOW_MASS_WINDK=WINDK
		0,				//HIGH_MASS_WINDK=SKIP
		1.84,			//LOW_MASS_WIND_SAT_W=WIND_SAT_W
		0,				//HIGH_MASS_WIND_SAT_W,
		NaN,			//LOW_MASS_WIND_SAT_P=WIND_SAT_P, 
		Inf,			//HIGH_MASS_WIND_SAT_P,
		5,				//CORE_ENV_COUPLING_TIMESCALE,
		6,				//LGQ
		1,				//MSTAR
		1,				//MPLANET
		1,				//RPLANET
		0,				//PLANET_FORMATION_AGE
		0.68,			//WDISK
		NaN,			//PDISK
		5,				//TDISK
		0.05,			//A_FORMATION
		NaN,			//P_FORMATION
		NaN,			//TSTART
		Inf,			//TEND
		NaN,			//START_WRAD, 
		NaN,			//START_WSURF, 
		NaN,			//START_LOCKED,	
		Inf,			//MAX_STEP
		6				//PRECISION
};

const std::string CommandLineOptions::__default_outfname="poet.evol",
	  CommandLineOptions::__default_serialized_evol="interp_state_data",
	  CommandLineOptions::__default_output_columns=
	  	"t,a,Lconv,Lrad,L,Iconv,Irad,I,mode";

void CommandLineOptions::define_options()
{
	std::ostringstream option_help;
	option_help << "The wind strength for low mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). To identify this quantity in "
			"--input-columns use '"
			<< __input_column_names[InCol::LOW_MASS_WINDK] << "' (identical "
			"to the --high-mass-K since only one K can be applicable to a "
			"run). Default: " << __defaults[InCol::LOW_MASS_WINDK] << ".";
	__direct_value_options[InCol::LOW_MASS_WINDK]=arg_dbl0("K", "low-mass-K",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The wind strength for low mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). To identify this quantity in "
			"--input-columns use '"
			<<  __input_column_names[InCol::HIGH_MASS_WINDK]
			<< "' (identical to the --low-mass-K since only one K can be "
			"applicable to a run). Default: "
			<< __defaults[InCol::LOW_MASS_WINDK] << ".";
	__direct_value_options[InCol::HIGH_MASS_WINDK]=arg_dbl0(NULL, "high-mass-K",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The frequency at which the wind saturates for low mass "
		"stars in units of rad/day. To identify this quantity in "
		"--input-columns use '"
		<<  __input_column_names[InCol::LOW_MASS_WIND_SAT_W]
		<< "' (identical to the --high-mass-wind-sat-w since only one "
		"saturation frequency can be applicable to a run). Default: "
		<< __defaults[InCol::LOW_MASS_WIND_SAT_W] << ".";
	__direct_value_options[InCol::LOW_MASS_WIND_SAT_W]=arg_dbl0(NULL,
			"low-mass-wind-sat-w", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The frequency at which the wind saturates for high mass "
		"stars in units of rad/day. To identify this quantity in "
		"--input-columns use '"
		<< __input_column_names[InCol::HIGH_MASS_WIND_SAT_W]
		<< "' (identical to the --high-mass-wind-sat-w since only one "
		"saturation frequency can be applicable to a run). Default: "
		<< __defaults[InCol::HIGH_MASS_WIND_SAT_W] << ".";
	__direct_value_options[InCol::HIGH_MASS_WIND_SAT_W]=arg_dbl0(NULL,
			"high-mass-wind-sat-w", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The period at which the wind saturates for low mass "
		"stars in days. This is an alternative to --low-mass-wind-sat-w. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::LOW_MASS_WIND_SAT_P] << "'.";
	__direct_value_options[InCol::LOW_MASS_WIND_SAT_P]=arg_dbl0(NULL,
			"low-mass-wind-sat-p", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The period at which the wind saturates for high mass "
		"stars in days. This is an alternative to --high-mass-wind-sat-w In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::HIGH_MASS_WIND_SAT_P] << "'.";
	__direct_value_options[InCol::HIGH_MASS_WIND_SAT_P]=arg_dbl0(NULL,
			"high-mass-wind-sat-p", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The timescale on which the core end envelope are coupled"
		" in Myr. In --input-columns identified by '"
		<< __input_column_names[InCol::CORE_ENV_COUPLING_TIMESCALE]
		<< "'. Default: " << __defaults[InCol::CORE_ENV_COUPLING_TIMESCALE]
		<< "5.";
	__direct_value_options[InCol::CORE_ENV_COUPLING_TIMESCALE]=arg_dbl0(NULL,
			"core-env-coupling-timescale", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "Log base 10 of the tidal quality factor of the star. In "
		"--input-columns identified by "
		<< __input_column_names[InCol::LGQ] << "'. Default: "
		<< __defaults[InCol::LGQ] << ".";
	__direct_value_options[InCol::LGQ]=arg_dbl0(NULL, "lgQ", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "Mass of the star in solar masses. In --input-columns "
		"identified by '"
		<< __input_column_names[InCol::MSTAR] << "'. Default: "
		<< __defaults[InCol::MSTAR] << "1.";
	__direct_value_options[InCol::MSTAR]=arg_dbl0("M", "Mstar", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "Mass of the planet in jupiter masses. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::MPLANET] << "'. Default: "
		<< __defaults[InCol::MPLANET] << ".";
	__direct_value_options[InCol::MPLANET]=arg_dbl0("m", "Mplanet",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "Radius of the planet in jupiter radii. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::RPLANET] << "'. Default: "
		<< __defaults[InCol::RPLANET] << ".";
	__direct_value_options[InCol::RPLANET]=arg_dbl0("r", "Rplanet",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The age at which the planet forms. If it is smaller "
		"than the disk dissipation age, the planet forms at the disk "
		"dissipation age. In --input-columns identifid by '"
		<< __input_column_names[InCol::PLANET_FORMATION_AGE] << "'. Default "
		<< __defaults[InCol::PLANET_FORMATION_AGE] << ".";
	__direct_value_options[InCol::PLANET_FORMATION_AGE]=arg_dbl0(NULL,
			"planet-formation-age", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The spin frequency at which the star is locked while the"
		" disk is present in rad/day. In --input-columns identifid by '"
		<< __input_column_names[InCol::WDISK] << "'. Default: "
		<< __defaults[InCol::WDISK] << ".";
	__direct_value_options[InCol::WDISK]=arg_dbl0(NULL, "w-disk", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "The spin period at which the star is locked while the "
		"disk is present in days. This is an alternative to --w-disk. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::PDISK] << "'.";
	__direct_value_options[InCol::PDISK]=arg_dbl0(NULL, "p-disk", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "The age at which the disk dissipates in Myr. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::TDISK] << "'. Default: "
		<< __defaults[InCol::TDISK] << ".";
	__direct_value_options[InCol::TDISK]=arg_dbl0(NULL, "t-disk", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "The semimajor axis at which the planet first "
			"appears in AU. In --input-columns identified by '"
			<< __input_column_names[InCol::A_FORMATION] << "'. Default: "
			<< __defaults[InCol::A_FORMATION] << ".";
	__direct_value_options[InCol::A_FORMATION]=arg_dbl0("a",
			"init-semimajor", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The orbital period at which the planet first "
			"appears in days. This is an alternative to --init-semimajor. In"
			"--input-columns identified by '"
			<< __input_column_names[InCol::P_FORMATION] << "'.";
	__direct_value_options[InCol::P_FORMATION]=arg_dbl0("p",
			"init-orbit-period", "<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The starting age for the evolution in Gyr. If this "
		"argument is used, --w-disk, --w-disk-col, --t-disk and --t-disk-col"
		" are ignored and --init-semimajor or --init-semimajor-col is the "
		"semimajor axis that the planet has at the age specified by this"
		" argument. Further, if this argument is  used, --init-wsurf or "
		"--init-wsurf-col must also be specified and if the evolution "
		"of low mass stars is to be computed --init-wrad or "
		"--init-wrad-col is also required. In --input-columns identified"
		" by '" << __input_column_names[InCol::TSTART] << "'.";
	__direct_value_options[InCol::TSTART]=arg_dbl0(NULL, "t0", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "The maximum end age for the evolution in Gyr. The "
		"evolution could stop earlier if the star moves off the main "
		"sequence before this age. In --input-columns identified by '"
		<< __input_column_names[InCol::TEND] << "'. Default: "
		<< __defaults[InCol::TEND] << ".";
	__direct_value_options[InCol::TEND]=arg_dbl0(NULL, "tmax", "<double>",
			option_help.str().c_str());

	option_help.clear();
	option_help << "Initial spin of the stellar core in rad/day for low"
		" mass stars. This argument is ignored, unless --t0 or --t0-col "
		"is also specified. In --input-columns identified by '"
		<< __input_column_names[InCol::START_WRAD] << "'.";
	__direct_value_options[InCol::START_WRAD]=arg_dbl0(NULL, "init-wrad",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "Initial spin of the stellar surface in rad/day. For"
		" low mass stars this is the rotation of the convective zone, "
		"and for high mass stars it is the unique rotation of the star. "
		"This argument is ignored, unless --t0 or --t0-col is also "
		"specified. In --input-columns identified by '"
		<< __input_column_names[InCol::START_WSURF] << "'.";
	__direct_value_options[InCol::START_WSURF]=arg_dbl0(NULL, "init-wsurf",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "An upper limit to impose on the sover timestep. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::MAX_STEP] << "'.";
	__direct_value_options[InCol::MAX_STEP]=arg_dbl0(NULL, "max-step",
			"<double>", option_help.str().c_str());

	option_help.clear();
	option_help << "The number of digits of precision to require of the "
		"solution. Need not be an integer. In --input-columns identified by "
		"'" << __input_column_names[InCol::PRECISION] << "'. Default: "
		<< __defaults[InCol::PRECISION] << ".";
	__direct_value_options[InCol::PRECISION]=arg_dbl0(NULL, "precision",
			"<double>", option_help.str().c_str());

	__start_locked=arg_lit0(NULL, "start-locked", "Whether the planet should"
		   " start locked to the star. If true, causes the value of "
		   "--init-wsurf-col to be ignored.");

	__output_fname=arg_file0("o", "output", "<file|index>", ("If nothing is "
				"read from an input file, then this argument should be the "
				"name of the file to use to output the evolution. If at "
				"least one quantity is read from a list file, that file "
				"should also contain a column specifying the input file name"
				" for each parameter set, and this column should be "
				"specified via the --input-columns option. In the latter "
				"case, this argument is ignored. Default: "+
				__default_outfname).c_str());

	__input_file_columns=arg_str0(NULL, "input-columns",
			"<comma separated list>", "Allows multiple evolutions to be "
			"calculated with a single run. This argument should be a list of"
			" the columns in an input file that includes all quantities that"
			" change, and also a column specifying the file to output the "
			"corresponding evolution to. Any argument not included in this "
			"list of columns is held constant at the value specified at the "
			"corresponding option value or its default. Columns that should "
			"be skipped should have no identifier listed (i.e. there should "
			"be two consecutive commas). Not all columns in the input file "
			"need to be listed. Any content past the last column identified "
			"in this argument on each line is ignored.");

	option_help.clear();
	option_help << "Specifies what to output. This "
			"argument should be a comma separated list containing some of "
			"the following columns:" << std::endl;
	for(int i=0; i<OutCol::NUM_OUTPUT_QUANTITIES; i++)
		option_help << "\t* " << __output_column_names[i] << ": "
			<< __output_column_descr[i] << std::endl;
	outpion_help << "Default: " << __default_output_columns;
	__output_file_columns=arg_str0(NULL, "output-columns",
			"<comma separated list>", option_help.str().c_str());

	__input_fname=arg_file0("i", "input", "<file>", "The file to read the "
			"parameters of the planet-star systems to calculate evolutions "
			"for. If omitted, standard input is used instead. Any lines "
			"beginning with '#' are ignored. The other lines should start "
			"with at least the number of columns specified in the "
			"--input-columns option with white space only between columns.");

	__serialized_stellar_evolution=arg_file0(NULL, "serialized-stellar-evol",
			"<file>", "The file to read previously serialized stellar "
			"evolution from. Default: 'interp_state_data'.");
}

void CommandLineOptions::set_defaults()
{
	for(int i=0; i<InCol::NUM_REAL_INPUT_QUANTITIES; i++)
		__direct_value_options[i]->dval[0]=__defaults[i];
	__output_fname->filename[0]=__default_outfname.c_str();
	__serialized_stellar_evolution->filename[0]=
		__default_serialized_evol.c_str();
}

template<typename COL_ID_TYPE>
void CommandLineOptions::parse_column_list(char *columns_str,
		const std::string[] column_names, int num_column_names,
		std::vector<COL_ID_TYPE> &columns, bool allow_noname)
{
	std::istringstream instream(columns_str);
	columns.clear();
	while(!instream.eof()) {
		std::string colname;
		std::getline(instream, colname, ',');
		for(int i=0; i<num_column_names && column_names[i]!=colname; i++);
		if(column_names[i]==colname) columns.push_back(i);
		else if(allow_noname && colname=="")
			columns.push_back(num_column_names);
		else throw Error::CommandLine("Unrecognized input column '"+colname+
				"'");
	}
}

void CommandLineOptions::postprocess()
{
	if(__input_file_columns->count>0)
		parse_column_list(__input_file_columns->sval[0],
				__input_column_names, NUM_INPUT_QUANTITIES,
				__input_file_format, true);
	parse_column_list(__output_file_columns->sval[0],
			__output_column_names, NUM_OUTPUT_QUANTITIES,
			__output_file_format);
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
		__parsed_ok=false;
		return;
	}
	postprocess();
	free(help_option);
	free(end);
	__parsed_ok=true;
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

void run(const CommandLineOptions &options,
		const StellarEvolution &stellar_evolution,
		std::istream &input)
{
	double star_mass=options.star_mass(),
		   lgQ=options.lgQ(),
		   high_mass_windK=options.high_mass_windK(),
		   low_mass_windK=options.low_mass_windK(),
		   high_mass_wsat=options.high_mass_wind_saturation(),
		   low_mass_wsat=options.low_mass_wind_saturation(),
		   coupling_timescale=options.core_env_coupling_timescale(),
		   wdisk=options.disk_lock_frequency(),
		   tdisk=options.disk_dissipation_age(),
		   Mplanet=options.planet_mass(),
		   Rplanet=options.planet_radius(),
		   planet_formation_age=options.planet_formation_age(),
		   a_formation=options.planet_formation_semimajor(),
		   start_wrad=options.start_wrad(),
		   start_wsurf=options.start_wsurf(),
		   tstart=options.start_age(),
		   tend=options.end_age();
	std::string outfname=options.output_filename();
	bool done=false;
	size_t input_lineno=0;
	while(!done) {
		if(options.input_from_list()) {
			std::string line("#"), word;
			while(line[0]=='#' && !input.eof()) {
				std::getline(input, line);
				++input_lineno;
			}
			if(input.eof()) return;
			else {
				std::istringstream line_stream(line);
				for(int column=0; !line_stream.eof(); column++) {
					if(column==options.star_mass_column())
						line_stream >> star_mass;
					else if(column==options.lgQ_column()) line_stream >> lgQ;
					else if(column==options.high_mass_windK_column())
						line_stream >> high_mass_windK; 
					else if(column==options.low_mass_windK_column())
						line_stream >> low_mass_windK; 
					else if(column==
							options.high_mass_wind_saturation_column())
						line_stream >> high_mass_wsat;
					else if(column==
							options.core_env_coupling_timescale_column())
						line_stream >> coupling_timescale;
					else if(column==options.disk_lock_frequency_column())
						line_stream >> wdisk;
					else if(column==options.disk_dissipation_age_column())
						line_stream >> tdisk;
					else if(column==options.planet_mass_column())
						line_stream >> Mplanet;
					else if(column==options.planet_radius_column())
						line_stream >> Rplanet;
					else if(column==options.planet_formation_age_column())
						line_stream >> planet_formation_age;
					else if(column==
							options.planet_formation_semimajor_column())
						line_stream >> a_formation;
					else if(column==options.start_wrad_column())
						line_stream >> start_wrad;
					else if(column==options.start_wsurf_column())
						line_stream >> start_wsurf;
					else if(column==options.start_age_column())
						line_stream >> tstart;
					else if(column==options.end_age_column())
						line_stream >> tend;
					else if(column==options.output_filename_column())
						line_stream >> outfname;
					else line_stream >> word;
					if(line_stream.fail()) {
						std::ostringstream msg;
						msg << "Error while parsing column " << column+1
							<< " on line " << input_lineno << " of the "
							"input file '" << options.input_filename()
							<< "'";
						throw Error::IO(msg.str());
					}
				}
			}
		} else done=true;
		double Kwind, wsat;
		if(options.star_mass()>stellar_evolution.get_mass_break()) {
			Kwind=high_mass_windK;
			wsat=high_mass_wsat;
		} else {
			Kwind=low_mass_windK;
			wsat=low_mass_wsat;
		}
		OrbitSolver solver(tstart, tend, options.precision());
		calculate_evolution(star_mass, std::pow(10.0, lgQ), Kwind, wsat,
				coupling_timescale, wdisk, tdisk, Mplanet, Rplanet,
				planet_formation_age, a_formation, start_wrad,
				start_wsurf, tstart, tend, options.max_timestep(),
				options.start_locked(), stellar_evolution, solver,
				outfname);
	}
}

///Calculates a realistic evolution chosen to be comlicated.
int main(int argc, char **argv)
{
	try {
		CommandLineOptions options(argc, argv);
		if(!options) return 1;
		YRECEvolution stellar_evolution;
		stellar_evolution.load_state(
				options.serialized_stellar_evolution());
		run(options, stellar_evolution, options.input());
	} catch(Error::General &err) {
		std::cerr << err.what() << ": " << err.get_message() << std::endl;
		return 2;
	}
}
