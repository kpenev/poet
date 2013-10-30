/**\file
 *
 * \brief Produces an executable that calculates orbital evolutions of
 * planet-star systems.
 *
 */

#include "poet.h"

const std::string CommandLineOptions::__input_column_names[]={
	"K",		//LOW_MASS_WINDK
	"K",		//HIGH_MASS_WINDK
	"wsat",		//LOW_MASS_WIND_SAT_W
	"wsat",		//HIGH_MASS_WIND_SAT_W
	"psat",		//LOW_MASS_WIND_SAT_P
	"psat",		//HIGH_MASS_WIND_SAT_P
	"tcoup",	//CORE_ENV_COUPLING_TIMESCALE
	"lgQ",		//LGQ
	"M",		//MSTAR
	"m",		//MPLANET
	"r",		//RPLANET
	"tform",	//PLANET_FORMATION_AGE
	"wdisk",	//WDISK
	"pdisk",	//PDISK
	"tdisk",	//TDISK
	"aform",	//A_FORMATION
	"pform",	//P_FORMATION
	"t0",		//TSTART
	"tmax",		//TEND
	"wrad0",	//START_WRAD
	"wsurf0",	//START_WSURF
	"maxdt",	//MAX_STEP
	"prec",		//PRECISION
	"outf",		//OUT_FNAME
	"locked",	//START_LOCKED
};

const std::string CommandLineOptions::__output_column_descr[]={
	//AGE
	"Age in Gyr",

	//SEMIMAJOR
	"Semimajor axis in AU",

	//WORB
	"Orbital frequency in rad/day",

	//PORB
	"Orbital period in days",

	//LCONV
	"Convective zone angular momentum (low mass stars only else NaN).",

	//LRAD
	"Radiative zone angular momentum (low mass stars only else NaN).",

	//LTOT
	"Total angular momentum of the star.",

	//ICONV
	"Convective zone moment of inertia (low mass stars only else NaN).",

	//IRAD
	"Radiative zone moment of inertia (low mass stars only else NaN).",

	//ITOT
	"Total moment of inertia of the star.",

	//WSURF
	"Stellar surface angular velocity of the star in rad/day.",

	//WRAD
	"Angular velocity of the radiative zone in rad/day (low mass stars only,"
		", else NaN)",

	//PSURF
	"Stellar surface spin period in days.",

	//PRAD
	"Spin period of the radiative zone in days (low mass stars only, else "
		"NaN)",

	//EVOL_MODE
	"Evolution mode of the step starting at this age.",

	//WIND_STATE
	"The wind state (saturated/not saturated) of the step starting at this "
		"age.",

	//RSTAR
	"Stellar radius",

	//LSTAR
	"Stellar luminosity",

	//RRAD
	"Radius of the stellar radiative zone (low mass stars only, else NaN).",

	//MRAD
	"Mass of the stellar radiative zone (low mass stars only, else NaN)."};

const double CommandLineOptions::__defaults[]={
		0.35,			//LOW_MASS_WINDK=WINDK
		0,				//HIGH_MASS_WINDK
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
		Inf,			//MAX_STEP
		6				//PRECISION
};

const std::string CommandLineOptions::__default_outfname="poet.evol",
	  CommandLineOptions::__default_serialized_evol="interp_state_data",
	  CommandLineOptions::__default_output_columns=
	  	"t,a,Lconv,Lrad,L,Iconv,Irad,I,mode";

char *CommandLineOptions::cstr_copy(const std::ostringstream &stream)
{
	char *cstr=new char[stream.str().length()+1];
	std::strcpy(cstr, stream.str().c_str());
	__option_help_copies.push_back(cstr);
	return cstr;
}

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
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The wind strength for high mass stars in units of "
			"Msun*Rsun^2*day^2/(rad^2*Gyr). To identify this quantity in "
			"--input-columns use '"
			<<  __input_column_names[InCol::HIGH_MASS_WINDK]
			<< "' (identical to the --low-mass-K since only one K can be "
			"applicable to a run). Default: "
			<< __defaults[InCol::LOW_MASS_WINDK] << ".";
	__direct_value_options[InCol::HIGH_MASS_WINDK]=arg_dbl0(NULL,
			"high-mass-K", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The frequency at which the wind saturates for low mass "
		"stars in units of rad/day. To identify this quantity in "
		"--input-columns use '"
		<<  __input_column_names[InCol::LOW_MASS_WIND_SAT_W]
		<< "' (identical to the --high-mass-wind-sat-w since only one "
		"saturation frequency can be applicable to a run). Default: "
		<< __defaults[InCol::LOW_MASS_WIND_SAT_W] << ".";
	__direct_value_options[InCol::LOW_MASS_WIND_SAT_W]=arg_dbl0(NULL,
			"low-mass-wind-sat-w", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The frequency at which the wind saturates for high mass "
		"stars in units of rad/day. To identify this quantity in "
		"--input-columns use '"
		<< __input_column_names[InCol::HIGH_MASS_WIND_SAT_W]
		<< "' (identical to the --high-mass-wind-sat-w since only one "
		"saturation frequency can be applicable to a run). Default: "
		<< __defaults[InCol::HIGH_MASS_WIND_SAT_W] << ".";
	__direct_value_options[InCol::HIGH_MASS_WIND_SAT_W]=arg_dbl0(NULL,
			"high-mass-wind-sat-w", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The period at which the wind saturates for low mass "
		"stars in days. This is an alternative to --low-mass-wind-sat-w. "
		"If both options are specified, --low-mass-wind-sat-w is used. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::LOW_MASS_WIND_SAT_P] << "'.";
	__direct_value_options[InCol::LOW_MASS_WIND_SAT_P]=arg_dbl0(NULL,
			"low-mass-wind-sat-p", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The period at which the wind saturates for high mass "
		"stars in days. This is an alternative to --high-mass-wind-sat-w. "
		"If both options are specified, --high-mass-wind-sat-w is used. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::HIGH_MASS_WIND_SAT_P] << "'.";
	__direct_value_options[InCol::HIGH_MASS_WIND_SAT_P]=arg_dbl0(NULL,
			"high-mass-wind-sat-p", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The timescale on which the core end envelope are coupled"
		" in Myr. In --input-columns identified by '"
		<< __input_column_names[InCol::CORE_ENV_COUPLING_TIMESCALE]
		<< "'. Default: " << __defaults[InCol::CORE_ENV_COUPLING_TIMESCALE]
		<< "5.";
	__direct_value_options[InCol::CORE_ENV_COUPLING_TIMESCALE]=arg_dbl0(NULL,
			"core-env-coupling-timescale", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "Log base 10 of the tidal quality factor of the star. In "
		"--input-columns identified by "
		<< __input_column_names[InCol::LGQ] << "'. Default: "
		<< __defaults[InCol::LGQ] << ".";
	__direct_value_options[InCol::LGQ]=arg_dbl0(NULL, "lgQ", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "Mass of the star in solar masses. In --input-columns "
		"identified by '"
		<< __input_column_names[InCol::MSTAR] << "'. Default: "
		<< __defaults[InCol::MSTAR] << "1.";
	__direct_value_options[InCol::MSTAR]=arg_dbl0("M", "Mstar", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "Mass of the planet in jupiter masses. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::MPLANET] << "'. Default: "
		<< __defaults[InCol::MPLANET] << ".";
	__direct_value_options[InCol::MPLANET]=arg_dbl0("m", "Mplanet",
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "Radius of the planet in jupiter radii. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::RPLANET] << "'. Default: "
		<< __defaults[InCol::RPLANET] << ".";
	__direct_value_options[InCol::RPLANET]=arg_dbl0("r", "Rplanet",
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The age at which the planet forms. If it is smaller "
		"than the disk dissipation age, the planet forms at the disk "
		"dissipation age. In --input-columns identifid by '"
		<< __input_column_names[InCol::PLANET_FORMATION_AGE] << "'. Default "
		<< __defaults[InCol::PLANET_FORMATION_AGE] << ".";
	__direct_value_options[InCol::PLANET_FORMATION_AGE]=arg_dbl0(NULL,
			"planet-formation-age", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The spin frequency at which the star is locked while the"
		" disk is present in rad/day. In --input-columns identifid by '"
		<< __input_column_names[InCol::WDISK] << "'. Default: "
		<< __defaults[InCol::WDISK] << ".";
	__direct_value_options[InCol::WDISK]=arg_dbl0(NULL, "w-disk", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "The spin period at which the star is locked while the "
		"disk is present in days. This is an alternative to --w-disk. "
		"If both options are specified, --w-disk is used. In --input-columns"
		" identified by '" << __input_column_names[InCol::PDISK] << "'.";
	__direct_value_options[InCol::PDISK]=arg_dbl0(NULL, "p-disk", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "The age at which the disk dissipates in Myr. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::TDISK] << "'. Default: "
		<< __defaults[InCol::TDISK] << ".";
	__direct_value_options[InCol::TDISK]=arg_dbl0(NULL, "t-disk", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "The semimajor axis at which the planet first "
			"appears in AU. In --input-columns identified by '"
			<< __input_column_names[InCol::A_FORMATION] << "'. Default: "
			<< __defaults[InCol::A_FORMATION] << ".";
	__direct_value_options[InCol::A_FORMATION]=arg_dbl0("a",
			"init-semimajor", "<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The orbital period at which the planet first "
			"appears in days. This is an alternative to --init-semimajor. In"
			"--input-columns identified by '"
			<< __input_column_names[InCol::P_FORMATION] << "'.";
	__direct_value_options[InCol::P_FORMATION]=arg_dbl0("p",
			"init-orbit-period", "<double>", cstr_copy(option_help));

	option_help.str("");
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
			cstr_copy(option_help));

	option_help.str("");
	option_help << "The maximum end age for the evolution in Gyr. The "
		"evolution could stop earlier if the star moves off the main "
		"sequence before this age. In --input-columns identified by '"
		<< __input_column_names[InCol::TEND] << "'. Default: "
		<< __defaults[InCol::TEND] << ".";
	__direct_value_options[InCol::TEND]=arg_dbl0(NULL, "tmax", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "Initial spin of the stellar core in rad/day for low"
		" mass stars. This argument is ignored, unless --t0 or --t0-col "
		"is also specified. In --input-columns identified by '"
		<< __input_column_names[InCol::START_WRAD] << "'.";
	__direct_value_options[InCol::START_WRAD]=arg_dbl0(NULL, "init-wrad",
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "Initial spin of the stellar surface in rad/day. For"
		" low mass stars this is the rotation of the convective zone, "
		"and for high mass stars it is the unique rotation of the star. "
		"This argument is ignored, unless --t0 or --t0-col is also "
		"specified. In --input-columns identified by '"
		<< __input_column_names[InCol::START_WSURF] << "'.";
	__direct_value_options[InCol::START_WSURF]=arg_dbl0(NULL, "init-wsurf",
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "An upper limit to impose on the sover timestep. In "
		"--input-columns identified by '"
		<< __input_column_names[InCol::MAX_STEP] << "'. Default: "
		<< __defaults[InCol::MAX_STEP] << ".";
	__direct_value_options[InCol::MAX_STEP]=arg_dbl0(NULL, "max-step",
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "The number of digits of precision to require of the "
		"solution. Need not be an integer. In --input-columns identified by "
		"'" << __input_column_names[InCol::PRECISION] << "'. Default: "
		<< __defaults[InCol::PRECISION] << ".";
	__direct_value_options[InCol::PRECISION]=arg_dbl0(NULL, "precision",
			"<double>", cstr_copy(option_help));

	__start_locked=arg_lit0(NULL, "start-locked", "Whether the planet should"
		   " start locked to the star. If true, causes the value of "
		   "--init-wsurf-col to be ignored.");

	option_help.str("");
	option_help << "If nothing is read from an input file, then this "
		"argument should be the name of the file to use to output the "
		"evolution. If at least one quantity is read from a list file, that "
		"file should also contain a column specifying the input file name"
		" for each parameter set, and this column should be specified via "
		"the --input-columns option. In the latter case, this argument is "
		"ignored. Default: " << __default_outfname;
	__output_fname=arg_file0("o", "output", "<file|index>",
			cstr_copy(option_help));

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

	option_help.str("");
	option_help << "Specifies what to output. This "
			"argument should be a comma separated list containing some of "
			"the following columns:" << std::endl;
	for(int i=0; i<OutCol::NUM_OUTPUT_QUANTITIES; i++)
		option_help << "\t* " << OUTPUT_COLUMN_NAMES[i] << ": "
			<< __output_column_descr[i] << std::endl;
	option_help << "Default: " << __default_output_columns;
	__output_file_columns=arg_str0(NULL, "output-columns",
			"<comma separated list>", cstr_copy(option_help));

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
	__output_file_columns->sval[0]=__default_output_columns.c_str();
}

template<typename COL_ID_TYPE>
void CommandLineOptions::parse_column_list(const char *columns_str,
		const std::string *column_names, int num_column_names,
		std::vector<COL_ID_TYPE> &columns, bool allow_noname)
{
	std::istringstream instream(columns_str);
	columns.clear();
	while(!instream.eof()) {
		std::string colname;
		std::getline(instream, colname, ',');
		int i=0;
		while(i<num_column_names && column_names[i]!=colname) i++;
		if(column_names[i]==colname)
			columns.push_back(static_cast<COL_ID_TYPE>(i));
		else if(allow_noname && colname=="")
			columns.push_back(static_cast<COL_ID_TYPE>(num_column_names));
		else throw Error::CommandLine("Unrecognized input column '"+colname+
				"'");
	}
}

void CommandLineOptions::postprocess()
{
	if(__input_file_columns->count>0)
		parse_column_list(__input_file_columns->sval[0],
				__input_column_names, InCol::NUM_INPUT_QUANTITIES,
				__input_file_format, true);
	parse_column_list(__output_file_columns->sval[0],
			OUTPUT_COLUMN_NAMES, OutCol::NUM_OUTPUT_QUANTITIES,
			__output_file_format);
	if(__input_fname->count) {
		__input_stream.open(__input_fname->filename[0]);
		__opened_stream=true;
	}
		
	if(__direct_value_options[InCol::LOW_MASS_WIND_SAT_W]->count==0 &&
			__direct_value_options[InCol::LOW_MASS_WIND_SAT_P]->count>0)
		__direct_value_options[InCol::LOW_MASS_WIND_SAT_W]->dval[0]=2.0*M_PI/
			__direct_value_options[InCol::LOW_MASS_WIND_SAT_P]->dval[0];
	else if(__direct_value_options[InCol::LOW_MASS_WIND_SAT_W]->count>0)
		__direct_value_options[InCol::LOW_MASS_WIND_SAT_P]->dval[0]=2.0*M_PI/
			__direct_value_options[InCol::LOW_MASS_WIND_SAT_W]->dval[0];

	if(__direct_value_options[InCol::WDISK]->count==0 &&
			__direct_value_options[InCol::PDISK]->count>0)
		__direct_value_options[InCol::WDISK]->dval[0]=2.0*M_PI/
			__direct_value_options[InCol::PDISK]->dval[0];
	else if(__direct_value_options[InCol::WDISK]->count>0)
		__direct_value_options[InCol::PDISK]->dval[0]=2.0*M_PI/
			__direct_value_options[InCol::WDISK]->dval[0];

	if(__direct_value_options[InCol::A_FORMATION]->count>0 ||
			__direct_value_options[InCol::P_FORMATION]->count>0) {
		double M=__direct_value_options[InCol::MSTAR]->dval[0]*
					AstroConst::solar_mass,
			m=__direct_value_options[InCol::MPLANET]->dval[0]*
					AstroConst::jupiter_mass;
		if(__direct_value_options[InCol::A_FORMATION]->count==0 &&
				__direct_value_options[InCol::P_FORMATION]->count>0) {
				double P2=std::pow(
						__direct_value_options[InCol::P_FORMATION]->dval[0]*
						AstroConst::day, 2);
				__direct_value_options[InCol::A_FORMATION]->dval[0]=
					std::pow(AstroConst::G*(M+m)*P2/4/M_PI/M_PI,
							1.0/3.0)/AstroConst::AU;
		} else if(__direct_value_options[InCol::A_FORMATION]->count>0) {
			double a3=std::pow(
					__direct_value_options[InCol::A_FORMATION]->dval[0]*
					AstroConst::AU, 3);
			__direct_value_options[InCol::P_FORMATION]->dval[0]=std::sqrt(
					4.0*M_PI*M_PI*a3/(M+m)/AstroConst::G)/AstroConst::day;
		}
	}
}

void CommandLineOptions::cleanup()
{
	if(__opened_stream) __input_stream.close();
	arg_freetable(__argtable, sizeof(__argtable)/sizeof(__argtable[0]));
	for(std::list<char*>::iterator i=__option_help_copies.begin();
			i!=__option_help_copies.end(); i++) delete[] *i;
}

CommandLineOptions::CommandLineOptions(int argc, char **argv) :
	__direct_value_options(InCol::NUM_REAL_INPUT_QUANTITIES),
	__opened_stream(false)
{
	define_options();
	arg_lit *help_option=arg_lit0("h", "help", "Print this help and exit.");
	struct arg_end *end = arg_end(100);
	for(int i=0; i<InCol::NUM_REAL_INPUT_QUANTITIES; i++)
		__argtable[i]=__direct_value_options[i];
	__argtable[InCol::OUT_FNAME]=__output_fname;
	__argtable[InCol::START_LOCKED]=__start_locked;
	__argtable[InCol::NUM_INPUT_QUANTITIES]=__input_file_columns;
	__argtable[InCol::NUM_INPUT_QUANTITIES+1]=__output_file_columns;
	__argtable[InCol::NUM_INPUT_QUANTITIES+2]=__input_fname;
	__argtable[InCol::NUM_INPUT_QUANTITIES+3]=__serialized_stellar_evolution;
	__argtable[InCol::NUM_INPUT_QUANTITIES+4]=help_option;
	__argtable[InCol::NUM_INPUT_QUANTITIES+5]=end;

	if(arg_nullcheck(__argtable) != 0) {
		cleanup();
		throw Error::CommandLine("Failed to allocate argument table.");
	}
	set_defaults();
	int nerrors=arg_parse(argc, argv, __argtable);
	if(help_option->count>0 || nerrors>0) {
        printf("Usage: %s", "SubPixPhot");
        arg_print_syntax(stdout, __argtable,"\n");
        arg_print_glossary(stdout, __argtable,"  %-25s %s\n");
		if(help_option->count==0)
			arg_print_errors(stdout, end, "poet");
		__parsed_ok=false;
		return;
	}
	postprocess();
	__parsed_ok=true;
}

double CommandLineOptions::get_real_value(InCol::InputColumns quantity) const
{
	if(quantity<0 || quantity>=InCol::NUM_REAL_INPUT_QUANTITIES)
		throw Error::BadFunctionArguments("Unrecognized real quantity in "
				"CommandLineOptions::get_real_value().");
	return __direct_value_options[quantity]->dval[0];
}





///AU/\f$\mathrm{R}_\odot\f$.
const double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;

void output_solution(const OrbitSolver &solver, const StellarSystem &system,
		const std::string &filename,
		const std::vector<OutCol::OutputColumns> &output_file_format)
{

	std::list<double>::const_iterator
		age_i=solver.get_tabulated_var(AGE)->begin(),
		a_i=solver.get_tabulated_var(SEMIMAJOR)->begin(),
		Lconv_i=solver.get_tabulated_var(LCONV)->begin(),
		Lrad_i=solver.get_tabulated_var(LRAD)->begin();
	const Star &star=system.get_star();
	const Planet &planet=system.get_planet();
	std::list<double>::const_iterator
		last_age=solver.get_tabulated_var(AGE)->end();

	std::list<EvolModeType>::const_iterator
		mode_i=solver.get_tabulated_evolution_mode()->begin();
	std::list<WindSaturationState>::const_iterator
		wind_i=solver.get_tabulated_wind_state()->begin();

	std::ofstream outf(filename.c_str());
	outf.precision(16);
	outf << "#";
	for(size_t i=0; i<output_file_format.size(); i++)
		outf << std::setw(i==0 ? 24 : 25)
            << OUTPUT_COLUMN_NAMES[output_file_format[i]];
	outf << std::endl;

	while(age_i!=last_age) {
		for(size_t i=0; i<output_file_format.size(); i++) {
			double Iconv=star.moment_of_inertia(*age_i, convective),
				   Irad=star.moment_of_inertia(*age_i, radiative);
			outf << std::setw(25);
			double semimajor=(*a_i)/AU_Rsun,
				   worb=planet.orbital_angular_velocity_semimajor(semimajor);
			switch(output_file_format[i]) {
				case OutCol::AGE : outf << *age_i; break;
				case OutCol::SEMIMAJOR : outf << semimajor; break;
				case OutCol::WORB : outf << worb; break;
				case OutCol::PORB : outf << 2.0*M_PI/worb; break;
				case OutCol::LCONV : outf << *Lconv_i; break;
				case OutCol::LRAD : outf << *Lrad_i; break;
				case OutCol::LTOT : outf << *Lconv_i + *Lrad_i; break;
				case OutCol::ICONV : outf << Iconv; break;
				case OutCol::IRAD : outf << Irad; break;
				case OutCol::ITOT : outf << Iconv+Irad; break;
				case OutCol::WSURF : outf << (*Lconv_i)/Iconv; break;
				case OutCol::WRAD : outf << (*Lrad_i)/Irad; break;
				case OutCol::PSURF : outf << 2.0*M_PI*Iconv/(*Lconv_i);
									 break;
				case OutCol::PRAD : outf << 2.0*M_PI*Irad/(*Lrad_i); break;
				case OutCol::EVOL_MODE : outf << *mode_i; break;
				case OutCol::WIND_STATE : outf << *wind_i; break;
				case OutCol::RSTAR : outf << star.get_radius(*age_i); break;
				case OutCol::LSTAR : outf << star.get_luminosity(*age_i);
									 break;
				case OutCol::RRAD : outf << star.get_rad_radius(*age_i);
									break;
				case OutCol::MRAD : outf << star.get_rad_mass(*age_i);
									break;
				default : throw Error::BadFunctionArguments(
								  "Unrecognized output column encountered in"
								  " output_file_format in "
								  "poet.cpp:output_solution.");
			}
		}
		outf << std::endl;
		++age_i; ++a_i; ++Lconv_i; ++Lrad_i; ++mode_i; ++wind_i;
	}
	outf.close();
}

void calculate_evolution(const std::vector<double> &real_parameters,
		bool start_locked, const StellarEvolution &stellar_evolution,
		const std::string &outfname,
		const std::vector<OutCol::OutputColumns> &output_file_format)
{
	Star star(real_parameters[InCol::MSTAR],
			std::pow(10.0, real_parameters[InCol::LGQ]),
			real_parameters[InCol::WINDK],
			real_parameters[InCol::WIND_SAT_W],
			real_parameters[InCol::CORE_ENV_COUPLING_TIMESCALE]*1e-3,
			1e-3, //HACK
			real_parameters[InCol::WDISK],
			real_parameters[InCol::TDISK]*1e-3,
			stellar_evolution);
	double tstart=real_parameters[InCol::TSTART];
	if(star.is_low_mass() && 
			star.core_formation_age()>star.get_disk_dissipation_age() &&
			tstart<star.core_formation_age())
		throw Error::Runtime("At present the case when the disk dissipates "
				"before the stellar core starts to form is not supported.");
	Planet planet(&star, real_parameters[InCol::MPLANET],
			real_parameters[InCol::RPLANET],
			real_parameters[InCol::A_FORMATION]);
	StellarSystem system(&star, &planet);
	std::valarray<double> start_orbit(0.0, 1);
	EvolModeType start_evol_mode;
	if(std::isnan(tstart) || tstart<star.get_disk_dissipation_age())
		start_evol_mode=LOCKED_TO_DISK;
	if(std::isnan(tstart) ||
			(star.is_low_mass() && tstart<star.core_formation_age()))
		tstart=NaN;
	else if(tstart<star.get_disk_dissipation_age()) {
		start_orbit[0]=(star.is_low_mass() ?
				real_parameters[InCol::START_WRAD] :
				real_parameters[InCol::START_WSURF])*
			star.moment_of_inertia(tstart, radiative);
	} else if(start_locked) {
		start_evol_mode=LOCKED_TO_PLANET;
		start_orbit.resize(2);
		start_orbit[0]=real_parameters[InCol::A_FORMATION]*AU_Rsun;
		start_orbit[1]=(star.is_low_mass() ?
				real_parameters[InCol::START_WRAD]*
				star.moment_of_inertia(tstart, radiative) : 0);
	} else {
		if(planet.orbital_angular_velocity_semimajor(
					real_parameters[InCol::A_FORMATION])<
				real_parameters[InCol::START_WSURF])
			start_evol_mode=SLOW_PLANET;
		else start_evol_mode=FAST_PLANET;
		start_orbit.resize(3);
		start_orbit[0]=
			std::pow(real_parameters[InCol::A_FORMATION]*AU_Rsun, 6.5);
		start_orbit[1]=real_parameters[InCol::START_WSURF]*
			star.moment_of_inertia(tstart, convective);
		start_orbit[2]=(star.is_low_mass() ?
				real_parameters[InCol::START_WRAD]*
				star.moment_of_inertia(tstart, radiative) : 0);
	}
    double tend=std::min(real_parameters[InCol::TEND], star.get_lifetime());
	OrbitSolver solver(tstart, tend,
			std::pow(10.0, -real_parameters[InCol::PRECISION]));
	solver(system, real_parameters[InCol::MAX_STEP],
			real_parameters[InCol::PLANET_FORMATION_AGE],
			real_parameters[InCol::A_FORMATION], tstart,
			start_evol_mode, start_orbit);
	output_solution(solver, system, outfname, output_file_format);
}

std::string update_run_parameters(std::vector<double> &real_parameters,
		bool &start_locked,
		const std::vector<InCol::InputColumns> &input_format,
		std::istringstream &line, size_t input_lineno)
{
	std::string outfname;
	bool read_wind_sat_w=false, read_wind_sat_p=false, read_wdisk=false,
		 read_pdisk=false, read_a_formation=false, read_p_formation=false;
	std::string word;
	for(size_t column=0; column<input_format.size(); column++) {
		InCol::InputColumns quantity=input_format[column];
		if(quantity<InCol::NUM_REAL_INPUT_QUANTITIES) {
			line >> real_parameters[quantity];
			if(quantity==InCol::WINDK)
				real_parameters[InCol::HIGH_MASS_WINDK]=
					real_parameters[InCol::WINDK];
			else if(quantity==InCol::WIND_SAT_W) {
				real_parameters[InCol::HIGH_MASS_WIND_SAT_W]=
					real_parameters[InCol::WIND_SAT_W];
				read_wind_sat_w=true;
			} else if(quantity==InCol::WIND_SAT_P) {
				real_parameters[InCol::HIGH_MASS_WIND_SAT_P]=
					real_parameters[InCol::WIND_SAT_P];
				read_wind_sat_p=true;
			} else if(quantity==InCol::WDISK) read_wdisk=true;
			else if(quantity==InCol::PDISK) read_pdisk=true;
			else if(quantity==InCol::A_FORMATION) read_a_formation=true;
			else if(quantity==InCol::P_FORMATION) read_p_formation=true;
		} else if(quantity==InCol::OUT_FNAME) line >> outfname;
		else if(quantity==InCol::START_LOCKED) line >> start_locked;
		else if(quantity==InCol::SKIP) line >> word;
		else throw Error::BadFunctionArguments(
				"Unrecognized input quantity tag in "
				"poet.cpp:run().");
		if(line.fail()) {
			std::ostringstream msg;
			msg << "Error while parsing column " << column+1
				<< " on line " << input_lineno << " of the "
				"input file.";
			throw Error::IO(msg.str());
		}
	}
	if(read_wind_sat_w) real_parameters[InCol::WIND_SAT_P]=
		2.0*M_PI/real_parameters[InCol::WIND_SAT_W];
	else if(read_wind_sat_p) real_parameters[InCol::WIND_SAT_W]=
		2.0*M_PI/real_parameters[InCol::WIND_SAT_P];
	if(read_wdisk) real_parameters[InCol::PDISK]=
		2.0*M_PI/real_parameters[InCol::WDISK];
	else if(read_pdisk) real_parameters[InCol::WDISK]=
		2.0*M_PI/real_parameters[InCol::PDISK];
	double M=real_parameters[InCol::MSTAR]*AstroConst::solar_mass,
		m=real_parameters[InCol::MPLANET]*AstroConst::jupiter_mass;
	if(read_a_formation) {
		double a3=std::pow(
                real_parameters[InCol::A_FORMATION]*AstroConst::AU, 3);
		real_parameters[InCol::P_FORMATION]=std::sqrt(
					4.0*M_PI*M_PI*a3/(M+m)/AstroConst::G)/AstroConst::day;
	} else if(read_p_formation) {
		double P2=std::pow(
                real_parameters[InCol::P_FORMATION]*AstroConst::day, 2);
		real_parameters[InCol::A_FORMATION]=
					std::pow(AstroConst::G*(M+m)*P2/4/M_PI/M_PI, 1.0/3.0)/
					AstroConst::AU;
	}
	return outfname;
}

void run(const CommandLineOptions &options,
		const StellarEvolution &stellar_evolution,
		std::istream &input)
{
	std::vector<double> real_parameters(InCol::NUM_REAL_INPUT_QUANTITIES);
	for(int i=0; i<InCol::NUM_REAL_INPUT_QUANTITIES; i++) 
		real_parameters[i]=options.get_real_value(
				static_cast<InCol::InputColumns>(i));

	std::string outfname=options.output_filename();
	bool done=false, start_locked=options.start_locked();
	size_t input_lineno=0;
	while(!done) {
		if(options.input_from_list()) {
			std::string line("#");
			while(line[0]=='#' && !input.eof()) {
				std::getline(input, line);
				++input_lineno;
			}
			if(input.eof()) return;
			else {
				std::istringstream line_str(line);
				outfname=update_run_parameters(real_parameters,
					start_locked, options.input_file_format(),
					line_str, input_lineno);
			}
		} else done=true;
		if(real_parameters[InCol::MSTAR]>stellar_evolution.get_mass_break()){
			real_parameters[InCol::WINDK]=
				real_parameters[InCol::HIGH_MASS_WINDK];
			real_parameters[InCol::WIND_SAT_W]=
				real_parameters[InCol::HIGH_MASS_WIND_SAT_W];
		} else {
			real_parameters[InCol::WINDK]=
				real_parameters[InCol::LOW_MASS_WINDK];
			real_parameters[InCol::WIND_SAT_W]=
				real_parameters[InCol::LOW_MASS_WIND_SAT_W];
		}
		calculate_evolution(real_parameters, start_locked, stellar_evolution,
				outfname, options.output_file_format());
	}
}

///Calculates a realistic evolution chosen to be comlicated.
int main(int argc, char **argv)
{
	try {
#ifdef DEBUG_STOPPING
		std::cerr.precision(16);
		std::cerr.setf(std::ios_base::scientific);
#endif
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
