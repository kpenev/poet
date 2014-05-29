/**\file
 *
 * \brief Produces an executable that calculates orbital evolutions of
 * planet-star systems.
 *
 */

#include "poet.h"

void CommandLineOptions::init_input_column_names()
{
	__input_column_names.resize(InCol::NUM_INPUT_QUANTITIES);
	__input_column_names[InCol::LOW_MASS_WINDK]="K";
	__input_column_names[InCol::HIGH_MASS_WINDK]="K";
	__input_column_names[InCol::LOW_MASS_WIND_SAT_W]="wsat";
	__input_column_names[InCol::HIGH_MASS_WIND_SAT_W]="wsat";
	__input_column_names[InCol::LOW_MASS_WIND_SAT_P]="psat";
	__input_column_names[InCol::HIGH_MASS_WIND_SAT_P]="psat";
	__input_column_names[InCol::CORE_ENV_COUPLING_TIMESCALE]="tcoup";
	__input_column_names[InCol::LGQ]="lgQ";
	__input_column_names[InCol::LGQ_INERTIAL]="lgQinr";
	__input_column_names[InCol::MSTAR]="M";
	__input_column_names[InCol::MPLANET]="m";
	__input_column_names[InCol::RPLANET]="r";
	__input_column_names[InCol::PLANET_FORMATION_AGE]="tform";
	__input_column_names[InCol::WDISK]="wdisk";
	__input_column_names[InCol::PDISK]="pdisk";
	__input_column_names[InCol::TDISK]="tdisk";
	__input_column_names[InCol::A_FORMATION]="aform";
	__input_column_names[InCol::P_FORMATION]="pform";
	__input_column_names[InCol::INCLINATION_FORMATION]="thetaform";
	__input_column_names[InCol::TSTART]="t0";
	__input_column_names[InCol::TEND]="tmax";
	__input_column_names[InCol::START_WRAD]="wrad0";
	__input_column_names[InCol::START_WSURF]="wsurf0";
	__input_column_names[InCol::MAX_STEP]="maxdt";
	__input_column_names[InCol::PRECISION]="prec";
	__input_column_names[InCol::OUT_FNAME]="outf";
	__input_column_names[InCol::START_LOCKED]="locked";
	__input_column_names[InCol::REQUIRED_AGES]="trequired";
#ifdef COLUMN_NAME_EMPHASIS
	for(int i=0; i<InCol::NUM_INPUT_QUANTITIES; ++i)
		__input_column_names[i]=std::string(COLUMN_NAME_EMPHASIS) +
								__input_column_names[i] +
								std::string(COLUMN_NAME_EMPHASIS);
#endif
}

void CommandLineOptions::init_output_column_descriptions()
{
	__output_column_descr.resize(OutCol::NUM_OUTPUT_QUANTITIES);
	__output_column_descr[OutCol::AGE]="Age in Gyr";
	__output_column_descr[OutCol::SEMIMAJOR]="Semimajor axis in AU";
	__output_column_descr[OutCol::INCLINATION]="The stellar spin-orbit "
		"inclination";
	__output_column_descr[OutCol::WORB]="Orbital frequency in rad/day";
	__output_column_descr[OutCol::PORB]="Orbital period in days";
	__output_column_descr[OutCol::LORB]="Orbital angular momentum in "
		"Msun Rsun^2 rad/day";
	__output_column_descr[OutCol::LCONV]=
		"Convective zone angular momentum in Msun Rsun^2 rad/day "
		"(low mass stars only else NaN).";
	__output_column_descr[OutCol::LRAD]=
		"Radiative zone angular momentum in Msun Rsun^2 rad/day "
		"(low mass stars only else NaN).";
	__output_column_descr[OutCol::LRAD_PAR]=
		"Radiative zone angular momentum component along the rotation axis "
		"of the convective zone in Msun Rsun^2 rad/day (low mass stars only "
		"else NaN).";
	__output_column_descr[OutCol::LRAD_PERP]=
		"Radiative zone angular momentum component perpendicular to the "
		"rotation axis of the convective zone in Msun Rsun^2 rad/day (low "
		"mass stars only else NaN).";
	__output_column_descr[OutCol::LTOT]=
		"Total angular momentum of the star in Msun Rsun^2 rad/day.";
	__output_column_descr[OutCol::ICONV]=
		"Convective zone moment of inertia Msun Rsun^2 "
		"(low mass stars only else NaN).";
	__output_column_descr[OutCol::IRAD]=
		"Radiative zone moment of inertia Msun Rsun^2 "
		"(low mass stars only else NaN).";
	__output_column_descr[OutCol::ITOT]=
		"Total moment of inertia of the star Msun Rsun^2.";
	__output_column_descr[OutCol::WSURF]=
		"Stellar surface angular velocity of the star in rad/day.";
	__output_column_descr[OutCol::WRAD]="Angular velocity of the radiative "
		"zone in rad/day (low mass stars only, else NaN)";
	__output_column_descr[OutCol::PSURF]=
		"Stellar surface spin period in days.";
	__output_column_descr[OutCol::PRAD]="Spin period of the radiative zone "
		"in days (low mass stars only, else NaN)";
	__output_column_descr[OutCol::EVOL_MODE]=
		"Evolution mode of the step starting at this age.";
	__output_column_descr[OutCol::WIND_STATE]="The wind state (saturated/not"
	   " saturated) of the step starting at this age.";
	__output_column_descr[OutCol::RSTAR]="Stellar radius in solar radii";
	__output_column_descr[OutCol::LSTAR]=
		"Stellar luminosity in solar luminosities";
	__output_column_descr[OutCol::RRAD]="Radius of the stellar radiative "
		"zone in solar radii (low mass stars only, else NaN).";
	__output_column_descr[OutCol::MRAD]="Mass of the stellar radiative zone "
		"in solar masses (low mass stars only, else NaN).";
	__output_column_descr[OutCol::ICONV_DERIV]="First order derivative with "
		"respect to age of the moment of inertia of the stallar surface "
		"convective zone in Msun Rsun^2/Gyr "
		"(low mass stars only, else NaN).";
	__output_column_descr[OutCol::IRAD_DERIV]="First order derivative with "
		"respect to age of the moment of inertia of the stallar radiative "
		"core in Msun Rsun^2/Gyr (low mass stars only, else NaN).";
	__output_column_descr[OutCol::ITOT_DERIV]="First order derivative with "
		"respect to age of the moment of inertia of the entire star in Msun "
		"Rsun^2/Gyr.";
	__output_column_descr[OutCol::RSTAR_DERIV]="First order derivative with "
		"respect to age of the stellar radius in solar radii per Gyr.";
	__output_column_descr[OutCol::RRAD_DERIV]="First order derivative with "
		"respect to age of the radius of the stellar radiative core in solar"
		" radii per Gyr (low mass stars only, else NaN).";
	__output_column_descr[OutCol::MRAD_DERIV]="First order derivative with "
		"respect to age of the mass of the stellar radiative core in solar "
		"masses per Gyr (low mass stars only, else NaN).";
	__output_column_descr[OutCol::ICONV_SECOND_DERIV]=
		"Second order derivative with respect to age of the moment of "
		"inertia of the stallar surface convective zone in Msun Rsun^2/Gyr^2"
		" (low mass stars only, else NaN).";
	__output_column_descr[OutCol::IRAD_SECOND_DERIV]=
		"Second order derivative with respect to age of the moment of "
		"inertia of the stallar radiative core in Msun Rsun^2/Gyr^2"
		" (low mass stars only, else NaN).";
	__output_column_descr[OutCol::ITOT_SECOND_DERIV]=
		"Second order derivative with respect to age of the moment of "
		"inertia of the entire star in Msun Rsun^2/Gyr^2";
	__output_column_descr[OutCol::RRAD_SECOND_DERIV]=
		"First order derivative with respect to age of the radius of the "
		"stellar radiative core in solar radii per Gyr^2.";
}

void CommandLineOptions::init_track_column_descriptions_and_units()
{
	__track_column_descr.resize(
			CustomStellarEvolution::NUM_TRACK_QUANTITIES);
	__track_column_units.resize(
			CustomStellarEvolution::NUM_TRACK_QUANTITIES);

	__track_column_descr[CustomStellarEvolution::ICONV]=
		"Convective zone moment of inertia";
	__track_column_units[CustomStellarEvolution::ICONV]="Msun Rsun^2";

	__track_column_descr[CustomStellarEvolution::IRAD]=
		"Radiative zone moment of inertia";
	__track_column_units[CustomStellarEvolution::IRAD]="Msun Rsun^2 ";

	__track_column_descr[CustomStellarEvolution::RSTAR]="Stellar radius";
	__track_column_units[CustomStellarEvolution::RSTAR]="solar radii";

	__track_column_descr[CustomStellarEvolution::LSTAR]="Stellar luminosity";
	__track_column_units[CustomStellarEvolution::LSTAR]="solar luminosities";

	__track_column_descr[CustomStellarEvolution::RRAD]=
		"Radius of the radiative zone";
	__track_column_units[CustomStellarEvolution::RRAD]="solar radii.";

	__track_column_descr[CustomStellarEvolution::MRAD]=
		"Mass of the radiative zone";
	__track_column_units[CustomStellarEvolution::MRAD]="solar masses.";

	__track_column_descr[CustomStellarEvolution::AGE]="Stellar age";
	__track_column_units[CustomStellarEvolution::AGE]="Gyr";
}

void CommandLineOptions::init_defaults()
{
	__defaults.resize(InCol::NUM_REAL_INPUT_QUANTITIES);
	__defaults[InCol::LOW_MASS_WINDK]=0.17;
	__defaults[InCol::HIGH_MASS_WINDK]=0;
	__defaults[InCol::LOW_MASS_WIND_SAT_W]=2.45;
	__defaults[InCol::HIGH_MASS_WIND_SAT_W]=0;
	__defaults[InCol::LOW_MASS_WIND_SAT_P]=NaN;
	__defaults[InCol::HIGH_MASS_WIND_SAT_P]=Inf;
	__defaults[InCol::CORE_ENV_COUPLING_TIMESCALE]=5;
	__defaults[InCol::LGQ]=8;
	__defaults[InCol::LGQ_INERTIAL]=6;
	__defaults[InCol::MSTAR]=1;
	__defaults[InCol::MPLANET]=1;
	__defaults[InCol::RPLANET]=1;
	__defaults[InCol::PLANET_FORMATION_AGE]=0;
	__defaults[InCol::WDISK]=0.898;
	__defaults[InCol::PDISK]=NaN;
	__defaults[InCol::TDISK]=5;
	__defaults[InCol::A_FORMATION]=0.05;
	__defaults[InCol::P_FORMATION]=NaN;
	__defaults[InCol::INCLINATION_FORMATION]=0.0;
	__defaults[InCol::TSTART]=NaN;
	__defaults[InCol::TEND]=-Inf;
	__defaults[InCol::START_WRAD]=NaN;
	__defaults[InCol::START_WSURF]=NaN;
	__defaults[InCol::MAX_STEP]=Inf;
	__defaults[InCol::PRECISION]=6;

	__default_track_smoothing.resize(CustomStellarEvolution::AGE);
	__default_track_smoothing[CustomStellarEvolution::ICONV]=NaN;
	__default_track_smoothing[CustomStellarEvolution::IRAD]=NaN;
	__default_track_smoothing[CustomStellarEvolution::RSTAR]=NaN;
	__default_track_smoothing[CustomStellarEvolution::LSTAR]=NaN;
	__default_track_smoothing[CustomStellarEvolution::RRAD]=NaN;
	__default_track_smoothing[CustomStellarEvolution::MRAD]=NaN;

	__default_track_nodes.resize(CustomStellarEvolution::AGE);
	__default_track_nodes[CustomStellarEvolution::ICONV]=-1000;
	__default_track_nodes[CustomStellarEvolution::IRAD]=-1000;
	__default_track_nodes[CustomStellarEvolution::RSTAR]=-1000;
	__default_track_nodes[CustomStellarEvolution::LSTAR]=-1000;
	__default_track_nodes[CustomStellarEvolution::RRAD]=-1000;
	__default_track_nodes[CustomStellarEvolution::MRAD]=-1000;
}

///Fills in all the static members.
void CommandLineOptions::setup()
{
	init_input_column_names();
	init_output_column_descriptions();
	init_track_column_descriptions_and_units();
	init_defaults();
}

void CommandLineOptions::verify_custom_stellar_evolution()
{
	std::vector<bool> found_quantity(
			CustomStellarEvolution::NUM_TRACK_QUANTITIES);
	for(size_t i=0; i<__track_format.size(); ++i)
		found_quantity[__track_format[i]]=true;
	for(int i=0; i<CustomStellarEvolution::NUM_TRACK_QUANTITIES; ++i) {
		if(!found_quantity[i] && 
				i!=CustomStellarEvolution::LSTAR && 
				i!=CustomStellarEvolution::SKIP)
			throw Error::CommandLine(TRACK_COLUMN_NAMES[i]+" not found in "
					"the list of columns contained in the custom stellar "
					"evolution track.");
		if(i<CustomStellarEvolution::AGE) {
			double smoothing=__custom_track_smoothing[i]->dval[0];
			if(!std::isnan(smoothing) && (smoothing<-15 || smoothing>15))
				throw Error::CommandLine("Smoothing for "+
						TRACK_COLUMN_NAMES[i]+
						" is outside the allowed range of [-15; 15]!");
		}
	}
}

void init_output_column_names()
{
	OUTPUT_COLUMN_NAMES[OutCol::AGE]="t";
	OUTPUT_COLUMN_NAMES[OutCol::SEMIMAJOR]="a";
	OUTPUT_COLUMN_NAMES[OutCol::INCLINATION]="theta";
	OUTPUT_COLUMN_NAMES[OutCol::WORB]="Worb";
	OUTPUT_COLUMN_NAMES[OutCol::PORB]="Porb";
	OUTPUT_COLUMN_NAMES[OutCol::LORB]="Lorb";
	OUTPUT_COLUMN_NAMES[OutCol::LCONV]="Lconv";
	OUTPUT_COLUMN_NAMES[OutCol::LRAD]="Lrad";
	OUTPUT_COLUMN_NAMES[OutCol::LRAD_PAR]="Lradpar";
	OUTPUT_COLUMN_NAMES[OutCol::LRAD_PERP]="Lradperp";
	OUTPUT_COLUMN_NAMES[OutCol::LTOT]="L";
	OUTPUT_COLUMN_NAMES[OutCol::ICONV]="Iconv";
	OUTPUT_COLUMN_NAMES[OutCol::IRAD]="Irad";
	OUTPUT_COLUMN_NAMES[OutCol::ITOT]="I";
	OUTPUT_COLUMN_NAMES[OutCol::WSURF]="Wsurf";
	OUTPUT_COLUMN_NAMES[OutCol::WRAD]="Wrad";
	OUTPUT_COLUMN_NAMES[OutCol::PSURF]="Psurf";
	OUTPUT_COLUMN_NAMES[OutCol::PRAD]="Prad";
	OUTPUT_COLUMN_NAMES[OutCol::EVOL_MODE]="mode";
	OUTPUT_COLUMN_NAMES[OutCol::WIND_STATE]="wind";
	OUTPUT_COLUMN_NAMES[OutCol::RSTAR]="R";
	OUTPUT_COLUMN_NAMES[OutCol::LSTAR]="Lum";
	OUTPUT_COLUMN_NAMES[OutCol::RRAD]="Rrad";
	OUTPUT_COLUMN_NAMES[OutCol::MRAD]="Mrad";
	OUTPUT_COLUMN_NAMES[OutCol::ICONV_DERIV]="DIconv";
	OUTPUT_COLUMN_NAMES[OutCol::IRAD_DERIV]="DIrad";
	OUTPUT_COLUMN_NAMES[OutCol::ITOT_DERIV]="DI";
	OUTPUT_COLUMN_NAMES[OutCol::RSTAR_DERIV]="DR";
	OUTPUT_COLUMN_NAMES[OutCol::RRAD_DERIV]="DRrad";
	OUTPUT_COLUMN_NAMES[OutCol::MRAD_DERIV]="DMrad";
	OUTPUT_COLUMN_NAMES[OutCol::ICONV_SECOND_DERIV]="DDIconv";
	OUTPUT_COLUMN_NAMES[OutCol::IRAD_SECOND_DERIV]="DDIrad";
	OUTPUT_COLUMN_NAMES[OutCol::ITOT_SECOND_DERIV]="DDI";
	OUTPUT_COLUMN_NAMES[OutCol::RRAD_SECOND_DERIV]="DDRrad";
#ifdef COLUMN_NAME_EMPHASIS
	for(int i=0; i<OutCol::NUM_OUTPUT_QUANTITIES; ++i)
		OUTPUT_COLUMN_NAMES[i]=std::string(COLUMN_NAME_EMPHASIS) +
								OUTPUT_COLUMN_NAMES[i] +
								std::string(COLUMN_NAME_EMPHASIS);
#endif
}

void init_track_column_names()
{
	TRACK_COLUMN_NAMES[CustomStellarEvolution::AGE]="t";
	TRACK_COLUMN_NAMES[CustomStellarEvolution::ICONV]="Iconv";
	TRACK_COLUMN_NAMES[CustomStellarEvolution::IRAD]="Irad";
	TRACK_COLUMN_NAMES[CustomStellarEvolution::RSTAR]="R";
	TRACK_COLUMN_NAMES[CustomStellarEvolution::LSTAR]="Lum";
	TRACK_COLUMN_NAMES[CustomStellarEvolution::RRAD]="Rrad";
	TRACK_COLUMN_NAMES[CustomStellarEvolution::MRAD]="Mrad";
#ifdef COLUMN_NAME_EMPHASIS
	for(int i=0; i<CustomStellarEvolution::NUM_TRACK_QUANTITIES; ++i)
		TRACK_COLUMN_NAMES[i]=std::string(COLUMN_NAME_EMPHASIS) +
								TRACK_COLUMN_NAMES[i] +
								std::string(COLUMN_NAME_EMPHASIS);
#endif
}

const std::string CommandLineOptions::__default_outfname="poet.evol",
	  CommandLineOptions::__default_serialized_evol="serialized_evolution",
	  CommandLineOptions::__default_output_columns=
	  	"t,a,theta,Lconv,Lradpar,Lradperp,L,Iconv,Irad,I,mode",
	  CommandLineOptions::__default_track_columns="t,R,Iconv,Irad,Rrad,Mrad";

char *CommandLineOptions::cstr_copy(const std::string &str)
{
	char *cstr=new char[str.length()+1];
	std::strcpy(cstr, str.c_str());
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
		<< ".";
	__direct_value_options[InCol::CORE_ENV_COUPLING_TIMESCALE]=arg_dbl0(NULL,
			"core-env-coupling-timescale", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "Log base 10 of the tidal quality factor of the star "
		"outside the inertial mode range. In --input-columns identified by "
		<< __input_column_names[InCol::LGQ] << "'. Default: "
		<< __defaults[InCol::LGQ] << ".";
	__direct_value_options[InCol::LGQ]=arg_dbl0(NULL, "lgQ", "<double>",
			cstr_copy(option_help));

	option_help.str("");
	option_help << "Log base 10 of the tidal quality factor of the star in "
		"the inertial mode range. In --input-columns identified by "
		<< __input_column_names[InCol::LGQ_INERTIAL] << "'. Default: "
		<< __defaults[InCol::LGQ_INERTIAL] << ".";
	__direct_value_options[InCol::LGQ_INERTIAL]=arg_dbl0(NULL, "lgQinr",
			"<double>", cstr_copy(option_help));

	option_help.str("");
	option_help << "Mass of the star in solar masses. In --input-columns "
		"identified by '"
		<< __input_column_names[InCol::MSTAR] << "'. Default: "
		<< __defaults[InCol::MSTAR];
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
	option_help << "The inclination which the planet first "
			"appears. In --input-columns identified by '"
			<< __input_column_names[InCol::INCLINATION_FORMATION]
			<< "'. Default: " << __defaults[InCol::INCLINATION_FORMATION]
			<< ".";
	__direct_value_options[InCol::INCLINATION_FORMATION]=arg_dbl0("i",
			"init-inclination", "<double>", cstr_copy(option_help));

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
	option_help << "The maximum end age for the evolution in Gyr. If a "
		"negative value is passed, the evolution stops at absolute value "
		" of this parameter of if the star moves off the main sequence "
		"In --input-columns identified by '"
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
	option_help << "A list of ages to include in the calculated evolution. "
		"This argument can be overwritten or appended to by the '"
		<< __input_column_names[InCol::REQUIRED_AGES] << "' column in the "
		"input file. If the entry in the file starts with a comma, the ages "
		"listed there are combined with the value specified by this "
		"argument, otherwise only the ages listed in the file are used.";
	__required_ages_option=arg_str0(NULL, "require-ages",
			"<comma separated list>", cstr_copy(option_help));

	option_help.str("");
	option_help << "If nothing is read from an input file, then this "
		"argument should be the name of the file to use to output the "
		"evolution. If at least one quantity is read from a list file, that "
		"file should also contain a column (named '"
		<< __input_column_names[InCol::OUT_FNAME]
		<< "') specifying the input file name"
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
	for(int i=0; i<OutCol::NUM_OUTPUT_QUANTITIES; ++i)
		option_help << "\t* " << OUTPUT_COLUMN_NAMES[i]
					<< ": " << __output_column_descr[i] << std::endl;
	option_help << "Default: " << __default_output_columns;
	__output_file_columns=arg_str0(NULL, "output-columns",
			"<comma separated list>", cstr_copy(option_help));

	__input_fname=arg_file0("i", "input", "<file>", "The file to read the "
			"parameters of the planet-star systems to calculate evolutions "
			"for. If omitted, standard input is used instead. Any lines "
			"beginning with '#' are ignored. The other lines should start "
			"with at least the number of columns specified in the "
			"--input-columns option with white space only between columns.");

	option_help.str("");
	option_help << "The file to read previously serialized stellar "
			"evolution from or write one if the file does not exist. "
			"Default: '" << data_directory() << __default_serialized_evol
			<< "'.";
	__serialized_stellar_evolution=arg_file0(NULL, "serialized-stellar-evol",
			"<file>", cstr_copy(option_help));

	__custom_stellar_evolution=arg_file0(NULL, "custom-stellar-evolution",
			"<file>", "A single stellar evolution track from which to "
			"construct the stellar evolution to use (assumed to apply to all"
			" input systems regardless of the stellar mass). It should "
			"contain the quantities specified by "
			"--custom-stellar-evolution-format as columns in the precise "
			"order specified there. If this option is used, all input stars "
			"are treated as low mass stars, since that way the high mass "
			"star behavior can be reproduced  by simply assigning the entire"
			"star to a surface  convective zone.");

	option_help.str("");
	option_help << "A comma separated list of the columns in the custom "
		"stellar evolution track specified by the --custom-stellar-evolution"
		" option. The recognized column names are:" << std::endl;
	for(int i=0; i<CustomStellarEvolution::NUM_TRACK_QUANTITIES; ++i)
		option_help << "\t* "  << TRACK_COLUMN_NAMES[i]
					<< ": " << __track_column_descr[i] << " in " 
					<< __track_column_units[i] << std::endl;
	__custom_stellar_evolution_format=arg_str0(NULL,
			"custom-stellar-evolution-format", "<col1,col2,...>",
			cstr_copy(option_help));
	for(int i=0; i<CustomStellarEvolution::AGE; ++i) {
		option_help.str("");
		option_help << "Smoothing to apply to the "
			<< __track_column_descr[i] << " when interpolating the custom "
			   "stellar evolution track. If no smoothing should be applied "
			   "use 'NaN'. Default: " << __default_track_smoothing[i]
			<< "This option is ignored, unless --custom-stellar-evolution is"
			   " specified.";
		int to_strip=
#ifdef COLUMN_NAME_EMPHASIS
			std::string(COLUMN_NAME_EMPHASIS).size();
#else 
			0;
#endif
		std::string colname=TRACK_COLUMN_NAMES[i].substr(to_strip,
				TRACK_COLUMN_NAMES[i].size()-2*to_strip);

		__custom_track_smoothing[i]=arg_dbl0(NULL,
				cstr_copy(colname+"-smoothing"), "<double>",
				cstr_copy(option_help));
		option_help.str("");
		option_help << "The number of nodes to use when smoothing the "
			<< __track_column_descr[i] << " of the custom stellar evolution "
			   "track. This value is ignored unless both "
			   "--custom-stellar-evolution is used and --" + 
			   (colname+"-smoothing") + " is not NaN. Negative values result"
			   " in using a number of nodes which is the smaller of the "
			   "absolute of the given value and three times the number of "
			   "age entries in the track. Default: " 
			   << __default_track_nodes[i] << ".";
		__custom_track_nodes[i]=arg_int0(NULL,
				cstr_copy(colname+"-nodes"),
				"<integer>", cstr_copy(option_help));
	}
}

void CommandLineOptions::set_defaults()
{
	for(int i=0; i<InCol::NUM_REAL_INPUT_QUANTITIES; ++i)
		__direct_value_options[i]->dval[0]=__defaults[i];
	__output_fname->filename[0]=__default_outfname.c_str();
	std::ostringstream default_serialized;
	default_serialized << data_directory() << __default_serialized_evol;
	__serialized_stellar_evolution->filename[0]=
		cstr_copy(default_serialized);
	__output_file_columns->sval[0]=__default_output_columns.c_str();
	__custom_stellar_evolution->filename[0]="";
	__custom_stellar_evolution_format->sval[0]=
		__default_track_columns.c_str();
	for(int i=0; i<CustomStellarEvolution::AGE; ++i) {
		__custom_track_smoothing[i]->dval[0]=__default_track_smoothing[i];
		__custom_track_nodes[i]->ival[0]=__default_track_nodes[i];
	}
}

template<typename COL_ID_TYPE>
void CommandLineOptions::parse_column_list(const char *columns_str,
		const std::vector<std::string> column_names, int num_column_names,
		std::vector<COL_ID_TYPE> &columns, bool allow_noname)
{
	std::istringstream instream(columns_str);
	columns.clear();
	while(!instream.eof()) {
		std::string colname;
		std::getline(instream, colname, ',');
		int i=0;
		if(allow_noname && colname=="")
			columns.push_back(static_cast<COL_ID_TYPE>(num_column_names));
		else {
			while(i<num_column_names && column_names[i]!=colname) ++i;
			if(i==num_column_names)
				throw Error::CommandLine("Unrecognized input column '"+colname+
						"'");
			columns.push_back(static_cast<COL_ID_TYPE>(i));
		}
	}
}

void parse_real_list(const char *values_str, std::list<double> &values)
{
	std::istringstream instream(values_str);
	if(values_str[0]==',') ++values_str;
	else values.clear();
	while(!instream.eof()) {
		std::string val_string;
		std::getline(instream, val_string, ',');
		std::istringstream val_stream(val_string);
		double v;
		val_stream >> v;
		if(val_stream) values.push_back(v);
		else throw Error::BadFunctionArguments("The entry '"+
				val_string+"' cannot be parsed to a real number");
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
	__need_orbit=false;
	for(size_t i=0; i<__output_file_format.size(); ++i)
		if(__output_file_format[i]>OutCol::LAST_NO_ORBIT) __need_orbit=true;
	parse_column_list(__custom_stellar_evolution_format->sval[0],
			TRACK_COLUMN_NAMES, CustomStellarEvolution::NUM_TRACK_QUANTITIES,
			__track_format);
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

	if(__required_ages_option->count>0) {
		parse_real_list(__required_ages_option->sval[0], __required_ages);
		__required_ages.sort();
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
	__custom_track_smoothing(CustomStellarEvolution::AGE),
	__custom_track_nodes(CustomStellarEvolution::AGE),
	__opened_stream(false)
{
	setup();
	define_options();
	arg_lit *help_option=arg_lit0("h", "help", "Print this help and exit.");
	arg_lit *doxyhelp_option=arg_lit0(NULL, "doxygen-help", "Print this help "
			"formatted suitably for Doxygen and exit.");
	struct arg_end *end = arg_end(100);
	for(int i=0; i<InCol::NUM_REAL_INPUT_QUANTITIES; i++)
		__argtable[i]=__direct_value_options[i];
	__argtable[InCol::OUT_FNAME]=__output_fname;
	__argtable[InCol::START_LOCKED]=__start_locked;
	__argtable[InCol::REQUIRED_AGES]=__required_ages_option;
	__argtable[InCol::NUM_INPUT_QUANTITIES]=__input_file_columns;
	__argtable[InCol::NUM_INPUT_QUANTITIES+1]=__output_file_columns;
	__argtable[InCol::NUM_INPUT_QUANTITIES+2]=__input_fname;
	__argtable[InCol::NUM_INPUT_QUANTITIES+3]=__serialized_stellar_evolution;
	__argtable[InCol::NUM_INPUT_QUANTITIES+4]=__custom_stellar_evolution;
	__argtable[InCol::NUM_INPUT_QUANTITIES+5]=__custom_stellar_evolution_format;
	for(int i=0; i<CustomStellarEvolution::AGE; ++i) {
		__argtable[InCol::NUM_INPUT_QUANTITIES+6+2*i]=
			__custom_track_smoothing[i];
		__argtable[InCol::NUM_INPUT_QUANTITIES+7+2*i]=
			__custom_track_nodes[i];
	}
	__argtable[InCol::NUM_INPUT_QUANTITIES+2*CustomStellarEvolution::AGE+6]=
		help_option;
	__argtable[InCol::NUM_INPUT_QUANTITIES+2*CustomStellarEvolution::AGE+7]=
		doxyhelp_option;
	__argtable[InCol::NUM_INPUT_QUANTITIES+2*CustomStellarEvolution::AGE+8]=
		end;

	if(arg_nullcheck(__argtable) != 0) {
		cleanup();
		throw Error::CommandLine("Failed to allocate argument table.");
	}
	set_defaults();
	int nerrors=arg_parse(argc, argv, __argtable);
	if(help_option->count>0 || nerrors>0) {
        printf("Usage: %s", "poet");
        arg_print_syntax(stdout, __argtable,"\n\n");
        arg_print_glossary(stdout, __argtable,"  %-25s %s\n");
		if(help_option->count==0)
			arg_print_errors(stdout, end, "poet");
		else exit(0);
		__parsed_ok=false;
		return;
	}
	if(doxyhelp_option->count>0) {
        printf("SubPixPhot [options]\n\n");
		printf("Supported Options\n");
		printf("-----------------\n");
		printf("<dl>");
        arg_print_glossary(stdout, __argtable,
				"<dt>%s</dt>\n <dd>%s</dd>\n\n");
		printf("</dl>");
		exit(0);
	}
	postprocess();
	if(__custom_stellar_evolution->count>0)
		verify_custom_stellar_evolution();
	__parsed_ok=true;
}

double CommandLineOptions::get_real_value(InCol::InputColumns quantity) const
{
	if(quantity<0 || quantity>=InCol::NUM_REAL_INPUT_QUANTITIES)
		throw Error::BadFunctionArguments("Unrecognized real quantity in "
				"CommandLineOptions::get_real_value().");
	return __direct_value_options[quantity]->dval[0];
}

double CommandLineOptions::custom_track_smoothing(
		CustomStellarEvolution::Columns column) const
{
	if(column==CustomStellarEvolution::AGE)
		throw Error::BadFunctionArguments(
				"Requesting the smoothing to apply to custom track ages!");
	else if(column>=CustomStellarEvolution::NUM_TRACK_QUANTITIES)
		throw Error::BadFunctionArguments("Requesting the smoothing for an "
				"unknown custom track quantity!");
	else return __custom_track_smoothing[column]->dval[0];
}

int CommandLineOptions::custom_track_nodes(
		CustomStellarEvolution::Columns column) const
{
	if(column==CustomStellarEvolution::AGE)
		throw Error::BadFunctionArguments("Requesting the number of "
				"smoothing nodes for custom track ages!");
	else if(column>=CustomStellarEvolution::NUM_TRACK_QUANTITIES)
		throw Error::BadFunctionArguments("Requesting the number of "
				"smoothing nodes for an unknown custom track quantity!");
	else return __custom_track_nodes[column]->ival[0];
}

///AU/\f$\mathrm{R}_\odot\f$.
const double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;

void output_solution(const OrbitSolver &solver, const StellarSystem &system,
		const std::string &filename,
		const std::vector<OutCol::OutputColumns> &output_file_format,
		double start_age, double end_age, double timestep,
		const std::list<double> &required_ages)
{

	std::list<double>::const_iterator
		age_i=solver.get_tabulated_var(AGE)->begin(),
		a_i=solver.get_tabulated_var(SEMIMAJOR)->begin(),
		inclination_i=solver.get_tabulated_var(INCLINATION)->begin(),
		Lconv_i=solver.get_tabulated_var(LCONV)->begin(),
		Lrad_parallel_i=solver.get_tabulated_var(LRAD_PAR)->begin(),
		Lrad_perpendicular_i=solver.get_tabulated_var(LRAD_PERP)->begin();
	const Star &star=system.get_star();
	const Planet &planet=system.get_planet();
	double mstar=star.mass(), mplanet=planet.mass();
	std::list<double>::const_iterator
		last_age=solver.get_tabulated_var(AGE)->end();
	std::list<double> age_list;
	if(age_i==last_age) {
		std::list<double>::const_iterator
			required_ages_iter=required_ages.begin();
		for(double age=start_age; age<=end_age; age+=timestep) {
			for(;required_ages_iter!=required_ages.end() &&
					*required_ages_iter<age; ++required_ages_iter)
				age_list.push_back(*required_ages_iter);
			age_list.push_back(age);
		}
		age_i=age_list.begin();
		last_age=age_list.end();
	}

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
				   Irad=star.moment_of_inertia(*age_i, radiative),
				   dIconv_dt=star.moment_of_inertia_deriv(*age_i,convective),
				   dIrad_dt=star.moment_of_inertia_deriv(*age_i, radiative),
				   d2Iconv_dt2=
					   star.moment_of_inertia_deriv(*age_i, convective, 2),
				   d2Irad_dt2=
					   star.moment_of_inertia_deriv(*age_i, radiative, 2),
				   Rstar=star.radius(*age_i),
				   Lrad=std::sqrt(std::pow(*Lrad_parallel_i, 2) +
						   std::pow(*Lrad_perpendicular_i, 2)),
				   Ltot=std::sqrt(std::pow(*Lconv_i + *Lrad_parallel_i, 2) +
						   std::pow(*Lrad_perpendicular_i, 2)),
				   semimajor=(*a_i)/AU_Rsun,
				   worb=planet.orbital_angular_velocity_semimajor(semimajor),
				   Lorb=orbital_angular_momentum(mstar, mplanet, semimajor,
						   0);
			outf << std::setw(25);
			switch(output_file_format[i]) {
				case OutCol::AGE : outf << *age_i; break;
				case OutCol::SEMIMAJOR : outf << semimajor; break;
				case OutCol::WORB : outf << worb; break;
				case OutCol::PORB : outf << 2.0*M_PI/worb; break;
				case OutCol::LORB : outf << Lorb; break;
				case OutCol::INCLINATION : outf << *inclination_i; break;
				case OutCol::LCONV : outf << *Lconv_i; break;
				case OutCol::LRAD : outf << Lrad; break;
				case OutCol::LRAD_PAR : outf << *Lrad_parallel_i; break;
				case OutCol::LRAD_PERP : outf << *Lrad_perpendicular_i;
										 break;
				case OutCol::LTOT : outf << Ltot; break;
				case OutCol::ICONV : outf << Iconv; break;
				case OutCol::IRAD : outf << Irad; break;
				case OutCol::ITOT : outf << Iconv+Irad; break;
				case OutCol::WSURF : outf << (*Lconv_i)/Iconv; break;
				case OutCol::WRAD : outf << Lrad/Irad; break;
				case OutCol::PSURF : outf << 2.0*M_PI*Iconv/(*Lconv_i);
									 break;
				case OutCol::PRAD : outf << 2.0*M_PI*Irad/Lrad; break;
				case OutCol::EVOL_MODE : outf << *mode_i; break;
				case OutCol::WIND_STATE : outf << *wind_i; break;
				case OutCol::RSTAR : outf << Rstar; break;
				case OutCol::LSTAR : outf << star.luminosity(*age_i);
									 break;
				case OutCol::RRAD : outf << star.rad_radius(*age_i);
									break;
				case OutCol::MRAD : outf << star.rad_mass(*age_i);
									break;
				case OutCol::ICONV_DERIV : outf << dIconv_dt; break;
				case OutCol::IRAD_DERIV : outf << dIrad_dt; break;
				case OutCol::ITOT_DERIV : outf << dIconv_dt+dIrad_dt; break;
				case OutCol::RSTAR_DERIV : outf << 
										   	star.logradius_deriv(*age_i)*
								   			Rstar; break;
				case OutCol::RRAD_DERIV : outf <<
										  	star.rad_radius_deriv(*age_i);
										  break;
				case OutCol::MRAD_DERIV : outf <<
										  	star.rad_mass_deriv(*age_i);
								  break;
				case OutCol::ICONV_SECOND_DERIV : outf << d2Iconv_dt2; break;
				case OutCol::IRAD_SECOND_DERIV : outf << d2Irad_dt2; break;
				case OutCol::ITOT_SECOND_DERIV : outf <<
												 	d2Iconv_dt2 + d2Irad_dt2;
												 break;
				case OutCol::RRAD_SECOND_DERIV : outf <<
												 star.rad_radius_deriv(
														 *age_i, 2); break;
				default : throw Error::BadFunctionArguments(
								  "Unrecognized output column encountered in"
								  " output_file_format in "
								  "poet.cpp:output_solution.");
			}
		}
		outf << std::endl;
		++age_i; ++a_i; ++inclination_i, ++Lconv_i; ++Lrad_parallel_i;
		++Lrad_perpendicular_i; ++mode_i; ++wind_i;
	}
	outf.close();
}

void calculate_evolution(const std::vector<double> &real_parameters,
		bool start_locked, const std::list<double> &required_ages,
		const StellarEvolution &stellar_evolution,
		const std::string &outfname,
		const std::vector<OutCol::OutputColumns> &output_file_format,
		bool need_orbit)
{
	Star star(real_parameters[InCol::MSTAR],
			std::pow(10.0, real_parameters[InCol::LGQ]),
			std::pow(10.0, real_parameters[InCol::LGQ_INERTIAL]),
			real_parameters[InCol::WINDK],
			real_parameters[InCol::WIND_SAT_W],
			real_parameters[InCol::CORE_ENV_COUPLING_TIMESCALE]*1e-3,
			0,
			real_parameters[InCol::WDISK],
			real_parameters[InCol::TDISK]*1e-3,
			stellar_evolution);
	double tstart=real_parameters[InCol::TSTART];
	if(star.is_low_mass() && 
			star.core_formation_age()>star.disk_dissipation_age() &&
			tstart<star.core_formation_age() && need_orbit)
		throw Error::Runtime("At present the case when the disk dissipates "
				"before the stellar core starts to form is not supported.");
	Planet planet(star, real_parameters[InCol::MPLANET],
			real_parameters[InCol::RPLANET],
			real_parameters[InCol::A_FORMATION]);
	StellarSystem system(star, planet);
	std::valarray<double> start_orbit(0.0, 1);
	EvolModeType start_evol_mode;
	SpinOrbitLockInfo start_star_lock;
	if(std::isnan(tstart) || tstart<star.disk_dissipation_age()) {
		start_evol_mode=LOCKED_TO_DISK;
		if(std::isnan(tstart) || tstart<star.core_formation_age()) {
			tstart=star.core_formation_age();
			start_orbit[0]=(star.is_low_mass() ?
					real_parameters[InCol::WDISK] : 0)*
				star.moment_of_inertia(tstart, radiative);
		} else {
			start_orbit[0]=(star.is_low_mass() ?
					real_parameters[InCol::START_WRAD] : 0)*
				star.moment_of_inertia(tstart, radiative);
		}
	} else if(start_locked) {
		throw Error::NotImplemented("Starting in locked orbit");
/*		start_evol_mode=LOCKED_TO_PLANET;
		start_orbit.resize(2);
		start_orbit[0]=real_parameters[InCol::A_FORMATION]*AU_Rsun;
		start_orbit[1]=(star.is_low_mass() ?
				real_parameters[InCol::START_WRAD]*
				star.moment_of_inertia(tstart, radiative) : 0);*/
	} else {
		if(planet.orbital_angular_velocity_semimajor(
					real_parameters[InCol::A_FORMATION])<
				real_parameters[InCol::START_WSURF])
			start_evol_mode=BINARY;
		else start_evol_mode=BINARY;
		start_orbit.resize(5);
		start_orbit[0]=
			std::pow(real_parameters[InCol::A_FORMATION]*AU_Rsun, 6.5);
		start_orbit[1]=real_parameters[InCol::INCLINATION_FORMATION];
		start_orbit[2]=real_parameters[InCol::START_WSURF]*
			star.moment_of_inertia(tstart, convective);
		start_orbit[3]=(star.is_low_mass() ?
				real_parameters[InCol::START_WRAD]*
				star.moment_of_inertia(tstart, radiative) : 0);
		start_orbit[4]=0.0; //Start the core aligned with the surface.
	}
	OrbitSolver solver(real_parameters[InCol::TEND],
			std::pow(10.0, -real_parameters[InCol::PRECISION]));
	if(need_orbit)
		solver(system, real_parameters[InCol::MAX_STEP],
				real_parameters[InCol::PLANET_FORMATION_AGE],
				real_parameters[InCol::A_FORMATION], 
				real_parameters[InCol::INCLINATION_FORMATION], tstart,
				start_evol_mode, start_star_lock, start_orbit,
				required_ages);
    double tend=(real_parameters[InCol::TEND]>0 ? 
					real_parameters[InCol::TEND] :
					std::min(-real_parameters[InCol::TEND], 
							 star.lifetime()));
	output_solution(solver, system, outfname, output_file_format,
			tstart, tend, real_parameters[InCol::MAX_STEP], required_ages);
}

std::string update_run_parameters(std::vector<double> &real_parameters,
		bool &start_locked, std::list<double> &required_ages,
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
		else if(quantity==InCol::REQUIRED_AGES) {
			line >> word;
			parse_real_list(word.c_str(), required_ages);
			required_ages.sort();
		} else if(quantity==InCol::SKIP) line >> word;
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
	std::list<double> required_ages=options.required_ages();
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
				required_ages=options.required_ages();
				outfname=update_run_parameters(real_parameters,
					start_locked, required_ages, options.input_file_format(),
					line_str, input_lineno);
				if(start_locked)
					throw Error::NotImplemented("Locked evolution");
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
		calculate_evolution(real_parameters, start_locked, required_ages,
				stellar_evolution, outfname, options.output_file_format(),
				options.need_orbit());
	}
}

///Returns a pointer to a newly allocated stellar evolution constructed
///according to options.
StellarEvolution *get_stellar_evolution(const CommandLineOptions &options)
{
	std::ifstream serialized_file(options.serialized_stellar_evolution());
	bool serialized_exists=serialized_file.good();
	serialized_file.close();

	StellarEvolution *stellar_evolution;
	if(!options.custom_stellar_evolution().empty()) {
		std::vector<double> smoothing(CustomStellarEvolution::AGE);
		std::vector<int> nodes(CustomStellarEvolution::AGE);
		for(int i=0; i<CustomStellarEvolution::AGE; ++i) {
			CustomStellarEvolution::Columns
				col=static_cast<CustomStellarEvolution::Columns>(i);
			smoothing[i]=options.custom_track_smoothing(col);
			nodes[i]=options.custom_track_nodes(col);
		}
		stellar_evolution=new CustomStellarEvolution::Evolution(
				options.custom_stellar_evolution(),
				options.custom_track_format(),
				smoothing, nodes);
	} else if(serialized_exists) {
		stellar_evolution=new StellarEvolution;
		stellar_evolution->load_state(
				options.serialized_stellar_evolution());
		return stellar_evolution;
	} else
		stellar_evolution=new YRECEvolution(data_directory()+"YREC", 0, 2.0,
				2.0);
	if(!serialized_exists) stellar_evolution->save_state(
			options.serialized_stellar_evolution());
	return stellar_evolution;
}

///Calculates a realistic evolution chosen to be comlicated.
int main(int argc, char **argv)
{
	try {
#ifdef DEBUG
		std::cerr.precision(16);
		std::cerr.setf(std::ios_base::scientific);
#endif
		init_output_column_names();
		init_track_column_names();
		CommandLineOptions options(argc, argv);
		if(!options) return 1;
		StellarEvolution *stellar_evolution=get_stellar_evolution(options);
		run(options, *stellar_evolution, options.input());
		delete stellar_evolution;
	} catch(Error::General &err) {
		std::cerr << err.what() << ": " << err.get_message() << std::endl;
		return 2;
	}
}
