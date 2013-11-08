#ifndef __POET_H
#define __POET_H

/**\file 
 *
 * \brief Defines the command line options class for the main
 * executable.
 */

#include "Common.h"
#include "AstronomicalConstants.h"
#include "OrbitSolver.h"
#include "YRECIO.h"
#include <argtable2.h>
#include <iostream>
#include <fstream>
#include <sstream>

///Isolates the tags for the input columns.
namespace InCol {
	///\brief Tags for the quantities required to fully specify the system to
	///evolve.
	enum InputColumns {
		///\brief Wind strength in 
		/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
		/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
		WINDK=0, 

		LOW_MASS_WINDK=WINDK,///< The wind strength for low mass stars
		HIGH_MASS_WINDK,///< The wind strength for high mass stars
		WIND_SAT_W, ///< Wind saturation frequency in rad/day.

		///Wind saturation frequency in rad/day for low mass stars.
		LOW_MASS_WIND_SAT_W=WIND_SAT_W,

		///Wind saturation frequency in rad/day for high mass stars.
		HIGH_MASS_WIND_SAT_W,

		WIND_SAT_P, ///< Wind saturation period in days.

		///Wind saturation period in days for low mass stars.
		LOW_MASS_WIND_SAT_P=WIND_SAT_P, 

		///Wind saturation period in days for high mass stars.
		HIGH_MASS_WIND_SAT_P,

		///\brief Core-envelope coupling timescale in Myr, ignore for high
		///mass stars.
		CORE_ENV_COUPLING_TIMESCALE,

		LGQ,///< Lg(Q*).
		MSTAR, ///< Mass of the star in \f$M_\odot\f$.
		MPLANET, ///< Mass of the planet in Jupiter masses.
		RPLANET, ///< Radius of the planet in Jupiter radii.
		PLANET_FORMATION_AGE,///< Age when planet appears in Gyr.
		WDISK,///< Stellar surface spin while disk is present in rad/day.

		///Stellar surface rotation period in days while disk is present.
		PDISK,
		
		TDISK, ///< Age in Myr when disk dissipates.
		A_FORMATION, ///< Semimajor axis at which the planet forms in AU.
		P_FORMATION, ///< Orbital period in days at which the planet forms.
		TSTART, ///< The minimum age to start evolution at in Gyr.
		TEND, ///< The maximum age to stop the evolution at in Gyr.

		///\brief Initial rotation of the radiative core in rad/day if the
		///evolution starts after the core has formed.
		START_WRAD, 

		///\brief Surface rotation of the star in rad/day if the evolution
		///starts after the disk has dissipated.
		START_WSURF, 

		MAX_STEP,///< The maximum timestep to take.

		///The number of significant figures require of the evolution.
		PRECISION,

		OUT_FNAME,///< The name of the file to write the evolution to.

		///The number of real values quantities.
		NUM_REAL_INPUT_QUANTITIES=OUT_FNAME,

		///\brief Should the evolution start with the stellar surface
		///spinning synchronously with the orbit? 
		START_LOCKED,	

		///A list of ages guaranteed to be included in the tabulated orbit.
		REQUIRED_AGES,

		SKIP,///< A column which is not needed to calculate the evolution.

		///The number of different input quantities supported.
		NUM_INPUT_QUANTITIES=SKIP
	};
};

///Isolates the tags for the output columns.
namespace OutCol {
	///Tags for the possible columns to output.
	enum OutputColumns {
		AGE,///< Age of the system in Gyr.
		SEMIMAJOR,///< Semimajor axis of the orbit in AU.
		WORB,///< The orbital frequency in rad/day.
		PORB,///< The orbital period days.

		///\brief Angular momentum of the convective zone of the star in
		/// \f$ M_\odot R_\odot^2 \mathrm{rad}/\mathrm{day}\f$ (low mass
		///stars only)
		LCONV,

		///\brief Angular momentum of the radiative zone of the star in
		// \f$ M_\odot R_\odot^2 \mathrm{rad}/\mathrm{day}\f$ (low mass
		///stars only)
		LRAD,

		///\brief Total angular momentum of the star in 
		/// \f$ M_\odot R_\odot^2 \mathrm{rad}/\mathrm{day}\f$.
		LTOT,

		///\brief Moment of inertia of the convective zone of the star (low
		///mass stars only)
		ICONV,

		///\brief Moment of inertia of the radiative zone of the star (low
		///mass stars only)
		IRAD,

		ITOT,///< Total moment of inertia of the star.
		WSURF, ///< Angular velocity of the stellar surface in rad/day.

		///\brief Angular velocity of the stellar core in rad/day (low mass
		///stars only).
		WRAD,

		PSURF, ///< Spin period of the stellar surface in days.

		///Spin period of the stellar core in days (low mass stars only).
		PRAD,

		///The evolution mode for the step that starts at this age.
		EVOL_MODE,

		WIND_STATE,///< The saturation state of the wind.
		RSTAR,///< Radius of the star in \f$R_\odot\f$.
		LSTAR,///< Luminosity of the star in \f$L_\odot\f$.

		///Radius of the stellar core in \f$R_\odot\f$ (low mass stars only).
		RRAD,

		///Mass of the stellar core in \f$M_\odot\f$ (low mass stars only).
		MRAD,

		///The number of different output quantities supported.
		NUM_OUTPUT_QUANTITIES
	};
};

///\brief The names to use for the output columns in the
///--output-columns option (indexed by the corresponding
///OutputColumns tag) as well as when labeling the output file.
std::vector<std::string> OUTPUT_COLUMN_NAMES(OutCol::NUM_OUTPUT_QUANTITIES);

///All command line options can be accessed through members.
class CommandLineOptions {
private:

	///\brief The names to use for the input columns in the
	///--input-columns option.
	///
	///Indexed by the corresponding InputColumns tag.
	std::vector<std::string> __input_column_names,

				 ///\brief Description of the output columns.
				 ///
				 ///Indexed by the corresponding OutputColumns tag.
				 __output_column_descr;

	///\brief The default values for the quantities defining the evolution to
	///calculate.
	std::vector<double> __defaults;

	///The default output filename.
	static const std::string __default_outfname,

		///\brief The default filename to read the serialized stellar
		///evoliton from.
		__default_serialized_evol,
		
		///The default output columns
		__default_output_columns;

	///The command line options which directly specify a value.
	std::vector<arg_dbl*> __direct_value_options;

	///The columns in the input file.
	arg_str *__input_file_columns,

			///The columns to write to the output file.
			*__output_file_columns,
			
			///\brief The comma separated list of ages to include in the
			///tabulated evolution.
			*__required_ages_option;

	///Whether the initial surfarce spin of the star should be
	///synchronous to the orbit.
	arg_lit *__start_locked;

	///The name of the file to read the evolution scenarios from.
	arg_file *__input_fname,

			 ///The name of the file to output the solution to.
			 *__output_fname,

			 ///\brief The name of the file to read pre-serialized stellar
			 ///evolution from.
			 *__serialized_stellar_evolution;

	void *__argtable[InCol::NUM_INPUT_QUANTITIES+6];

	///A list of the columns in the input file.
	std::vector<InCol::InputColumns> __input_file_format;

		///A list of the columns in the output file.
	std::vector<OutCol::OutputColumns> __output_file_format;

	///The stream to the input filename if stdin is not being used.
	std::ifstream __input_stream;

	///The copies of the option help strings made when creating the options.
	std::list<char*> __option_help_copies;

	///\brief The sorted list of ages to make sure are included in the
	///tabulated evolution.
	std::list<double> __required_ages;

	///Returns a copy of the c-string content of the stream.
	char *cstr_copy(const std::ostringstream &stream);

	///Defines the command line options.
	void define_options();

	///Sets default values to the appropriate options.
	void set_defaults();

	///Parses a comma separated list of column names.
	template<typename COL_ID_TYPE>
	void parse_column_list(
			///The comma separated list of column names.
			const char *columns_str,

			///The allowed column names 
			const std::vector<std::string> column_names,

			///The number of allowed column names
			int num_column_names,

			///This is updated to contain the columns in the correct
			///order.
			std::vector<COL_ID_TYPE> &columns,
			
			///Whether to allow zero length column names.
			bool allow_noname=false);

	///\brief Updates some command line options after parsing and
	///parses column lists.
	///
	///Opens the input file if it is being used.
	///
	///Converts the precision required from number of significant
	///figures to an actual value.
	///
	///Calculates values of options which were specified using an
	///alternative.
	void postprocess();

	///Free all manually allocated memory and close open streams.
	void cleanup();

	///Did parsing the command line succeed.
	bool __parsed_ok,

		 ///Set to true only if and when the #__input_stream is opened.
		 __opened_stream;

	///Fills in all the static members.
	void setup();
public:
	///Parse the command line.
	CommandLineOptions(int argc, char **argv);

	///\brief Returns the value of the quantity, if it is not overwritten by
	///the input list.
	double get_real_value(InCol::InputColumns quantity) const;

	///\brief The stream to read the parameters of the planet-star systems 
	///for which to calculate evolution.
	std::istream &input()
	{if(__input_fname->count) return __input_stream; else return std::cin;}

	///The name of the file to read in the various evolution scenarios.
	std::string input_filename() const {return __input_fname->filename[0];}

	///The name of the file to output the solution to.
	std::string output_filename() const {return __output_fname->filename[0];}

	///The name of the file to read pre-serialized stellar evolution from.
	const char *serialized_stellar_evolution() const
	{return __serialized_stellar_evolution->filename[0];}

	///\brief Whether the planet should start locked to the star.
	///
	///If true, causes the value of __start_wsurf to be ignored.
	bool start_locked() const {return __start_locked->count>0;}

	///Are any quantities to be read from a list file?
	bool input_from_list() const {return __input_file_format.size();}

	const std::vector<InCol::InputColumns> &input_file_format() const
	{return __input_file_format;}

	const std::vector<OutCol::OutputColumns> &output_file_format() const
	{return __output_file_format;}

	const std::list<double> &required_ages() const
	{return __required_ages;}

	///Did parsing the command line succeed.
	operator bool() {return __parsed_ok;}

	///Closes the input  filename if it was opened.
	~CommandLineOptions() {cleanup();}
};

///\brief Parses a comma separated list of real values.
///
///If the list starts with a comma, the values are appended to the output
///list, otherwise the output list is overwritten with the new values.
void parse_real_list(
		///The comma separated list of values.
		const char *values_str,

		///The values to append to or overwrite.
		std::list<double> &values);

///Outputs the solution calculated by the given solver.
void output_solution(
		///The solver that contains the solution. sover.operator() should
		///already have been called.
		const OrbitSolver &solver,
		
		///The planet-star system for which solution was derived.
		const StellarSystem &system,

		///The name of the file to output the solution to.
		const std::string &filename,

		///The columns to include in the output file in the desired order.
		const std::vector<OutCol::OutputColumns> &output_file_format);

///Calculates the evolution for a set of parameters.
void calculate_evolution(
		///The real-valued parameters describing the system indexed by
		/// #InCol::InputColumns.
		const std::vector<double> &real_parameters,
		
		///Whether to start the system with the stellar surface rotation
		///synchronized to the orbit.
		bool start_locked,

		///The ages required to be included in the evolution (sorted).
		const std::list<double> &required_ages,

		///The stellar evolution interpolator.
		const StellarEvolution &stellar_evolution,

		///The filename to write the evolution to.
		const std::string &outfname,

		///The columns to include in the output file in the desired order.
		const std::vector<OutCol::OutputColumns> &output_file_format);

///\brief Updates the evolution parameters as indicated on the next line of
///the input stream and returns filename to output the evolution to.
std::string update_run_parameters(
		///The array of real valued parameters to update, indexed by
		///InCol::InputColumns.
		std::vector<double> &real_parameters,

		///Whether to start the system with the stellar surface spin
		///synchronous with the orbit.
		bool &start_locked,

		///The ages required to be included in the evolution (sorted).
		std::list<double> &required_ages,

		///List of the columns to read from the input stream.
		const std::vector<InCol::InputColumns> &input_format,

		///The line from the input stream to read parameters from.
		std::istringstream &line,

		///The line number being processed (used only for error messages).
		size_t input_lineno);

///Actually calculates the orbital evolutions.
void run(
		///All the configuration from the command line.
		const CommandLineOptions &options,

		///A fully functional stellar evolution interpolator.
		const StellarEvolution &stellar_evolution);
		

#endif
