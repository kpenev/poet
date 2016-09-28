#ifndef __POET_H
#define __POET_H

/**\file 
 *
 * \brief Defines the command line options class for the main
 * executable.
 *
 * \todo Add command line option to chose between YREC and MESA
 */

#include "Common.h"
#include "AstronomicalConstants.h"
#include "OrbitSolver.h"
#include "YRECIO.h"
#include "YRECStar.h"
#include "LockedPlanet.h"
#include "CustomStellarEvolution.h"
#include "DiskPlanetSystem.h"
#include "ConstSolutionIterator.h"
#include "IOColumns.h"
#include <argtable2.h>
#include <iostream>
#include <fstream>
#include <sstream>

std::vector<std::string>
	TRACK_COLUMN_NAMES(CustomStellarEvolution::NUM_TRACK_QUANTITIES);

std::string data_directory()
{
	std::string datadir=DATADIR;
	if(datadir[datadir.size()-1]!='/') datadir+='/';
	return datadir;
}

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
				 __output_column_descr,
				 
				 ///\brief Description of the columns in the input stellar
				 ///evolution track.
				 ///
				 ///Indexed by the corresponding
				 ///CustomStellarEvolution::Columns tag.
				 __track_column_descr,
				 
				 ///\brief Descriptions of the units expected of the stellar
				 ///evolution track quantities.
				 __track_column_units;

	///\brief The default values for the quantities defining the evolution to
	///calculate.
	std::vector<double> __defaults;

	///\brief Default smoothing to apply to custom stellar evolution track
	///quantities.
	///
	///Smoothing value of NaN signifies no smoothing.
	std::vector<double> __default_track_smoothing;

	///\brief Default numbef or spline nodes for custom stellar evolution
	///track quantities.
	///
	///For non-smoothed quantities (NaN entries in __default_track_smoothing)
	///this value is ignored.
	std::vector<int> __default_track_nodes;

	///The default output filename.
	static const std::string __default_outfname,

		///\brief The default filename to read the serialized stellar
		///evoliton from.
		__default_serialized_evol,
		
		///The default output columns
		__default_output_columns,
		
		///The default columns in a custom stellar evolution track.
		__default_track_columns,
		
		__default_eccentricity_expansion;

	///The command line options which directly specify a value.
	std::vector<arg_dbl*> __direct_value_options;

	///\brief The command line options which specify the smoothing for custom
	///stellar evolution tracks.
	std::vector<arg_dbl*> __custom_track_smoothing;

	///\brief The command line options which specify the number of nodes to
	///use if a custom track quantity is to be smoothed.
	std::vector<arg_int*> __custom_track_nodes;

	///The columns in the input file.
	arg_str *__input_file_columns,

			///The columns to write to the output file.
			*__output_file_columns,
			
			///\brief The comma separated list of ages to include in the
			///tabulated evolution.
			*__required_ages_option,
			
			*__custom_stellar_evolution_format;

	///Whether the initial surfarce spin of the star should be
	///synchronous to the orbit.
	arg_lit *__start_locked;

	///The name of the file to read the evolution scenarios from.
	arg_file *__input_fname,

			 ///The name of the file to output the solution to.
			 *__output_fname,

			 ///\brief The name of the file to read pre-serialized stellar
			 ///evolution from.
			 *__serialized_stellar_evolution,
			 
			 ///The filename of the custom stellar evolution track.
			 *__custom_stellar_evolution,
			 
			 ///\brief The name of the file to read eccentricity expansion
			 ///coefficients from.
			 *__eccentricity_expansion;

	void *__argtable[InCol::NUM_INPUT_QUANTITIES
                     +
                     2 * CustomStellarEvolution::AGE
                     +
                     10];

	///A list of the columns in the input file.
	std::vector<InCol::InputColumns> __input_file_format;

	///A list of the columns in the output file.
	std::vector<OutCol::OutputColumns> __output_file_format;

	///A list of the columns in the custom stellar evolution track.
	std::vector<CustomStellarEvolution::Columns> __track_format;

	///The stream to the input filename if stdin is not being used.
	std::ifstream __input_stream;

	///The copies of the option help strings made when creating the options.
	std::list<char*> __option_help_copies;

	///\brief The sorted list of ages to make sure are included in the
	///tabulated evolution.
	std::list<double> __required_ages;

	///Did parsing the command line succeed.
	bool __parsed_ok,

		 ///Set to true only if and when the #__input_stream is opened.
		 __opened_stream,
		 
		 ///Do we need to calculate the evolution in order to have all
		 ///required output quantities.
		 __need_orbit;

	///Returns a copy of the c-string content of the stream.
	char *cstr_copy(const std::ostringstream &stream)
	{return cstr_copy(stream.str());}

	///Returns a copy of the c-string content of the stream.
	char *cstr_copy(const std::string &str);

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

	///Fills in the names of the input columns in __input_column_names
	void init_input_column_names();

	///\brief Fills is the descriptions of the output columns in
	///__output_column_descr.
	void init_output_column_descriptions();

	///\brief Fills is the descriptions of the columns in a custom stellar
	///evolution track in __track_column_descr.
	void init_track_column_descriptions_and_units();

	///\brief Fills in default values for all possible real valued input
	///quantities and the smoothing parameters for custom stellar evolution.
	void init_defaults();

	///Fills in all the static members.
	void setup();

	///Applies common sense checks that all custom stellar evolution options.
	///
	///Throws an exception if something is wrong.
	void verify_custom_stellar_evolution();
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

	///\brief The name of the file to read/write a serialized stellar
	///evolution from/to.
	const char *serialized_stellar_evolution() const
	{return __serialized_stellar_evolution->filename[0];}

	///\brief The name of the file to read eccentricity expansion
	///coefficients from.
	const char *eccentricity_expansion() const
	{return __eccentricity_expansion->filename[0];}

	///\brief Whether the planet should start locked to the star.
	///
	///If true, causes the value of __start_wsurf to be ignored.
	bool start_locked() const {return __start_locked->count>0;}

	///Are any quantities to be read from a list file?
	bool input_from_list() const {return __input_file_format.size();}

	///List of the columns expected in the input file.
	const std::vector<InCol::InputColumns> &input_file_format() const
	{return __input_file_format;}

	///The columns to output.
	const std::vector<OutCol::OutputColumns> &output_file_format() const
	{return __output_file_format;}

	///Ages at which the evolution should definitely step.
	const std::list<double> &required_ages() const
	{return __required_ages;}

	///The filename from which to read a custom stellar evolution track.
	///
	///Empty string designates that a default evolution should be used.
	const std::string custom_stellar_evolution() const
	{return __custom_stellar_evolution->filename[0];}

	///A list of the columns in the custom stellar evolution track.
	const std::vector<CustomStellarEvolution::Columns>&
		custom_track_format() const {return __track_format;}

	///\brief The smoothing to apply to the given column from the custom 
	///stellar evolution track.
	double custom_track_smoothing(CustomStellarEvolution::Columns column)
		const;

	///\brief The nodes to use for the given column from the custom stellar
	///evolution track.
	int custom_track_nodes(CustomStellarEvolution::Columns column)
		const;

	///Do we need to calculate the evolution in order to have all required
	///output quantities.
	bool need_orbit() const {return __need_orbit;}

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
		const BinarySystem &system,

		///The name of the file to output the solution to.
		const std::string &filename,

		///The columns to include in the output file in the desired order.
		const std::vector<OutCol::OutputColumns> &output_file_format,

		///The starting age if no orbit was calculated. Ignored if solver 
		///contains an orbit.
		double start_age,

		///The starting age if no orbit was calculated. Ignored if solver 
		///contains an orbit.
		double end_age,
		
		///The time step if no orbit was calculated. Ignored if solver 
		///contains an orbit.
		double timestep,
		
		///A list of ages for which an output line must be written. Ignored
		///if solver contains an orbit
		const std::list<double> &required_ages=std::list<double>());

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
		const std::vector<OutCol::OutputColumns> &output_file_format,

		///Is calculating the orbit requried by the output.
		bool need_orbit);

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
