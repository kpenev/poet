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
		MPLANET, ///< Mass of the planet in \f$M_\Jupiter\f$.
		RPLANET, ///< Radius of the planet in \f$R_\Jupiter\f$.
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

		///\brief Should the evolution start with the stellar surface
		///spinning synchronously with the orbit? 
		START_LOCKED,	

		MAX_STEP,///< The maximum timestep to take.
		PRECISION,///< The precision to require of the evolution.
		OUT_FNAME,///< The name of the file to write the evolution to.

		///The number of real values quantities.
		NUM_REAL_INPUT_QUANTITIES=OUT_FNAME,

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

		///\brief Moment of inertia of the convective zone of the star (low mass
		///stars only)
		ICONV,

		///\brief Moment of inertia of the radiative zone of the star (low mass
		///stars only)
		IRAD,

		ITOT,///< Total moment of inertia of the star.
		WSURF, ///< Angular velocity of the stellar surface in rad/day.

		///Angular velocity of the stellar core in rad/day (low mass stars only).
		WRAD,

		PSURF, ///< Spin period of the stellar surface in days.
		PRAD,///< Spin period of the stellar core in days (low mass stars only).

		EVOL_MODE,///< The evolution mode for the step that starts at this age.

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

///All command line options can be accessed through members.
class CommandLineOptions {
private:

	///\brief The names to use for the input columns in the --input-columns
	///option (indexed by the corresponding InputColumns tag).
	static const std::string __input_column_names[],

				 ///\brief The names to use for the output columns in thei
				 ///--output-columns option (indexed by the corresponding
				 ///OutputColumns tag).
				 __output_column_names[],
				 
				 ///Description of the output columns,
				 __output_column_descr[];

	///\brief The default values for the quantities defining the evolution to
	///calculate.
	static const double __defaults[];

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
			*__output_file_columns;

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

	///Defines the command line options.
	void define_options();

	///Sets default values to the appropriate options.
	void set_defaults();

	///Parses a comma separated list of column names.
	template<typename COL_ID_TYPE>
	void parse_column_list(
			///The comma separated list of column names.
			char *columns_str,

			///The allowed column names 
			const std::string[] column_names,

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
	///Converts the various ages that are specified in Myrs to Gyrs.
	///
	///Calculates values of options which were specified using an
	///alternative.
	void postprocess();

	///Free all manually allocated memory and close open streams.
	void cleanup();

	///Did parsing the command line succeed.
	bool __parsed_ok,

		 ///Are any quantities to be read from a list file?
		 __input_from_list;
public:
	///Parse the command line.
	CommandLineOptions(int argc, char **argv);

	///\brief Returns the value of the quantity, if it is not overwritten by
	///the input list.
	double get_value(InCol::InputColumns quantity);

	///\brief The wind strength for low mass stars in units of
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	double low_mass_windK() const {return __low_mass_windK->dval[0];}

	///\brief The wind strength for low mass stars in units of
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	double high_mass_windK() const {return __high_mass_windK->dval[0];}

	///\brief The frequency at which the wind saturates for
	///low mass stars in units of rad/day.
	double low_mass_wind_saturation() const
	{return __low_mass_wind_saturation->dval[0];}

	///\brief The frequency at which the wind saturates for
	///high mass stars in units of rad/day.
	double high_mass_wind_saturation() const
	{return __high_mass_wind_saturation->dval[0];}

	///The timescale on which the core end envelope are coupled.
	double core_env_coupling_timescale() const
	{return __core_env_coupling_timescale->dval[0];}

	///Lg(tidal quality factor of the star).
	double lgQ() const {return __lgQ->dval[0];}

	///Mass of the star in \f$M_\odot\f$.
	double star_mass() const {return __star_mass->dval[0];}

	///Mass of the planet in \f$M_\Jupiter\f$
	double planet_mass() const {return __planet_mass->dval[0];}

	///Radius of the planet in \f$R_\Jupiter\f$
	double planet_radius() const {return __planet_radius->dval[0];}

	///\brief The age at which the planet forms.
	///
	///If it is smaller than the disk dissipation age, the planet
	///forms at the disk dissipation age.
	double planet_formation_age() const
	{return __planet_formation_age->dval[0];}

	///\brief The spin at which the star is locked while the
	///disk is present in rad/day.
	double disk_lock_frequency() const
	{return __disk_lock_frequency->dval[0];}

	///The age at which the disk dissipates in Gyr.
	double disk_dissipation_age() const 
	{return __disk_dissipation_age->dval[0];}

	///The semimajor axis at which the planet first appears in AU.
	double planet_formation_semimajor() const 
	{return __planet_formation_semimajor->dval[0];}

	///The starting age for the evolution in Gyr.
	double start_age() const {return __start_age->dval[0];}

	///The maximum end age for the evolution in Gyr.
	double end_age() const {return __end_age->dval[0];}

	///Initial spin of the stellar core in rad/day for low mass stars.
	double start_wrad() const {return __start_wrad->dval[0];}

	///Initial spin of the stellar surface in rad/day.
	double start_wsurf() const {return __start_wsurf->dval[0];}

	///A limit to impose on the ODE timestep
	double max_timestep() const {return __max_timestep->dval[0];}

	///The precision to require of the solution
	double precision() const {return __precision->dval[0];}

	///\brief The wind strength column for low mass stars in units of
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	int low_mass_windK_column() const 
	{return __low_mass_windK_column->ival[0];}

	///\brief The wind strength column for high mass stars in units
	///of \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	int high_mass_windK_column() const
	{return __high_mass_windK_column->ival[0];}

	///\brief The column of the frequency at which the wind
	///saturates for low mass stars in units of rad/day.
	int low_mass_wind_saturation_column() const
	{return __low_mass_wind_saturation_column->ival[0];}

	///\brief The column of the frequency at which the wind
	///saturates for high mass stars in units of rad/day.
	int high_mass_wind_saturation_column() const
	{return __high_mass_wind_saturation_column->ival[0];}

	///\brief The column of the timescale on which the core end
	///envelope are coupled.
	int core_env_coupling_timescale_column() const
	{return __core_env_coupling_timescale_column->ival[0];}

	///The Lg(tidal quality factor of the star) column.
	int lgQ_column() const {return __lgQ_column->ival[0];}

	///The column of the mass of the star in \f$M_\odot\f$.
	int star_mass_column() const {return __star_mass_column->ival[0];}

	///The column of the mass of the planet in \f$M_\Jupiter\f$
	int planet_mass_column() const {return __planet_mass_column->ival[0];}

	///The column of the radius of the planet in \f$R_\Jupiter\f$
	int planet_radius_column() const
	{return __planet_radius_column->ival[0];}

	///\brief The column of the age at which the planet forms in Gyr.
	///
	///If it is smaller than the disk dissipation age, the planet
	///forms at the disk dissipation age.
	int planet_formation_age_column() const
	{return __planet_formation_age_column->ival[0];}

	///\brief The column of the spin at which the star is locked
	///while the disk is present in rad/day.
	int disk_lock_frequency_column() const
	{return __disk_lock_frequency_column->ival[0];}

	///The column of the age at which the disk dissipates in Gyr.
	int disk_dissipation_age_column() const
	{return __disk_dissipation_age_column->ival[0];}

	///The column of the semimajor axis at which the planet first
	///appears in AU.
	int planet_formation_semimajor_column() const
	{return __planet_formation_semimajor_column->ival[0];}

	///The column of the starting age for the evolution in Gyr.
	int start_age_column() const {return __start_age_column->ival[0];}

	///The column of the maximum end age for the evolution in Gyr.
	int end_age_column() const {return __end_age_column->ival[0];}

	///\brief The column of the initial spin of the stellar core in
	///rad/day for low mass stars.
	int start_wrad_column() const {return __start_wrad_column->ival[0];}

	///\brief The column of the initial spin of the stellar surface
	///in rad/day.
	int start_wsurf_column() const {return __start_wsurf_column->ival[0];}

	///\brief The stream to read the parameters of the planet-star systems 
	///for which to calculate evolution.
	std::istream &input()
	{if(__input_fname->count) return __input_stream; else return std::cin;}

	///The name of the file to read in the various evolution scenarios.
	std::string input_filename() const {return __input_fname->filename[0];}

	///The name of the file to output the solution to.
	std::string output_filename() const {return __output_fname->filename[0];}

	///The index of the column in the input file containing the output
	///filename.
	int output_filename_column() const {return __output_fname_column;}

	///The name of the file to read pre-serialized stellar evolution from.
	const char *serialized_stellar_evolution() const
	{return __serialized_stellar_evolution->filename[0];}

	///\brief Whether the planet should start locked to the star.
	///
	///If true, causes the value of __start_wsurf to be ignored.
	bool start_locked() const {return __start_locked->count>0;}

	///Are any quantities to be read from a list file?
	bool input_from_list() const {return __input_from_list;}

	///Did parsing the command line succeed.
	operator bool() {return __parsed_ok;}

	///Closes the input  filename if it was opened.
	~CommandLineOptions() {cleanup();}
};

///Calculates the evolution for a system with the given parameters.
void calculate_evolution(
		///Mass of the star
		double Mstar,
		
		///Tidal quality factor of the star
		double Q,
		
		///\brief Magnetic wind strength in
		/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
		/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
		double Kwind,
		
		///Wind saturation frequency in rad/day
		double wsat, 

		///Core-envelope coupling timescare in Gyr
		double coupling_timescale,
		
		///Rotation rate of the stellar surface while locked to the
		///disk in rad/day.
		double wdisk,
		
		///The age at which the disk dissipates in Gyr.
		double tdisk,

		///The mass of the planet in \f$M_\Jupiter\f$.
		double Mplanet,
		
		///The radius of the planet in \f$R_\Jupiter\f$.
		double Rplanet,

		///The age at which the planet forms.
		///
		///If it is smaller than the disk dissipation age, the planet
		///forms at the disk dissipation age.
		double planet_formation_age,

		///The semiamjor axis at which the planet first appears in AU.
		double a_formation,

		///The initial spin of the radiative core for a low mass
		///star in rad/day. Ignored if the star is high mass, or if the
		//evolution starts before the radiative core has formed.
		double start_wrad,

		///The initial spin of the stellar surface in rad/day. Ignored if the
		///evolution starts before the disk has dissipated.
		double start_wsurf,

		///The minimum age to start the evolution from.
		///If this is smaller than both the disk dissipation age and
		///the age at which the stellar core forms, it is substituted
		///by the smaller of these two values.
		double tstart,

		///The maximum age to end the evolution on. It is replaced by
		///the star's lifetime if the latter is shorter.
		double tend,

		///The maximum time step the solver is allowed to take in Gyr.
		double max_time_step,

		///If true, and the start age is after the disk has
		///dissipated, the evolution is started with the planet locked
		///to the star, ignoring the value of start_wsurf.
		bool start_locked,
		
		///A stellar evolution interpolator.
		const StellarEvolution &stellar_evolution,

		///An orbit solver variable which is used to calculate the
		///evolution and as a result, on output contains it.
		OrbitSolver &solver);

///Outputs the solution calculated by the given solver.
void output_solution(
		///The solver that contains the solution. sover.operator() should
		///already have been called.
		const OrbitSolver &solver,
		
		///The planet-star system for which solution was derived.
		const StellarSystem &system,

		///The name of the file to output the solution to.
		const std::string &filename);

///Actually calculates the orbital evolutions.
void run(
		///All the configuration from the command line.
		const CommandLineOptions &options,

		///A fully functional stellar evolution interpolator.
		const StellarEvolution &stellar_evolution);
		

#endif
