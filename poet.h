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

///All command line options can be accessed through members.
class CommandLineOptions {
private:
	///\brief The wind strength for low mass stars in units of
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	arg_dbl *__low_mass_windK,

			///\brief The wind strength for low mass stars in units of
			/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
			/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
			*__high_mass_windK,

			///\brief The frequency at which the wind saturates for
			///low mass stars in units of rad/day.
			*__low_mass_wind_saturation,

			///\brief The frequency at which the wind saturates for
			///high mass stars in units of rad/day.
			*__high_mass_wind_saturation,

			///The timescale on which the core end envelope are coupled.
			*__core_env_coupling_timescale,

			///Lg(tidal quality factor of the star).
			*__lgQ,

			///Mass of the star in \f$M_\odot\f$.
			*__star_mass,

			///Mass of the planet in \f$M_\Jupiter\f$
			*__planet_mass,

			///Radius of the planet in \f$R_\Jupiter\f$
			*__planet_radius,

			///The age at which the planet forms.
			///
			///If it is smaller than the disk dissipation age, the planet
			///forms at the disk dissipation age.
			*__planet_formation_age,

			///\brief The spin at which the star is locked while the
			///disk is present in rad/day.
			*__disk_lock_frequency,

			///The age at which the disk dissipates in Gyr.
			*__disk_dissipation_age,

			///The semimajor axis at which the planet first appears in AU.
			*__planet_formation_semimajor,

			///The starting age for the evolution in Gyr.
			*__start_age,

			///The maximum end age for the evolution in Gyr.
			*__end_age,

			///\brief Initial spin of the stellar core in rad/day for low
			///mass stars.
			*__start_wrad,

			///\brief Initial spin of the stellar surface in rad/day.
			*__start_wsurf,
			
			///A limit to impose on the ODE timestep in Gyr.
			*__max_timestep,
			
			///The precision to require of the solution
			*__precision;
			
	///\brief The wind strength column for low mass stars in units of
	/// \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
	/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
	arg_int *__low_mass_windK_column,

			///\brief The wind strength column for high mass stars in units
			///of \f$\frac{M_\odot \cdot R_\odot^2 \cdot \mathrm{day}^2}
			/// {\mathrm{rad}^2\cdot\mathrm{Gyr}}\f$.
			*__high_mass_windK_column,

			///\brief The column of the frequency at which the wind
			///saturates for low mass stars in units of rad/day.
			*__low_mass_wind_saturation_column,

			///\brief The column of the frequency at which the wind
			///saturates for high mass stars in units of rad/day.
			*__high_mass_wind_saturation_column,

			///\brief The column of the timescale on which the core end
			///envelope are coupled.
			*__core_env_coupling_timescale_column,

			///The Lg(tidal quality factor of the star) column.
			*__lgQ_column,

			///The column of the mass of the star in \f$M_\odot\f$.
			*__star_mass_column,

			///The column of the mass of the planet in \f$M_\Jupiter\f$
			*__planet_mass_column,

			///The column of the radius of the planet in \f$R_\Jupiter\f$
			*__planet_radius_column,

			///The column of the age at which the planet forms.
			///
			///If it is smaller than the disk dissipation age, the planet
			///forms at the disk dissipation age.
			*__planet_formation_age_column,

			///\brief The column of the spin at which the star is locked
			///while the disk is present in rad/day.
			*__disk_lock_frequency_column,

			///The column of the age at which the disk dissipates in Gyr.
			*__disk_dissipation_age_column,

			///The column of the semimajor axis at which the planet first
			///appears in AU.
			*__planet_formation_semimajor_column,

			///The column of the starting age for the evolution in Gyr.
			*__start_age_column,

			///The column of the maximum end age for the evolution in Gyr.
			*__end_age_column,

			///\brief The column of the initial spin of the stellar core in
			///rad/day for low mass stars.
			*__start_wrad_column,

			///\brief The column of the initial spin of the stellar surface
			///in rad/day.
			*__start_wsurf_column;

	///The name of the file to read the evolution scenarios from.
	arg_file *__input_fname,

			 ///The name of the file to output the solution to.
			 *__output_fname,

			 ///\brief The name of the file to read pre-serialized stellar
			 ///evolution from.
			 *__serialized_stellar_evolution;

	///Whether the planet should start locked to the star.
	///
	///If true, causes the value of __start_wsurf to be ignored.
	arg_lit *__start_locked;

	///The stream to the input filename if stdin is not being used.
	std::ifstream __input_stream;

	///\brief The column number in the input file containing the name of the
	///output filename.
	///
	///Not initialized if no quantities are read from a list file.
	int __output_fname_column;

	///Defines the command line options.
	void define_options();

	///Sets default values to the appropriate options.
	void set_defaults();

	///\brief Updates some of the values of the command line options after
	///parsing.
	///
	///On the command line, column indices are counted from 1, while it is
	///more convenient to count from 0 in the code.
	///
	///Also opens the input file if it is being used.
	///
	///Converst the precision required from number of significant
	///figures to an actual value.
	///
	///Converts the various ages that are specified in Myrs to Gyrs.
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
