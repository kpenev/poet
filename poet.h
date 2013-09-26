#ifndef __POET_H
#define __POET_H

/**\file 
 *
 * \brief Defines the command line options class for the main
 * executable.
 */

#include <argtable2.h>

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

			///\brief The spin at which the star is locked while the
			///disk is present in rad/day.
			*__disk_lock_frequency,

			///The age at which the disk dissipates in Gyr.
			*__disk_dissipation_age,

			///The semimajor axis at which the planet first appears in AU.
			*__planet_formation_semimajor,

			///The starting age for the evolution in Gyr.
			*__start_age,

			///\brief Initial spin of the stellar core in rad/day for low
			///mass stars.
			*__start_wrad,

			///\brief Initial spin of the stellar surface in rad/day.
			*__start_wsurf;

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

			///\brief The column of the initial spin of the stellar core in
			///rad/day for low mass stars.
			*__start_wrad_column,

			///\brief The column of the initial spin of the stellar surface
			///in rad/day.
			*__start_wsurf_column;

	///Did parsing the command line succeed.
	bool __parsed_ok;
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
	double low_mass_wind_saturatior() const
	{return __low_mass_wind_saturation->dval[0];}

	///\brief The frequency at which the wind saturates for
	///high mass stars in units of rad/day.
	double high_mass_wind_saturatior() const
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

	///Initial spin of the stellar core in rad/day for low mass stars.
	double start_wrad() const {return __start_wrad->dval[0];}

	///Initial spin of the stellar surface in rad/day.
	double start_wsurf() const {return __start_wsurf->dval[0];}

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

	///\brief The column of the initial spin of the stellar core in
	///rad/day for low mass stars.
	int start_wrad_column() const {return __start_wrad_column->ival[0];}

	///\brief The column of the initial spin of the stellar surface
	///in rad/day.
	int start_wsurf_column() const {return __start_wsurf_column->ival[0];}

	///Did parsing the command line succeed.
	operator bool() {return __parsed_ok;}
};

#endif
