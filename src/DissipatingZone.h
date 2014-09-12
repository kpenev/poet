#ifndef __DISSIPATING_ZONE_H

/**\file
 *
 * \brief Declares a class representing one zone of a body dissipative to
 * tidal distortions.
 */

#include "EccentricityExpansionCoefficients.h"
#include "DissipationQuantities.h"
#include "Common.h"
#include <valarray>

///\brief A layer of a system body for which the tidal bulge is not exactly
///in phase with the tidal potential.
class DissipatingZone {
private:
	///\brief The constant coefficiients in \f$\mathcal{U}_{m,m'}\f$ of Lai
	///(2012).
	///
	///The first index is m+2 (since m starts from -2) and the second index
	///is m'/2+1 since the only allowed values are -2, 0 and 1.
	static const double __Umm_coef[][3],

				 ///\brief \f$\kappa_{m,m'}^+/\kappa_{m,m'}\f$ as a function
				 ///of \f$m=-2 \ldots 2\f$.
				 __torque_x_plus_coef[],

				 ///\brief \f$\kappa_{m,m'}^-/\kappa_{m,m'}\f$ as a function
				 ///of \f$m=-2 \ldots 2\f$.
				 __torque_x_minus_coef[];

	///The eccentricity expansion of \f$p_{m,s}\f$.
	static EccentricityExpansionCoefficients __pms;

	///The inclination with which __Ummp was last filled.
	double __Ummp_inclination;

	///The \f$\mathcal{U}_{m,m'}\f$ quantities defined in Lai (2012).
	std::valarray< std::valarray<double> > __Ummp,

		///The derivatives of the \f$\mathcal{U}_{m,m'}\f$ quantities w.r.t.
		///the inclination.
		__Ummp_deriv;

	///\brief The dimensionless tidal power and its derivatives.
	///
	///Consists of pairs of numbers one for each derivative. The first number
	///of each pair is always filled and if the zone is in a lock it is the
	///tidal power calculated assuming the zone spin frequency approaches the
	///lock from below. The second number is filled only if the zone is in a
	///spin-orbit lock and is the tidal power assuming the zone spin
	///frequency approaches the lock from above.
	std::valarray<double> __power,

		///\brief The dimensionless tidal torque in the x direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_x,

		///\brief The dimensionless tidal torque in the y direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_y,

		///\brief The dimensionless tidal torque in the y direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_y,

		///\brief The dimensionless tidal torque in the z direction and its
		///derivatives.
		///
		///See description of __power for details on the content.
		__torque_z;

	///Computes the \f$\mathcal{U}_{m,m'}\f$ values and their derivatives. 
	void fill_Umm();

	///\brief Calculates \f$\sum_s W_{2,s}D_{m,s}(\Theta)p_{s,m'}\f$ (see
	///documentation) and its derivatives w.r.t. e and \f$\Theta\f$.
	///
	///fill_Umm should already have been called with the appropriate
	///inclination.
	void potential_term(
			///The eccentricity.
			double e,

			///The m index.
			int m,

			///The m' index.
			int mp,
			
			///The highest order of the eccentricity to include.
			unsigned e_order,

			///Set to the undifferentiated value.
			double &no_deriv,

			///Set to the inclination derivative.
			double &inclination_deriv,

			///Set to the eccentricity_derivative.
			double &eccentricity_deriv);

protected:
	///To what order should eccentricity expansion be performed for the given
	///value of the eccentricity.
	virtual unsigned eccentricity_order(double e) const =0;

public:
	DissipatingZone();

	///\brief Defines the current orbit, triggering re-calculation of all
	///quantities.
	void set_orbit(
			///The angular velocity of the orbit in rad/day.
			double orbital_frequency,

			///The eccentricity of the orbit
			double eccentricity);

	///\brief Last setting for the angle between the angular momenta of the
	///zone and the orbit.
	virtual double inclination() const =0;

	///\brief Should return true iff the given term is presently locked.
	virtual bool locked(int orbital_frequency_multiplier,
			int spin_frequency_multiplier) const =0;

	///\brief Should return the tidal phase lag time the love number for the
	///given tidal term (or one of its derivatives).
	///
	///In case the specified term is in a lock, it should return the phase
	///lag for the case of the spin frequency approaching the lock from
	///below. The lag for spin frequency approaching from above should be
	///written to above_lock_value. If the term is not locked 
	///leave above_lock_value untouched.
	virtual double modified_phase_lag(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,
			
			///The current orbital spin frequency in rad/day
			double orbital_frequency,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::Derivative deriv,

			///If the lag of a locked term is calculated this should be set
			///to the lag assuming the spin frequency is just above the lock.
			///Otherwise, leave untouched.
			double &above_lock_value) const =0;

	///\brief Should return the corresponding component of the love
	///coefficient (Lai 2012 Equation 24).
	virtual double love_coefficient(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::Derivative deriv) const =0;

	///\brief The dimensionless tidal power or one of its derivatives.
	double tidal_power(
			///What to return
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,

			///If a spin-orbit lock is in effect and the time-lag is
			///discontinuous near zero forcing frequency, two possible values
			///can be calculated, assuming that the zone spin frequency
			///approaches the lock from below (default) or from above.
			bool above=false) const
	{return __power[2*deriv+(above? 1 : 0)];}

	///\brief The dimensionless tidal torque along x.
	///
	///See tidal_power() for a description of the arguments.
	double tidal_torque_x(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			bool above=false) const
	{return __torque_x[2*deriv+(above? 1 : 0)];}

	///\brief The dimensionless torque along y.
	///
	///See tidal_power() for a description of the arguments.
	double tidal_torque_y(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			bool above=false) const
	{return __torque_y[2*deriv+(above? 1 : 0)];}

	///\brief The dimensionless tidal torque along z.
	///
	///See tidal_power() for a description of the arguments.
	double tidal_torque_z(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			bool above=false) const
	{return __torque_z[2*deriv+(above? 1 : 0)];}

	///\brief Reads the eccentricity expansion coefficients of \f$p_{m,s}\f$.
	///
	///The given file should have been generated by
	///tabulate_eccentricity_expansion_coefficients.py.
	static void read_eccentricity_expansion(const std::string &fname)
	{__pms.read(fname);}
};

#endif
