#include "DissipatingBody.h"
#include "AstronomicalConstants.h"
#include "Common.h"
#include <cmath>

enum DissipationComponent {POWER, TORQUEX, TORQUEZ,
	DPOWER_DINCLINATION, DTORQUEX_DINCLINATION, DTORQUEZ_DINCLINATION,
	NUM_DISSIPATION_QUANTITIES};

///\brief The rates of change of various quantities due to tidal dissipation.
class TidalDissipation {
private:
	///\brief The constant coefficiients in \f$\mathcal{U}_{m,m'}\f$ of Lai
	///(2012).
	///
	///The first index is m+2 (since m starts from -2) and the second index
	///is m'/2+1 since the only allowed values are -2, 0 and 1.
	static const double __Umm1_coef[][],

				 ///\brief \f$\kappa_{m,m'}^+/\kappa_{m,m'}\f$ as a function of
				 /// \f$m=-2 \ldots 2\f$.
				 __torque_x_plus_coef[],

				 ///\brief \f$\kappa_{m,m'}^-/\kappa_{m,m'}\f$ as a function of
				 /// \f$m=-2 \ldots 2\f$.
				 __torque_x_minus_coef[];

	///The orbital frequency in rad/day.
	double __orbital_frequency;

	///The \f$\mathcal{U}_{m,m'}\f$ quantities defined in Lai (2012).
	std::valarray< std::valarray<double> > __Umm;

	///\brief Rates of change of the various quantities and derivatives due
	///to tidal dissipation
	std::valarray<double> __dissipation_rates;

	///Computes the \f$\mathcal{U}_{m,m'}\f$ values. 
	void fill_Umm(
			///The angle between the spin angular momentum of the body and
			///the orbital angular momentum in radians.
			double inclination
			
			///Whether to use the derivatives w.r.t. \f$\Theta\f$ instead of
			///the values
			bool deriv);

	///\brief Calculates the dimensionless x and y torques and power due to
	///tidal dissipation.
	///
	///The \f$\mathcal{U}_{m,m'}\f$ constants must already be filled.
	void calculate_torque_power(
			///The body doing the dissipating.
			DissipatingBody &body,

			///If this values is not zero, the spin frequency of is /assumed
			//to approach the orbital frequency from below/above if
			///forcing_sign is +1/-1 instead, regardless of the currently set
			///spin frequency of the body. Result is undefined if the spin
			///frequency is not very close to the orbital frequency and yet
			///this value is non-zero.
			short forcing_sign,

			///Set to the tidal dissipation power in the same units as
			///power_scale.
			double &power,

			///Set to the tidal dissipation torque in the x direction, in the
			///same units as torque_scale.
			double &torque_x,

			///Set to the tidal dissipation torque in the z direction, in the
			///same units as torque_scale.
			double &torque_y);
public:
	///\brief Calculates the rates of change of various quantities due to
	///tidal dissipation.
	///
	///For now only works for zero eccentricity, throws and ecception
	///otherwise.
	TidalDissipation(
			///The first dissipating body.
			const DissipatingBody &body1,

			///The second dissipating body.
			const DissipatingBody &body2,

			///The semimajor axis of the orbit in AU.
			double semimajor,

			///The eccentricity of the orbit.
			double eccentricity,
				
			///If this values is not zero, the spin frequency of body1 is
			///assumed to approach the orbital frequency from below/above if
			///forcing_sign is +1/-1 instead, regardless of the currently set
			///spin frequency of body1. Result is undefined if the spin
			///frequency is not very close to the orbital frequency and yet
			///this value is non-zero.
			short forcing_sign1=0,

			///If this values is not zero, the spin frequency of body2 is
			///assumed to approach the orbital frequency from below/above if
			///forcing_sign is +/- instead, regardless of the currently set
			///spin frequency of body2. Result is undefined if the spin
			///frequency is not very close to the orbital frequency and yet
			///this value is non-zero.
			short forcing_sign2=0);

	///\brief Retuns the rate of change of some quantity due to tidal
	///dissipation.
	///
	///The units are:  
	///power: \f$\frac{M_\odot R_\odot^2 rad^2}{day^2\,Gyr}\f$
	///torque: \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
	double operator()(DissipationComponent component)
	{
		switch(component) {
			case POWER :: return __power;
			case TORQUE_X :: return __torque_x;
			case TORQUE_Z :: return __torque_z;
			default : throw Error::BadFunctionArguments("Urecognized tidal "
							  "dissipation component in "
							  "Dissipation::operator()!");
		}
	}
};

