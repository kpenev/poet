#include "DissipatingBody.h"
#include "AstronomicalConstants.h"
#include "Common.h"
#include <cmath>

enum DissipationComponent {POWER, TORQUE_X, TORQUE_Z};

///\brief The rates of change of various quantities due to tidal dissipation.
class TidalDissipation {
private:
	///\brief The constant coefficiients in \f$\mathcal{U}_{m,m'}\f$ of Lai
	///(2012).
	///
	///The first index is m+2 (since m starts from -2) and the second index
	///is m'/2+1 since the only allowed values are -2, 0 and 1.
	static const double __Umm1_coef[][],
				 __torque_x_plus_coef[],
				 __torque_x_minus_coef[];

	///\brief The tidal dissipation power in body1 in units of 
	/// \f$\frac{M_\odot R_\odot^2 rad^2}{day^2\,Gyr}\f$
	double __power1,

		   ///\brief The tidal dissipation power in body2 in units of 
		   /// \f$\frac{M_\odot R_\odot^2 rad^2}{day^2\,Gyr}\f$
		   __power2

		   ///\brief The x tidal dissipation torque of body 1 in units of
		   /// \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
		   __torque1_x,

		   ///\brief The x tidal dissipation torque of body 2 in units of
		   /// \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
		   __torque2_x,

		   ///\brief The z tidal dissipation torque of body 1 in units of
		   /// \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
		   __torque1_z
		  
		   ///\brief The z tidal dissipation torque of body 2 in units of
		   /// \f$\frac{M_\odot R_\odot^2 rad}{day\,Gyr}\f$.
		   __torque2_z;

	///Computes the \f$\mathcal{U}_{m,m'}\f$ values. 
	void fill_Umm(
			///The angle between the spin angular momentum of the body and
			///the orbital angular momentum in radians.
			double inclination,
			
			std::valarray< std::valarray<double> > &Umm);
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

