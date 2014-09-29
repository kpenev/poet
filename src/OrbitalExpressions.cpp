/**\file
 *
 * \brief The definitions of the orbital expression functions.
 */

#include "OrbitalExpressions.h"

double orbital_angular_velocity(double m1, double m2, double semimajor,
		bool deriv)
{
	return (deriv ? -1.5 : 1.0)*std::sqrt(
				AstroConst::G*(m1+m2)*AstroConst::solar_mass/
				std::pow(semimajor*AstroConst::solar_radius, (deriv ? 5 :3))
			)*AstroConst::day;
}

double orbital_energy(double m1, double m2, double semimajor)
{
	return -m1*m2/(2.0*semimajor)*
			AstroConst::G*AstroConst::solar_mass*
			std::pow(AstroConst::day/AstroConst::solar_radius, 2)/
			AstroConst::AU;
}

double orbital_angular_momentum(double m1, double m2, double semimajor,
		double eccentricity)
{
	return m1*m2*std::sqrt(AstroConst::G*semimajor*
			(1.0-std::pow(eccentricity, 2))/(m1+m2)*AstroConst::AU*
			AstroConst::solar_mass)*AstroConst::day/
		std::pow(AstroConst::solar_radius, 2);
}
