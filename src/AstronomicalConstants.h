/**\file
 * 
 * \brief Defines various astronomical constants
 *
 * \ingroup Utilities_group
 */

#ifndef __ASTRONOMICAL_CONSTANTS_H
#define __ASTRONOMICAL_CONSTANTS_H

/**\brief Namespace to isolate physical constants in.
 *
 * \ingroup Utilities_group
 */
namespace AstroConst {
	const double 
		AU=149597870700.0, ///< Astronomical unit [m]
		solar_mass=1.98892e30, ///< Mass of the sun [kg]
		solar_radius=6.955e8, ///< Radius of the sun [m]
		solar_age=4.57, ///< Age of the sun in [Gyr]
		Gyr=1e9*24.0*3600.0*365.25636, ///< Gyr [s]
		day=24.0*3600.0, ///< day [s]
		G=6.674e-11, ///< Gravitational constant in SI
		jupiter_mass=1.8987e27, ///< Mass of Jupiter [kg]
		jupiter_radius=7.1492e7; ///< Radius of Jupiter [m]
};

#endif
