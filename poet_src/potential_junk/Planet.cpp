/**\file
 *
 * \brief The definition of some of the methods of the Planet class.
 *
 * \ingroup StellarSystem_group
 */

#include "Planet.h"
#include "AstronomicalConstants.h"
#include <cmath>
#include <limits>

Planet::Planet(const Star &star, double mass, double radius,double semimajor)
: DissipatingBody(0, 0),
	__star(star), __semimajor(NULL),
	__mass(mass*AstroConst::jupiter_mass/AstroConst::solar_mass),
	__radius(radius*AstroConst::jupiter_radius/AstroConst::solar_radius),
	__current_semimajor(semimajor)
{
	__mstar=star.mass()*AstroConst::solar_mass;
	__rplanet=radius*AstroConst::jupiter_radius;
	__mplanet=mass*AstroConst::jupiter_mass;
	if(__rplanet==0) __rroche=0;
	else __rroche=(__mplanet>0 ? 2.44*__rplanet*std::pow(__mstar/__mplanet,
				1.0/3.0) : Inf);
}

double Planet::semimajor_derivative(double age) const
{
	const CubicSplineDerivatives *deriv=__semimajor->deriv(age);
	double result=deriv->order(1); 
	delete deriv; 
	return result;
}

double Planet::orbital_angular_velocity_semimajor(double semimajor) const
{
	return std::sqrt(AstroConst::G*(__mstar+__mplanet)/
			std::pow(AstroConst::AU*semimajor, 3.0))*
			AstroConst::day;
}

double Planet::orbital_angular_velocity_semimajor_deriv(double semimajor)
	const
{
	return -1.5*orbital_angular_velocity_semimajor(semimajor)/semimajor;
}

void Planet::set_semimajor_evolution(const std::valarray<double> &ages,
		const std::valarray<double> &semimajor_values,
		const std::valarray<double> &semimajor_derivatives)
{
	__semimajor=new InterpolatingFunctionALGLIB(ages, 
		semimajor_values, semimajor_derivatives);
	if(semimajor_values[semimajor_values.size()-1] < 0) {
		__lifetime=*(__semimajor->crossings());
	}
	else __lifetime=__star.lifetime();
}

double Planet::minimum_semimajor(double age) const
{
	return std::max(
		__star.radius(age)*AstroConst::solar_radius,
		__rroche)/AstroConst::AU;
}

double Planet::orbital_angular_momentum(double semimajor) const
{
	using namespace std;
	using namespace AstroConst;
	double a = AU*semimajor;
	return sqrt(G*a/(__mplanet+__mstar))*__mplanet*__mstar/
		(solar_mass*std::pow(solar_radius,2)/day);
}

double Planet::orbital_angmom_deriv(double semimajor, 
		double semimajor_deriv) const
{
	double a=semimajor*AstroConst::AU;
	return 0.5*__mplanet*__mstar*std::sqrt(AstroConst::G/
			((__mplanet+__mstar)*a))*
		semimajor_deriv*AstroConst::AU*AstroConst::day/
		(AstroConst::solar_mass*std::pow(AstroConst::solar_radius,2));
}


double Planet::orbital_angmom_deriv(double age) const
{
	return orbital_angmom_deriv((*__semimajor)(age), 
			semimajor_derivative(age));
}

