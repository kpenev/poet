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
#include <assert.h>

Planet::Planet(Star* star,
		double observed_mass, double observed_radius,
		double observed_semimajor)
	:star(star), semimajor_axis(NULL), mass(observed_mass),
	radius(observed_radius), current_semimajor(observed_semimajor)
{
	mstar=star->get_mass()*AstroConst::solar_mass;
	rplanet=radius*AstroConst::jupiter_radius;
	mplanet=mass*AstroConst::jupiter_mass;
	if(rplanet==0) rroche=0;
	else rroche=(mplanet>0 ? 2.44*rplanet*std::pow(mstar/mplanet, 1.0/3.0) :
				 Inf);
}

void Planet::transform_into_earth() {
	mass = 0.003;
	mplanet = mass*AstroConst::jupiter_mass;
	rroche = 0;
	current_semimajor = 1;
}

double Planet::get_semimajor_derivative(double age) const
{
	const CubicSplineDerivatives *deriv=semimajor_axis->deriv(age);
	double result=deriv->order(1); 
	delete deriv; 
	return result;
}

double Planet::orbital_angular_velocity_semimajor(double semimajor) const
{
	return std::sqrt(AstroConst::G*(mstar+mplanet)/
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
	semimajor_axis=new InterpolatingFunctionALGLIB(ages, 
		semimajor_values, semimajor_derivatives);
	if(semimajor_values[semimajor_values.size()-1] < 0) {
		lifetime=*(semimajor_axis->crossings());
	}
	else lifetime=star->get_lifetime();
}

double Planet::minimum_semimajor(double age) const
{
	return std::max(
		star->get_radius(age)*AstroConst::solar_radius,
		rroche)/AstroConst::AU;
}

double Planet::tidal_decay(double age, double semimajor, 
		double stellar_spin_frequency, bool use_a6p5) const
{
	double tidal_w=(stellar_spin_frequency==Inf ||
			stellar_spin_frequency==-Inf ? -stellar_spin_frequency :
			orbital_angular_velocity_semimajor(semimajor)-
			stellar_spin_frequency),
		   Q=star->get_tidal_Q(tidal_w),
		   rstar=star->get_radius(age)*AstroConst::solar_radius;
	if(use_a6p5)
		return -4.5*6.5*std::sqrt(AstroConst::G/mstar)*
			std::pow(rstar, 5.0)*
			mplanet/Q*AstroConst::Gyr/
			std::pow(AstroConst::AU, 6.5);
	else {
		double a=semimajor*AstroConst::AU;
		return -4.5*std::sqrt(AstroConst::G/(a*mstar))*
			std::pow(rstar/a, 5.0)*
			mplanet/Q*AstroConst::Gyr/AstroConst::AU;
	}
}

double Planet::tidal_decay_semimajor_deriv(double age, double semimajor,
	double stellar_spin_frequency, bool use_a6p5) const
{
	double tidal_freq=orbital_angular_velocity_semimajor(semimajor)-
			stellar_spin_frequency,
		   Q=star->get_tidal_Q(tidal_freq),
		   dQ_da=star->get_tidal_Q_deriv(tidal_freq)*
				   orbital_angular_velocity_semimajor_deriv(semimajor),
		   da_dt=tidal_decay(age, semimajor, stellar_spin_frequency, use_a6p5);
	if(use_a6p5) return -dQ_da/Q*da_dt/6.5/std::pow(semimajor, 5.5);
	else return -(5.5/semimajor+dQ_da/Q)*da_dt;
}

double Planet::tidal_decay_star_spin_deriv(double age, double semimajor,
					   double stellar_spin_frequency, bool use_a6p5) const
{
	double tidal_freq=orbital_angular_velocity_semimajor(semimajor)-
			stellar_spin_frequency,
		   Q=star->get_tidal_Q(tidal_freq),
		   dQ_dw= star->get_tidal_Q_deriv(tidal_freq),
		   da_dt=tidal_decay(age, semimajor, stellar_spin_frequency,
				     use_a6p5);
	return -dQ_dw/Q*da_dt;
}

double Planet::tidal_torque_star_spin_deriv(double, double semimajor,
		double semi_deriv, double stellar_spin_frequency) const
{
	double tidal_freq=orbital_angular_velocity_semimajor(semimajor)-
			stellar_spin_frequency,
		   Q=star->get_tidal_Q(tidal_freq),
		   dQ_dw= star->get_tidal_Q_deriv(tidal_freq),
		   dLc_dt = orbital_angmom_deriv(semimajor, semi_deriv);
	return dQ_dw/Q*dLc_dt;
}

double Planet::tidal_decay_age_deriv(double age, double semimajor, 
				     double stellar_spin_frequency,
				     bool use_a6p5) const
{
	//double worb = std::sqrt(G)
	double Lc = star->moment_of_inertia(age, convective)*stellar_spin_frequency;
	double tidal_freq = orbital_angular_velocity_semimajor(semimajor) -
			stellar_spin_frequency;
	double Q = star->get_tidal_Q(tidal_freq);
	double Q_deriv = get_tidal_Q_age_deriv(age,semimajor,Lc);
	double da_dt =
			tidal_decay(age, semimajor, stellar_spin_frequency, use_a6p5);
	return da_dt*(5.0*star->get_logradius_deriv(age) - Q_deriv/Q);
}

double Planet::orbital_angular_momentum(double semimajor) const
{
	using namespace std;
	using namespace AstroConst;
	double a = AU*semimajor;
	return sqrt(G*a/(mplanet+mstar))*mplanet*mstar/
		(solar_mass*std::pow(solar_radius,2)/day);
}

double Planet::orbital_angmom_deriv(double semimajor, 
		double semimajor_deriv) const
{
	double a=semimajor*AstroConst::AU;
	return 0.5*mplanet*mstar*std::sqrt(AstroConst::G/
			((mplanet+mstar)*a))*
		semimajor_deriv*AstroConst::AU*AstroConst::day/
		(AstroConst::solar_mass*std::pow(AstroConst::solar_radius,2));
}

double Planet::orbit_angmom_deriv_age_deriv(double age, double semimajor,
		double semimajor_deriv, double stellar_spin_frequency) const
{
	double Lc = star->moment_of_inertia(age, convective)*stellar_spin_frequency;
	double tidal_freq = orbital_angular_velocity_semimajor(semimajor) -
			stellar_spin_frequency;
	double Q = star->get_tidal_Q(tidal_freq);
	double Q_deriv = get_tidal_Q_age_deriv(age,semimajor,Lc);
	double dLc_dt_tide = orbital_angmom_deriv(semimajor, semimajor_deriv);
	return dLc_dt_tide*(5*star->get_logradius_deriv(age) - Q_deriv/Q);
}

double Planet::orbital_angmom_deriv_semimajor_deriv(double semimajor, 
		double semimajor_deriv, double stellar_spin_frequency) const
{
	double tidal_freq=orbital_angular_velocity_semimajor(semimajor)-
			stellar_spin_frequency,
		   Q=star->get_tidal_Q(tidal_freq),
		   dQ_da=star->get_tidal_Q_deriv(tidal_freq)*
				   orbital_angular_velocity_semimajor_deriv(semimajor);
	return orbital_angmom_deriv(semimajor, semimajor_deriv)*
			(-6/semimajor - dQ_da/Q);
}

double Planet::orbital_angmom_deriv(double age) const
{
	return orbital_angmom_deriv((*semimajor_axis)(age), 
			get_semimajor_derivative(age));
}

double Planet::get_tidal_Q_age_deriv(double age, double semi, double Lc) const
{
	if(!std::isfinite(Lc)) return 0;
	double tidal_freq = orbital_angular_velocity_semimajor(semi) -
			star->spin_frequency(age, convective, Lc);
	double wc_deriv = star->spin_frequency_age_deriv(age, convective, Lc);
	double deriv = star->get_tidal_Q_deriv(tidal_freq)*(-wc_deriv);
	return deriv;
}
