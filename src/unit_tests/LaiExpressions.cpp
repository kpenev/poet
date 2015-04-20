#include "LaiExpressions.h"

double dimensionless_torque_x_Lai(double inclination, const Lags &lags)
{
	double c=std::cos(inclination), s=std::sin(inclination);
	return 0.15*M_PI*(s*std::pow(1+c, 3)*lags(2,2)/2
			+ s*std::pow(1+c, 2)*(2.0-c)*lags(1,2)
			+ 3.0*std::pow(s, 3)*lags(0,2)
			+ s*std::pow(1-c, 2)*(2.0+c)*lags(-1, 2)
			+ s*std::pow(1-c, 3)*lags(-2, 2)/2)
		- 3.0*M_PI/5.0*(std::pow(s, 3)*c*lags(2,0)/2
				+ s*std::pow(c, 3)*lags(1,0));
}

double dimensionless_torque_z_Lai(double inclination, const Lags &lags)
{
	double c=std::cos(inclination), s=std::sin(inclination);
	return 0.15*M_PI*(
			std::pow(1+c, 4)/2*lags(2, 2)
			+ std::pow(s*(1+c), 2)*lags(1, 2)
			- std::pow(s*(1-c), 2)*lags(-1, 2)
			- std::pow(1-c, 4)/2*lags(-2, 2))
		+ 3.0*M_PI/5.0*(std::pow(s, 4)/2*lags(2, 0)
				+ std::pow(s*c, 2)*lags(1, 0));
}

double dimensionless_power_Lai(double inclination, const Lags &lags)
{
	double c=std::cos(inclination), s=std::sin(inclination);
	return 0.15*M_PI*(
			std::pow(1+c, 4)*lags(2,2)/2
			+ 2.0*std::pow(s*(1+c), 2)*lags(1,2)
			+ 3.0*std::pow(s, 4)*lags(0, 2)
			+ 2.0*std::pow(s*(1-c), 2)*lags(-1, 2)
			+ std::pow(1-c, 4)*lags(-2, 2)/2);
}

double torque_norm_Lai(double perturber_mass, double dissipator_radius,
		double semimajor)
{
	return AstroConst::G*std::pow(perturber_mass/std::pow(semimajor, 3), 2)
		   *std::pow(dissipator_radius, 5)*AstroConst::solar_mass
		   /std::pow(AstroConst::solar_radius, 3)
		   *AstroConst::day*AstroConst::Gyr;
}

double power_norm_Lai(double dissipator_mass, double perturber_mass,
					  double dissipator_radius, double semimajor)
{
	return torque_norm_Lai(perturber_mass, dissipator_radius, semimajor)
		   *orbital_angular_velocity(dissipator_mass, perturber_mass,
				   					 semimajor);
}
