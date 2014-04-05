/**\file
 *
 * \brief The definition of some of the methods of the StarBase class.
 *
 * \ingroup StellarSystem_group
 */

#include "Star.h"
#include "Error.h"
#include <iostream>
#include <assert.h>

std::ostream &operator<<(std::ostream &os,
		const WindSaturationState &wind_state)
{
	switch(wind_state) {
		case NOT_SATURATED: os << "NOT_SATURATED"; break;
		case UNKNOWN: os << "UNKNOWN"; break;
		case SATURATED: os << "SATURATED"; break;
	}
	return os;
}

StarBase::StarBase(double current_mass, double wind_strength,
		double wind_saturation, double coupling_timescale,
		double disk_lock_ang_vel, double disk_lock_time,
		const StellarEvolution &evolution, double current_age,
		double current_conv_spin, double current_rad_spin) :
	mass(current_mass), age(current_age),
	disk_lock_frequency(disk_lock_ang_vel),
	disk_dissipation_age(disk_lock_time),
	conv_angular_momentum(NULL), 
	rad_angular_momentum(NULL),
	magnetic_wind_strength(wind_strength), 
	magnetic_wind_saturation_freq(wind_saturation), 
	core_env_coupling_timescale(coupling_timescale),
	core_formation(evolution.core_formation_age())
{
	lifetime = 9*std::pow(mass, -3);
	low_mass=(mass<evolution.get_mass_break());
	radius=evolution.interpolate_radius(mass);
	luminosity=evolution.interpolate_luminosity(mass);
	if(low_mass) {
		conv_moment_of_inertia=evolution.interpolate_moment_of_inertia(
				mass, convective);
		rad_moment_of_inertia=evolution.interpolate_moment_of_inertia(
				mass, radiative);
		rad_mass=evolution.interpolate_zone_mass(mass, radiative);
		rad_radius=evolution.interpolate_core_boundary(mass);
	} else {
		conv_moment_of_inertia=evolution.interpolate_moment_of_inertia(
				mass, total);
		rad_moment_of_inertia=new ZeroQuantity();
		rad_mass=new ZeroQuantity();
		rad_radius=new ZeroQuantity();
	}
	if(!std::isnan(current_age)) {
		current_conv_angular_momentum=current_conv_spin*
			(*conv_moment_of_inertia)(current_age);
		current_rad_angular_momentum=current_rad_spin*
			(*rad_moment_of_inertia)(current_age);
	}
}

StarBase::~StarBase()
{
	delete radius;
	delete luminosity;
	delete conv_moment_of_inertia;
	delete rad_moment_of_inertia;
	delete rad_mass;
	delete rad_radius;
}

/*
//Same as above, but without unnecessary arguments
StarBase::StarBase(double current_mass, double tidal_quality,
		double wind_strength, double wind_saturation,
		double coupling_timescale, double dissipation_transition_width,
		const StellarEvolution &evolution) :
	mass(current_mass),
	Q_transition_width(dissipation_transition_width),
	conv_angular_momentum(NULL),
	rad_angular_momentum(NULL),
	tidal_Q(tidal_quality),
	magnetic_wind_strength(wind_strength),
	magnetic_wind_saturation_freq(wind_saturation),
	core_env_coupling_timescale(coupling_timescale),
	current_conv_angular_momentum(0),
	current_rad_angular_momentum(0),
	core_formation(evolution.core_formation_age())
{
	const double FAKE_AGE = 5;
	radius=evolution.interpolate_radius(mass, FAKE_AGE);
	conv_moment_of_inertia=evolution.interpolate_moment_of_inertia(
			mass, convective, FAKE_AGE);
	tot_moment_of_inertia=evolution.interpolate_moment_of_inertia(
			mass, total, FAKE_AGE);
	rad_mass=evolution.interpolate_zone_mass(mass, radiative);
	rad_radius=evolution.interpolate_core_boundary(mass, FAKE_AGE);
	lifetime = 9*std::pow(mass, -3);
}*/

double StarBase::current_age() const
{
	return age;
}

double StarBase::get_mass() const
{
	return mass;
}

double StarBase::get_radius(double age) const
{
	return (*radius)(age);
}

double StarBase::get_luminosity(double age) const
{
	if(luminosity==NULL) throw Error::Runtime("Asking for the luminosity of "
			"a star for which the luminosity is not defined.");
	return (*luminosity)(age);
}

double StarBase::core_formation_age() const
{
	return core_formation;
}

double StarBase::get_disk_lock_frequency() const
{
	return disk_lock_frequency;
}

void StarBase::set_disk_lock_frequency(double value)
{
	disk_lock_frequency=value;
}

double StarBase::get_disk_dissipation_age() const
{
	return disk_dissipation_age;
}

double StarBase::get_rad_radius(double age) const
{
	if(age<core_formation) return 0;
	return (*rad_radius)(age);
}

double StarBase::get_rad_radius_deriv(double age, unsigned order) const
{
	if(age<core_formation) return 0;
	const FunctionDerivatives *rad_radius_derivatives =
		rad_radius->deriv(age);
	double rad_radius_deriv=rad_radius_derivatives->order(order);
	delete rad_radius_derivatives;
	return rad_radius_deriv;
}

double StarBase::get_rad_mass(double age) const
{
	if(age<core_formation) return 0;
	return (*rad_mass)(age);
}

double StarBase::get_rad_mass_deriv(double age) const
{
	if(age<core_formation) return 0;
	const FunctionDerivatives *rad_mass_derivatives=
		rad_mass->deriv(age);
	double rad_mass_deriv=rad_mass_derivatives->order(1);
	delete rad_mass_derivatives;
	return rad_mass_deriv;
}

double StarBase::get_logradius_deriv(double age) const
{
	const FunctionDerivatives *deriv=radius->deriv(age);
	double result=deriv->order(1)/deriv->order(0);
	delete deriv;
	return result;
}

///Returns the present radius of the star (in solar radii).
/*double StarBase::get_radius() const
{
	return (*radius)(curr_age);
}*/

double StarBase::moment_of_inertia(double age, StellarZone zone) const
{
	switch(zone) {
		case total : return (*conv_moment_of_inertia)(age) +
					 (age<core_formation ? 0 :
					  (*rad_moment_of_inertia)(age));
		case radiative : return (age<core_formation ? 0 :
								 (*rad_moment_of_inertia)(age));
		case convective : return (*conv_moment_of_inertia)(age);
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::moment_of_inertia."
		  );
	}
}

double StarBase::moment_of_inertia_deriv(double age, StellarZone zone,
		int order) const
{
	double result;
	if(zone==total) {
		const FunctionDerivatives *dIconv=conv_moment_of_inertia->deriv(age);
		result=dIconv->order(order);
		delete dIconv;
		if(age>core_formation) {
			const FunctionDerivatives 
				*dIrad=rad_moment_of_inertia->deriv(age);
			result+=dIrad->order(order);
			delete dIrad;
		}
	} else {
		if(zone==radiative && age<core_formation) return 0;
		const FunctionDerivatives *dI=(zone==convective ?
				conv_moment_of_inertia->deriv(age) :
				rad_moment_of_inertia->deriv(age));
		result=dI->order(order);
		delete dI;
	}
	return result;
}

double StarBase::get_angular_momentum(double age, StellarZone zone) const
{
	switch(zone) {
		case total : return (*conv_angular_momentum)(age)+
			     (age<core_formation ? 0 : (*rad_angular_momentum)(age));
		case radiative : return (age<core_formation ? 0 :
								 (*rad_angular_momentum)(age));
		case convective : return (*conv_angular_momentum)(age);
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::get_angular_momentum."
		  );
	}
}

double StarBase::angular_momentum_deriv(double age, StellarZone zone, int order) {
	if(zone==radiative && age<core_formation) return 0;
	const FunctionDerivatives* deriv;
	switch(zone) {
		case radiative: deriv = rad_angular_momentum->deriv(age);
			break;
		case convective: deriv = conv_angular_momentum->deriv(age);
			break;
		default: throw Error::BadStellarZone(
				"Unrecognized zone in StarBase::angular_momentum");
	}
	double result = deriv->order(order);
	delete deriv;
	return result;
}

double StarBase::spin_frequency(double age, StellarZone zone) const
{
	switch(zone) {
		case total : throw Error::BadStellarZone(
		     "Spin of the entire star is not defined");
		case radiative : return spin_frequency(age, zone,
					(*rad_angular_momentum)(age));
		case convective : return spin_frequency(age, zone,
				 	(*conv_angular_momentum)(age));
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::spin_frequency.");
	}
}

double StarBase::spin_frequency(double age, StellarZone zone,
		double angular_momentum) const
{
	return angular_momentum/moment_of_inertia(age, zone);
}

double StarBase::spin_frequency_age_deriv(double age, StellarZone zone,
		double angular_momentum) const
{
	if(zone==radiative && age<core_formation) return 0;
	double I, I_deriv;
	if(zone==total) {
		const FunctionDerivatives
			*Ic_deriv=conv_moment_of_inertia->deriv(age),
			*Ir_deriv=rad_moment_of_inertia->deriv(age);
		I= Ic_deriv->order(0) + Ir_deriv->order(0);
		I_deriv=Ic_deriv->order(1) + Ir_deriv->order(1);
		delete Ic_deriv;
		delete Ir_deriv;
	} else {
		const FunctionDerivatives *Ideriv_func =
				(zone==convective ?
				 conv_moment_of_inertia->deriv(age) :
				 rad_moment_of_inertia->deriv(age));
		I=Ideriv_func->order(0);
		I_deriv=Ideriv_func->order(1);
		delete Ideriv_func;
	}
	return -(angular_momentum/std::pow(I,2)*I_deriv);
}

double StarBase::spin_frequency_angmom_deriv(double age, StellarZone zone,
		double) const
{
	return 1.0/moment_of_inertia(age, zone);
}

double StarBase::spin_period(double age, StellarZone zone) const
{
	return 2.0*M_PI/spin_frequency(age, zone);
}

void StarBase::set_angular_momentum_evolution(
		const std::valarray<double> &ages, 
		const std::valarray<double> &L_values,
		const std::valarray<double> &L_derivatives,
		StellarZone zone)
{
	if(zone==convective) conv_angular_momentum=new 
		InterpolatingFunctionALGLIB(ages, L_values, L_derivatives);
	else if(zone==radiative) rad_angular_momentum=new 
		InterpolatingFunctionALGLIB(ages, L_values, L_derivatives);
	else throw Error::BadStellarZone("Attempting to set total angular "
			"momentum of the star.");
}

double StarBase::get_zone_mass(double age, StellarZone zone) const
{
	switch(zone) {
		case total : return mass;
		case radiative : return (age<core_formation ? 0 : (*rad_mass)(age));
		case convective :
					 return mass-(age<core_formation ? 0 : (*rad_mass)(age));
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::get_zone_mass.");
	}
}

double StarBase::wind_torque(double age, double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	double coef=magnetic_wind_strength*std::sqrt((*radius)(age)/mass);
	if(assume_wind_saturation!=SATURATED &&
			(conv_frequency<magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED)) 
		return coef*std::pow(conv_frequency, 3);
	else {
		return coef*conv_frequency*std::pow(
			magnetic_wind_saturation_freq, 2);
	}
}

double StarBase::wind_torque_freq_deriv(double age, double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	double coef=magnetic_wind_strength*std::sqrt((*radius)(age)/mass);
	if(assume_wind_saturation!=SATURATED &&
			(conv_frequency<magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED)) 
		return 3.0*coef*std::pow(conv_frequency, 2);
	else return coef*std::pow(magnetic_wind_saturation_freq, 2);
}

double StarBase::wind_torque_age_deriv(double age, double angular_momentum,
		bool const_angular_momentum,
		WindSaturationState assume_wind_saturation) const
{
	using namespace std;
	const FunctionDerivatives *radius_deriv=radius->deriv(age);
	const FunctionDerivatives  *Ic_deriv =
			conv_moment_of_inertia->deriv(age);

	double r=radius_deriv->order(0), r_deriv=radius_deriv->order(1);
	double wc = angular_momentum/Ic_deriv->order(0);
	double wc_deriv = (const_angular_momentum ?
			spin_frequency_age_deriv(age, convective, angular_momentum) :
			0);

	delete radius_deriv;
	delete Ic_deriv;
	double coef = magnetic_wind_strength/sqrt(mass);
	if (assume_wind_saturation!=SATURATED &&
			(wc<magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED))
		return coef*wc*wc*(3*sqrt(r)*wc_deriv + 0.5*wc*r_deriv/sqrt(r));
	else
		return coef*pow(magnetic_wind_saturation_freq, 2)*
				(wc_deriv*sqrt(r) + 0.5*wc*r_deriv/sqrt(r));

}

double StarBase::wind_torque(double age) const
{
	return wind_torque(age, spin_frequency(age, envelope));
}

double StarBase::differential_rotation_torque_angmom(double age, 
		double Lconv,
		double Lrad) const
{
	return differential_rotation_torque(age, 
		differential_rotation(age, Lconv, Lrad),
		spin_frequency(age, convective, Lconv));
}

double StarBase::differential_rotation_torque_deriv(double age, double Lconv,
		double Lrad, StellarZone with_respect_to) const
{
	double w_deriv;
	if (with_respect_to == total) {
		w_deriv = spin_frequency_age_deriv(age, convective, Lconv);
	}
	else if (with_respect_to == convective)
		w_deriv = spin_frequency_angmom_deriv(age, convective, Lconv);
	else w_deriv = 0.0;
	bool with_respect_to_age = (with_respect_to == total);

	return differential_rotation_torque_deriv(age, 
			differential_rotation(age, Lconv, Lrad),
			differential_rotation_deriv(age, Lconv, Lrad, 
				with_respect_to),
			spin_frequency(age, convective, Lconv), w_deriv,
			with_respect_to_age);
}

double StarBase::differential_rotation_torque(double age, 
		double differential_rotation_amount,
		double conv_frequency) const
{
	return differential_rotation_amount/core_env_coupling_timescale-
		core_inertia_gain(age)*conv_frequency;
}

double StarBase::differential_rotation_torque_deriv(double age, 
		double, 
		double differential_rotation_deriv,
		double conv_frequency,
		double conv_frequency_deriv, 
		bool with_respect_to_age) const
{
	return differential_rotation_deriv/core_env_coupling_timescale-
		(with_respect_to_age ? core_inertia_gain_deriv(age)*conv_frequency :
		 0.0) - core_inertia_gain(age)*conv_frequency_deriv;
}

double StarBase::differential_rotation_torque(double age) const
{
	double diff_rot = differential_rotation(age);
	double conv_freq = spin_frequency(age, convective);
	return differential_rotation_torque(age, diff_rot, conv_freq);

}

double StarBase::differential_rotation(double age, double Lconv, double Lrad) 
	const
{
	double Iconv=(*conv_moment_of_inertia)(age),
	       Irad=(age<core_formation ? 0 : (*rad_moment_of_inertia)(age));
	return (Iconv*Lrad-Irad*Lconv)/(Iconv+Irad);
}

double StarBase::differential_rotation_deriv(double age, double Lconv, 
	double Lrad, StellarZone with_respect_to) const
{
	using namespace std;
	if(with_respect_to==total) {
		if(age<core_formation) return 0;
		const FunctionDerivatives 
			*Iconv_deriv=conv_moment_of_inertia->deriv(age),
			*Irad_deriv=rad_moment_of_inertia->deriv(age);
		double Ic=Iconv_deriv->order(0),
		       Ic_deriv=Iconv_deriv->order(1);

		double Ir = Irad_deriv->order(0),
				Ir_deriv = Irad_deriv->order(1);
		delete Iconv_deriv;
		delete Irad_deriv;
		return (Ic_deriv*Lrad-Ir_deriv*Lconv)/(Ic+Ir)-
			 (Ic*Lrad-Ir*Lconv)/std::pow(Ic+Ir, 2)*
			 (Ic_deriv+Ir_deriv);
	} else {
		double Iconv=(*conv_moment_of_inertia)(age),
			   Irad=(age<core_formation ? 0 : (*rad_moment_of_inertia)(age));
		if(with_respect_to==convective) return -Irad/(Iconv+Irad);
		else return Iconv/(Iconv+Irad);
	}
}

double StarBase::differential_rotation(double age) const
{
	return differential_rotation(age, (*conv_angular_momentum)(age), 
			(*rad_angular_momentum)(age));
}

double StarBase::core_inertia_gain(double age) const
{
	if(age<core_formation) return 0.0;
	const FunctionDerivatives *rad_mass_derivatives=
		rad_mass->deriv(age);
	double rad_mass_deriv=rad_mass_derivatives->order(1);
	delete rad_mass_derivatives;
	return 2.0/3.0*std::pow((*rad_radius)(age), 2)*rad_mass_deriv;
}

double StarBase::core_inertia_gain_deriv(double age) const
{
	if(age<core_formation) return 0;
	const FunctionDerivatives *Rrad_deriv=rad_radius->deriv(age),
			    *Mrad_deriv=rad_mass->deriv(age);
	double r_rad=Rrad_deriv->order(0), r_rad_deriv=Rrad_deriv->order(1),
	       m_rad_deriv=Mrad_deriv->order(1),
	       m_rad_deriv2=Mrad_deriv->order(2);
	delete Rrad_deriv;
	delete Mrad_deriv;

	return 2.0/3.0*(2.0*r_rad*m_rad_deriv*r_rad_deriv+
			r_rad*r_rad*m_rad_deriv2);
}

double StarBase::get_lifetime() const
{
	return lifetime;
}

double StarBase::get_current_angular_momentum(StellarZone zone) const
{
	if(zone==convective) return current_conv_angular_momentum;
	else if(zone==radiative) return current_rad_angular_momentum;
	else return current_conv_angular_momentum+
		current_rad_angular_momentum;
}
