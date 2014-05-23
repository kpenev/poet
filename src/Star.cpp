/**\file
 *
 * \brief The definition of some of the methods of the StarBase class.
 *
 * \ingroup StellarSystem_group
 */

#include "Star.h"
#include "Error.h"
#include <iostream>

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

StarBase::StarBase(double mass, double wind_strength,
		double wind_saturation, double coupling_timescale,
		double disk_lock_ang_vel, double disk_lock_time,
		const StellarEvolution &evolution, double age,
		double conv_spin, double rad_spin) :
	__mass(mass), __age(age),
	__disk_lock_frequency(disk_lock_ang_vel),
	__disk_dissipation_age(disk_lock_time),
	__low_mass(mass<evolution.get_mass_break()),
	__radius(evolution.interpolate_radius(mass)),
	__luminosity(evolution.interpolate_luminosity(mass)),
	__conv_angular_momentum(NULL),
	__rad_angular_momentum(NULL),
	__lifetime(9*std::pow(mass, -3)),
	__magnetic_wind_strength(wind_strength), 
	__magnetic_wind_saturation_freq(wind_saturation), 
	__core_env_coupling_timescale(coupling_timescale),
	__core_formation(evolution.core_formation_age())
{
	if(__low_mass) {
		__conv_moment_of_inertia=evolution.interpolate_moment_of_inertia(
				mass, convective);
		__rad_moment_of_inertia=evolution.interpolate_moment_of_inertia(
				mass, radiative);
		__rad_mass=evolution.interpolate_zone_mass(mass, radiative);
		__rad_radius=evolution.interpolate_core_boundary(mass);
	} else {
		__conv_moment_of_inertia=evolution.interpolate_moment_of_inertia(
				mass, total);
		__rad_moment_of_inertia=new ZeroQuantity();
		__rad_mass=new ZeroQuantity();
		__rad_radius=new ZeroQuantity();
	}
	if(!std::isnan(age)) {
		__current_conv_angular_momentum=conv_spin*
			(*__conv_moment_of_inertia)(age);
		__current_rad_angular_momentum=rad_spin*
			(*__rad_moment_of_inertia)(age);
	}
}

StarBase::~StarBase()
{
	delete __radius;
	delete __luminosity;
	delete __conv_moment_of_inertia;
	delete __rad_moment_of_inertia;
	delete __rad_mass;
	delete __rad_radius;
}

double StarBase::luminosity(double age) const
{
	if(__luminosity==NULL) throw Error::Runtime("Asking for the luminosity "
			"of a star for which the luminosity is not defined.");
	return (*__luminosity)(age);
}

double StarBase::rad_radius(double age) const
{
	if(age<__core_formation) return 0;
	return (*__rad_radius)(age);
}

double StarBase::rad_mass(double age) const
{
	if(age<__core_formation) return 0;
	return (*__rad_mass)(age);
}

double StarBase::moment_of_inertia(double age, StellarZone zone) const
{
	switch(zone) {
		case total : return (*__conv_moment_of_inertia)(age) +
					 (age<__core_formation ? 0 :
					  (*__rad_moment_of_inertia)(age));
		case radiative : return (age<__core_formation ? 0 :
								 (*__rad_moment_of_inertia)(age));
		case convective : return (*__conv_moment_of_inertia)(age);
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::moment_of_inertia."
		  );
	}
}

double StarBase::logradius_deriv(double age) const
{
	const FunctionDerivatives *deriv=__radius->deriv(age);
	double result=deriv->order(1)/deriv->order(0);
	delete deriv;
	return result;
}

double StarBase::rad_radius_deriv(double age, unsigned order) const
{
	if(age<__core_formation) return 0;
	const FunctionDerivatives *rad_radius_derivatives =
		__rad_radius->deriv(age);
	double rad_radius_deriv=rad_radius_derivatives->order(order);
	delete rad_radius_derivatives;
	return rad_radius_deriv;
}

double StarBase::rad_mass_deriv(double age) const
{
	if(age<__core_formation) return 0;
	const FunctionDerivatives *rad_mass_derivatives=__rad_mass->deriv(age);
	double rad_mass_deriv=rad_mass_derivatives->order(1);
	delete rad_mass_derivatives;
	return rad_mass_deriv;
}

double StarBase::moment_of_inertia_deriv(double age, StellarZone zone,
		int order) const
{
	double result;
	if(zone==total) {
		const FunctionDerivatives
			*dIconv=__conv_moment_of_inertia->deriv(age);
		result=dIconv->order(order);
		delete dIconv;
		if(age>__core_formation) {
			const FunctionDerivatives 
				*dIrad=__rad_moment_of_inertia->deriv(age);
			result+=dIrad->order(order);
			delete dIrad;
		}
	} else {
		if(zone==radiative && age<__core_formation) return 0;
		const FunctionDerivatives *dI=(zone==convective ?
				__conv_moment_of_inertia->deriv(age) :
				__rad_moment_of_inertia->deriv(age));
		result=dI->order(order);
		delete dI;
	}
	return result;
}

double StarBase::angular_momentum(double age, StellarZone zone) const
{
	switch(zone) {
		case total : return (*__conv_angular_momentum)(age)+
			     (age<__core_formation ? 0 : (*__rad_angular_momentum)(age));
		case radiative : return (age<__core_formation ? 0 :
								 (*__rad_angular_momentum)(age));
		case convective : return (*__conv_angular_momentum)(age);
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::get_angular_momentum."
		  );
	}
}

double StarBase::angular_momentum_deriv(double age, StellarZone zone,
		int order) {
	if(zone==radiative && age<__core_formation) return 0;
	const FunctionDerivatives* deriv;
	switch(zone) {
		case radiative: deriv = __rad_angular_momentum->deriv(age);
			break;
		case convective: deriv = __conv_angular_momentum->deriv(age);
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
					(*__rad_angular_momentum)(age));
		case convective : return spin_frequency(age, zone,
				 	(*__conv_angular_momentum)(age));
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
	if(zone==radiative && age<__core_formation) return 0;
	double I, I_deriv;
	if(zone==total) {
		const FunctionDerivatives
			*Ic_deriv=__conv_moment_of_inertia->deriv(age),
			*Ir_deriv=__rad_moment_of_inertia->deriv(age);
		I= Ic_deriv->order(0) + Ir_deriv->order(0);
		I_deriv=Ic_deriv->order(1) + Ir_deriv->order(1);
		delete Ic_deriv;
		delete Ir_deriv;
	} else {
		const FunctionDerivatives *Ideriv_func =
				(zone==convective ?
				 __conv_moment_of_inertia->deriv(age) :
				 __rad_moment_of_inertia->deriv(age));
		I=Ideriv_func->order(0);
		I_deriv=Ideriv_func->order(1);
		delete Ideriv_func;
	}
	return -(angular_momentum/std::pow(I,2)*I_deriv);
}

void StarBase::set_angular_momentum_evolution(
		const std::valarray<double> &ages, 
		const std::valarray<double> &L_values,
		const std::valarray<double> &L_derivatives,
		StellarZone zone)
{
	if(zone==convective) __conv_angular_momentum=new 
		InterpolatingFunctionALGLIB(ages, L_values, L_derivatives);
	else if(zone==radiative) __rad_angular_momentum=new 
		InterpolatingFunctionALGLIB(ages, L_values, L_derivatives);
	else throw Error::BadStellarZone("Attempting to set total angular "
			"momentum of the star.");
}

double StarBase::zone_mass(double age, StellarZone zone) const
{
	switch(zone) {
		case total : return __mass;
		case radiative : return (age<__core_formation ? 0 :
								 (*__rad_mass)(age));
		case convective :
					 return __mass-(age<__core_formation ? 0 :
							 (*__rad_mass)(age));
		default : throw Error::BadStellarZone(
		  "Unrecognized stellar zone in StarBase::get_zone_mass.");
	}
}

double StarBase::wind_torque(double age, double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	double coef=__magnetic_wind_strength*std::sqrt((*__radius)(age)/__mass);
	if(assume_wind_saturation!=SATURATED &&
			(conv_frequency<__magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED)) 
		return coef*std::pow(conv_frequency, 3);
	else {
		return coef*conv_frequency*std::pow(
			__magnetic_wind_saturation_freq, 2);
	}
}

double StarBase::wind_torque_freq_deriv(double age, double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	double coef=__magnetic_wind_strength*std::sqrt((*__radius)(age)/__mass);
	if(assume_wind_saturation!=SATURATED &&
			(conv_frequency<__magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED)) 
		return 3.0*coef*std::pow(conv_frequency, 2);
	else return coef*std::pow(__magnetic_wind_saturation_freq, 2);
}

double StarBase::wind_torque_age_deriv(double age, double angular_momentum,
		bool const_angular_momentum,
		WindSaturationState assume_wind_saturation) const
{
	using namespace std;
	const FunctionDerivatives *radius_deriv=__radius->deriv(age);
	const FunctionDerivatives  *Ic_deriv =
			__conv_moment_of_inertia->deriv(age);

	double r=radius_deriv->order(0), r_deriv=radius_deriv->order(1);
	double wc = angular_momentum/Ic_deriv->order(0);
	double wc_deriv = (const_angular_momentum ?
			spin_frequency_age_deriv(age, convective, angular_momentum) :
			0);

	delete radius_deriv;
	delete Ic_deriv;
	double coef = __magnetic_wind_strength/sqrt(__mass);
	if (assume_wind_saturation!=SATURATED &&
			(wc<__magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED))
		return coef*wc*wc*(3*sqrt(r)*wc_deriv + 0.5*wc*r_deriv/sqrt(r));
	else
		return coef*pow(__magnetic_wind_saturation_freq, 2)*
				(wc_deriv*sqrt(r) + 0.5*wc*r_deriv/sqrt(r));
}

std::complex<double> StarBase::differential_rotation_torque_angmom(
		double age, double Lconv, std::complex<double> Lrad) const
{
	return differential_rotation_torque(age, 
		differential_rotation(age, Lconv, Lrad),
		spin_frequency(age, convective, Lconv));
}

std::complex<double> StarBase::differential_rotation_torque_deriv(double age,
		double Lconv, std::complex<double> Lrad,
		StellarZone with_respect_to) const
{
	throw Error::NotImplemented("differential_rotation_torque_deriv");
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

std::complex<double> StarBase::differential_rotation_torque(double age, 
		std::complex<double> differential_rotation_amount,
		double conv_frequency) const
{
	return differential_rotation_amount/__core_env_coupling_timescale-
		core_inertia_gain(age)*conv_frequency;
}

std::complex<double> StarBase::differential_rotation_torque_deriv(double age,
		std::complex<double>, 
		std::complex<double> differential_rotation_deriv,
		double conv_frequency,
		double conv_frequency_deriv, 
		bool with_respect_to_age) const
{
	throw Error::NotImplemented("differential_rotation_torque_deriv");
	return differential_rotation_deriv/__core_env_coupling_timescale-
		(with_respect_to_age ? core_inertia_gain_deriv(age)*conv_frequency :
		 0.0) - core_inertia_gain(age)*conv_frequency_deriv;
}

std::complex<double> StarBase::differential_rotation_torque(double age) const
{
	throw Error::NotImplemented(
			"differential_rotation_torque interpolation");
	std::complex<double> diff_rot = differential_rotation(age);
	double conv_freq = spin_frequency(age, convective);
	return differential_rotation_torque(age, diff_rot, conv_freq);

}

std::complex<double> StarBase::differential_rotation(double age,
		double Lconv, std::complex<double> Lrad) 
	const
{
	double Iconv=(*__conv_moment_of_inertia)(age),
	       Irad=(age<__core_formation ? 0 : (*__rad_moment_of_inertia)(age));
	return (Iconv*Lrad-Irad*Lconv)/(Iconv+Irad);
}

std::complex<double> StarBase::differential_rotation_deriv(double age,
		double Lconv, std::complex<double> Lrad, StellarZone with_respect_to)
	const
{
	throw Error::NotImplemented("differential_rotation_deriv");
	using namespace std;
	if(with_respect_to==total) {
		if(age<__core_formation) return 0;
		const FunctionDerivatives 
			*Iconv_deriv=__conv_moment_of_inertia->deriv(age),
			*Irad_deriv=__rad_moment_of_inertia->deriv(age);
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
		double Iconv=(*__conv_moment_of_inertia)(age),
			   Irad=(age<__core_formation ? 0 :
					   (*__rad_moment_of_inertia)(age));
		if(with_respect_to==convective) return -Irad/(Iconv+Irad);
		else return Iconv/(Iconv+Irad);
	}
}

std::complex<double> StarBase::differential_rotation(double age) const
{
	throw Error::NotImplemented("Differential rotation interpolation");
	return differential_rotation(age, (*__conv_angular_momentum)(age), 
			(*__rad_angular_momentum)(age));
}

double StarBase::core_inertia_gain(double age) const
{
	if(age<__core_formation) return 0.0;
	const FunctionDerivatives *rad_mass_derivatives=
		__rad_mass->deriv(age);
	double rad_mass_deriv=rad_mass_derivatives->order(1);
	delete rad_mass_derivatives;
	return 2.0/3.0*std::pow((*__rad_radius)(age), 2)*rad_mass_deriv;
}

double StarBase::core_inertia_gain_deriv(double age) const
{
	if(age<__core_formation) return 0;
	const FunctionDerivatives *Rrad_deriv=__rad_radius->deriv(age),
			    *Mrad_deriv=__rad_mass->deriv(age);
	double r_rad=Rrad_deriv->order(0), r_rad_deriv=Rrad_deriv->order(1),
	       m_rad_deriv=Mrad_deriv->order(1),
	       m_rad_deriv2=Mrad_deriv->order(2);
	delete Rrad_deriv;
	delete Mrad_deriv;

	return 2.0/3.0*(2.0*r_rad*m_rad_deriv*r_rad_deriv+
			r_rad*r_rad*m_rad_deriv2);
}
