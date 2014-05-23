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

void StarBase::delete_current_age_quantities()
{
	for(unsigned i=0; i<__current_age_quantities.size(); ++i)
		if(__current_age_quantities[i]) {
			std::cout << "Deleting quantity #" << i << std::endl;
			delete __current_age_quantities[i];
			__current_age_quantities[i]=NULL;
		}
}

double StarBase::current_age_quantity(CurrentAgeQuantities quantity,
		unsigned deriv_order) const
{
#ifdef DEBUG
	if(quantity==CORE_INERTIA_GAIN) assert(deriv_order<=1);
#endif

	if(__current_age_quantities[quantity]==NULL) {
		if(quantity==MCONV || quantity==MRAD) 
			std::cout << "Creating quantity #" << MCONV << std::endl
				<< "Creating quantity #" << MRAD << std::endl;
		else if(quantity==ICONV || quantity==IRAD || quantity==ITOT)
			std::cout << "Creating quantity #" << ICONV << std::endl
				<< "Creating quantity #" << IRAD << std::endl
				<< "Creating quantity #" << ITOT << std::endl;
		else std::cout << "Creating quantity #" << quantity << std::endl;
		switch(quantity) {
			case RADIUS :
				__current_age_quantities[RADIUS]=__radius->deriv(__age);
				break;
			case LUMINOSITY :
				__current_age_quantities[LUMINOSITY]=
					__luminosity->deriv(__age);
				break;
			case RRAD :
				__current_age_quantities[RRAD]= (__age<__core_formation ?
						new CubicSplineDerivatives(0, 0, 0) : 
						__rad_radius->deriv(__age));
				break;
			case MCONV : case MRAD :
				__current_age_quantities[MRAD]=(
					__age<__core_formation ?
					new CubicSplineDerivatives(0, 0, 0) :
					__rad_mass->deriv(__age));
				__current_age_quantities[MCONV]=new CubicSplineDerivatives(
						__mass-__current_age_quantities[MRAD]->order(0),
						__current_age_quantities[MRAD]->order(1),
						__current_age_quantities[MRAD]->order(2));
				break;
			case ICONV : case IRAD : case ITOT :
				__current_age_quantities[ICONV]=
					__conv_moment_of_inertia->deriv(__age);
				__current_age_quantities[IRAD]=(__age<__core_formation ?
					new CubicSplineDerivatives(0, 0, 0) :
					__rad_moment_of_inertia->deriv(__age));
				__current_age_quantities[ITOT]=new CubicSplineDerivatives(
						__current_age_quantities[ICONV]->order(0)+
						__current_age_quantities[IRAD]->order(0),

						__current_age_quantities[ICONV]->order(1)+
						__current_age_quantities[IRAD]->order(1),

						__current_age_quantities[ICONV]->order(2)+
						__current_age_quantities[IRAD]->order(2));
				break;
			case CORE_INERTIA_GAIN :
				__current_age_quantities[CORE_INERTIA_GAIN]=
					new CubicSplineDerivatives(core_inertia_gain(__age), 
							core_inertia_gain_deriv(__age), NaN);
				break;
			default :
#ifdef DEBUG
				assert(false)
#endif
					;
		}
	}
	return __current_age_quantities[quantity]->order(deriv_order);
}

StarBase::StarBase(double mass, double wind_strength,
		double wind_saturation, double coupling_timescale,
		double disk_lock_ang_vel, double disk_lock_time,
		const StellarEvolution &evolution, double age,
		double conv_spin, double rad_spin) :
	__mass(mass), __age(NaN),
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
	__core_formation(evolution.core_formation_age()),
	__current_age_quantities(NUM_CURRENT_AGE_QUANTITIES, NULL)
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
	delete_current_age_quantities();
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

double StarBase::moment_of_inertia(StellarZone zone) const
{
	switch(zone) {
		case total : return current_age_quantity(ITOT);
		case convective : return current_age_quantity(ICONV);
		case radiative : return current_age_quantity(IRAD);
		default : throw Error::BadStellarZone("in current age "
						  "moment_of_inertia");
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

double StarBase::moment_of_inertia_deriv(StellarZone zone,
		int order) const
{
	switch(zone) {
		case total : return current_age_quantity(ITOT, order);
		case convective : return current_age_quantity(ICONV, order);
		case radiative : return current_age_quantity(IRAD, order);
		default : throw Error::BadStellarZone("in current age "
						  "moment_of_inertia_deriv");
	}
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

double StarBase::zone_mass(StellarZone zone) const
{
	switch(zone) {
		case total : return __mass;
		case convective : return current_age_quantity(MCONV);
		case radiative : return current_age_quantity(MRAD);
		default : throw Error::BadStellarZone("in current age zone_mass");
	}
}

double StarBase::wind_torque(double age, double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	double coef=__magnetic_wind_strength*std::sqrt(
			(age==__age ? current_age_quantity(RADIUS) : (*__radius)(age))/
			__mass);
	if(assume_wind_saturation!=SATURATED &&
			(conv_frequency<__magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED)) 
		return coef*std::pow(conv_frequency, 3);
	else {
		return coef*conv_frequency*std::pow(
			__magnetic_wind_saturation_freq, 2);
	}
}

double StarBase::wind_torque(double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	return wind_torque(__age, conv_frequency, assume_wind_saturation);
}

double StarBase::wind_torque_freq_deriv(double age, double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	double coef=__magnetic_wind_strength*std::sqrt(
			(age==__age ? current_age_quantity(RADIUS) : (*__radius)(age))/
			__mass);
	if(assume_wind_saturation!=SATURATED &&
			(conv_frequency<__magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED)) 
		return 3.0*coef*std::pow(conv_frequency, 2);
	else return coef*std::pow(__magnetic_wind_saturation_freq, 2);
}

double StarBase::wind_torque_freq_deriv(double conv_frequency,
		WindSaturationState assume_wind_saturation) const
{
	return wind_torque_freq_deriv(__age, conv_frequency,
			assume_wind_saturation);
}

double StarBase::wind_torque_age_deriv(double age, double angular_momentum,
		bool const_angular_momentum,
		WindSaturationState assume_wind_saturation) const
{
	using namespace std;
	double r, r_deriv, wc, wc_deriv;
	if(age==__age) {
		r=current_age_quantity(RADIUS, 0);
		r_deriv=current_age_quantity(RADIUS, 1);
		wc = angular_momentum/current_age_quantity(ICONV, 0);
		wc_deriv = (const_angular_momentum ?
				spin_frequency_age_deriv(convective, angular_momentum) : 0);

	} else {
		const FunctionDerivatives *radius_deriv=__radius->deriv(age);
		const FunctionDerivatives  *Ic_deriv =
			__conv_moment_of_inertia->deriv(age);

		r=radius_deriv->order(0);
		r_deriv=radius_deriv->order(1);
		wc = angular_momentum/Ic_deriv->order(0);
		wc_deriv = (const_angular_momentum ?
				spin_frequency_age_deriv(age, convective,
					angular_momentum) : 0);
		delete radius_deriv;
		delete Ic_deriv;
	}
	double coef = __magnetic_wind_strength/sqrt(__mass);
	if (assume_wind_saturation!=SATURATED &&
			(wc<__magnetic_wind_saturation_freq ||
			 assume_wind_saturation==NOT_SATURATED))
		return coef*wc*wc*(3*sqrt(r)*wc_deriv + 0.5*wc*r_deriv/sqrt(r));
	else
		return coef*pow(__magnetic_wind_saturation_freq, 2)*
				(wc_deriv*sqrt(r) + 0.5*wc*r_deriv/sqrt(r));
}

double StarBase::wind_torque_age_deriv(double angular_momentum,
		bool const_angular_momentum,
		WindSaturationState assume_wind_saturation) const
{
	return wind_torque_age_deriv(__age, angular_momentum,
			const_angular_momentum, assume_wind_saturation);
}

std::complex<double> StarBase::differential_rotation_torque_angmom(
		double age, double Lconv, std::complex<double> Lrad) const
{
	return differential_rotation_torque(age, 
		differential_rotation(age, Lconv, Lrad),
		spin_frequency(age, convective, Lconv));
}

std::complex<double> StarBase::differential_rotation_torque_angmom(
		double Lconv, std::complex<double> Lrad) const
{
	return differential_rotation_torque(differential_rotation(Lconv, Lrad),
			spin_frequency(convective, Lconv));
}

std::complex<double> StarBase::differential_rotation_torque_deriv(double age,
		double Lconv, std::complex<double> Lrad,
		StellarZone with_respect_to) const
{
	throw Error::NotImplemented("differential_rotation_torque_deriv");
	double w_deriv;
	if (with_respect_to == total) {
		w_deriv = (age==__age ?
				spin_frequency_age_deriv(convective, Lconv) :
				spin_frequency_age_deriv(age, convective, Lconv));
	}
	else if (with_respect_to == convective)
		w_deriv = (age==__age ?
				spin_frequency_angmom_deriv(convective, Lconv) :
				spin_frequency_angmom_deriv(age, convective, Lconv));
	else w_deriv = 0.0;
	bool with_respect_to_age = (with_respect_to == total);

	if(age==__age)
		return differential_rotation_torque_deriv(age, 
				differential_rotation(age, Lconv, Lrad),
				differential_rotation_deriv(age, Lconv, Lrad, 
					with_respect_to),
				spin_frequency(age, convective, Lconv), w_deriv,
				with_respect_to_age);
	else 
		return differential_rotation_torque_deriv(
				differential_rotation(age, Lconv, Lrad),
				differential_rotation_deriv(age, Lconv, Lrad, 
					with_respect_to),
				spin_frequency(age, convective, Lconv), w_deriv,
				with_respect_to_age);
}

std::complex<double> StarBase::differential_rotation_torque_deriv(
		double Lconv, std::complex<double> Lrad,
		StellarZone with_respect_to) const
{
	return differential_rotation_torque_deriv(__age, Lconv, Lrad,
			with_respect_to);
}

std::complex<double> StarBase::differential_rotation_torque(double age, 
		std::complex<double> differential_rotation_amount,
		double conv_frequency) const
{
	return differential_rotation_amount/__core_env_coupling_timescale-
		(age==__age ? core_inertia_gain() : core_inertia_gain(age))*
		conv_frequency;
}

std::complex<double> StarBase::differential_rotation_torque(
		std::complex<double> differential_rotation_amount,
		double conv_frequency) const
{
	return differential_rotation_torque(__age, differential_rotation_amount,
		conv_frequency);
}

std::complex<double> StarBase::differential_rotation_torque_deriv(double age,
		std::complex<double>, 
		std::complex<double> differential_rotation_deriv,
		double conv_frequency,
		double conv_frequency_deriv, 
		bool with_respect_to_age) const
{
	throw Error::NotImplemented("differential_rotation_torque_deriv");
	double dIcore=(age==__age ? core_inertia_gain():core_inertia_gain(age)),
		   d2Icore=(age==__age ? core_inertia_gain_deriv() :
				   core_inertia_gain_deriv(age));
	return differential_rotation_deriv/__core_env_coupling_timescale -
		(with_respect_to_age ? d2Icore*conv_frequency : 0.0) - 
		dIcore*conv_frequency_deriv;
}

std::complex<double> StarBase::differential_rotation_torque_deriv(
		std::complex<double> differential_rotation_amount, 
		std::complex<double> differential_rotation_deriv,
		double conv_frequency,
		double conv_frequency_deriv, 
		bool with_respect_to_age) const
{
	return differential_rotation_torque_deriv(__age,
			differential_rotation_amount, differential_rotation_deriv,
			conv_frequency, conv_frequency_deriv, with_respect_to_age);
}

std::complex<double> StarBase::differential_rotation_torque(double age) const
{
	throw Error::NotImplemented(
			"differential_rotation_torque interpolation");
	std::complex<double> diff_rot = (age==__age ? differential_rotation() :
			differential_rotation(age));
	double conv_freq = (age==__age ? spin_frequency(convective) :
			spin_frequency(age, convective));
	if(age==__age)
		return differential_rotation_torque(diff_rot, conv_freq);
	else 
		return differential_rotation_torque(age, diff_rot, conv_freq);

}

std::complex<double> StarBase::differential_rotation(double age,
		double Lconv, std::complex<double> Lrad) 
	const
{
	double Iconv, Irad;
	if(age==__age) {
		Iconv=current_age_quantity(ICONV);
		Irad=current_age_quantity(IRAD);
	} else {
		Iconv=(*__conv_moment_of_inertia)(age);
		Irad=(age<__core_formation ? 0 : (*__rad_moment_of_inertia)(age));
	}
	return (Iconv*Lrad-Irad*Lconv)/(Iconv+Irad);
}

std::complex<double> StarBase::differential_rotation(double Lconv,
		std::complex<double> Lrad) const
{
	return differential_rotation(__age, Lconv, Lrad);
}

std::complex<double> StarBase::differential_rotation_deriv(double age,
		double Lconv, std::complex<double> Lrad, StellarZone with_respect_to)
	const
{
	throw Error::NotImplemented("differential_rotation_deriv");
	using namespace std;
	if(with_respect_to==total) {
		if(age<__core_formation) return 0;
		double Ic, Ic_deriv, Ir, Ir_deriv;
		if(age==__age) {
			Ic=current_age_quantity(ICONV);
			Ic_deriv=current_age_quantity(ICONV, 1);
			Ir=current_age_quantity(IRAD);
			Ir_deriv=current_age_quantity(IRAD, 1);
		} else {
			const FunctionDerivatives 
				*Iconv_deriv=__conv_moment_of_inertia->deriv(age),
				*Irad_deriv=__rad_moment_of_inertia->deriv(age);
			Ic=Iconv_deriv->order(0);
			Ic_deriv=Iconv_deriv->order(1);
			Ir = Irad_deriv->order(0);
			Ir_deriv = Irad_deriv->order(1);
			delete Iconv_deriv;
			delete Irad_deriv;
		}
		return (Ic_deriv*Lrad-Ir_deriv*Lconv)/(Ic+Ir)-
			 (Ic*Lrad-Ir*Lconv)/std::pow(Ic+Ir, 2)*
			 (Ic_deriv+Ir_deriv);
	} else {
		double Iconv=(age==__age ? current_age_quantity(ICONV) :
				(*__conv_moment_of_inertia)(age)),
			   Irad=(age<__core_formation ? 0 :
					   (age==__age ? current_age_quantity(IRAD) :
						(*__rad_moment_of_inertia)(age)));
		if(with_respect_to==convective) return -Irad/(Iconv+Irad);
		else return Iconv/(Iconv+Irad);
	}
}

std::complex<double> StarBase::differential_rotation_deriv(double Lconv,
		std::complex<double> Lrad, StellarZone with_respect_to) const
{
	return differential_rotation_deriv(__age, Lconv, Lrad,with_respect_to);
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
