/**\file
 *
 * \brief Declarations of some methods of the TidalDissipation class.
 *
 * \ingroup StellarSystem_group
 */

#include "TidalDissipation.h"

const double TidalDissipation::__Umm_coef[][3]={
	{std::sqrt(3.0*M_PI/10.0)/4.0,
	 -std::sqrt(6.0*M_PI/5.0)/4.0,
	 std::sqrt(3.0*M_PI/10.0)/4.0},

	{-std::sqrt(3.0*M_PI/10.0)/2.0,
	 std::sqrt(6.0*M_PI/5.0)/2.0,
	 std::sqrt(3.0*M_PI/10.0)/2.0},

	{3.0*std::sqrt(M_PI/5.0)/4.0,
	 -std::sqrt(M_PI/5.0)/2.0,
	 3.0*std::sqrt(M_PI/5.0)/4.0},

	{-std::sqrt(3.0*M_PI/10.0)/2.0,
	 std::sqrt(6.0*M_PI/5.0)/2.0,
	 std::sqrt(3.0*M_PI/10.0)/2.0},
									   
	{std::sqrt(3.0*M_PI/10.0)/4.0,
	 -std::sqrt(6.0*M_PI/5.0)/4.0,
	 std::sqrt(3.0*M_PI/10.0)/4.0}};

const double TidalDissipation::__torque_x_plus_coef[]={1.0, 
													   std::sqrt(1.5),
													   std::sqrt(1.5),
													   1.0,
													   0.0};

const double TidalDissipation::__torque_x_minus_coef[]={0.0, 
													    1.0,
													    std::sqrt(1.5),
													    std::sqrt(1.5),
													    1.0};

void TidalDissipation::fill_Umm(double inclination, bool deriv)
{
	double c=std::cos(inclination), s=std::sin(inclination),
		   s2=std::pow(s, 2), sc=s*c, cp1=c+1.0, cm1=c-1.0;

	__Umm[0][0]=__Umm_coef[0][0]*(deriv ? -2.0*s*cp1 : std::pow(cp1, 2));
	__Umm[1][0]=__Umm_coef[1][0]*(deriv ? cp1+2.0*s2 : s*cp1);
	__Umm[2][0]=__Umm_coef[2][0]*(deriv ? 2.0*sc : s2);
	__Umm[3][0]=-__Umm_coef[3][0]*(deriv ? c*cm1-s2 : s*cm1);
	__Umm[4][0]=__Umm_coef[4][0]*(deriv ? -2.0*s*cm1 : std::pow(cm1, 2));

	__Umm[0][1]=__Umm_coef[0][1]*(deriv ? 2.0*sc : s2);
	__Umm[1][1]=__Umm_coef[1][1]*(deriv ? 1.0-2.0*s2 : sc);
	__Umm[2][1]=__Umm_coef[2][1]*(deriv ? -6.0*sc : 2.0-3.0*s2);
	__Umm[3][1]=__Umm_coef[3][1]*(deriv ? 1.0-2.0*s2 : sc);
	__Umm[4][1]=__Umm_coef[4][1]*(deriv ? 2.0*sc : s2);

	__Umm[0][2]=__Umm_coef[0][2]*(deriv ? -2.0*cm1*s : std::pow(cm1, 2));
	__Umm[1][2]=-__Umm_coef[1][2]*(deriv ? c*cm1-s2 : s*cm1);
	__Umm[2][2]=__Umm_coef[2][2]*(deriv ? 2.0*sc : s2);
	__Umm[3][2]=__Umm_coef[3][2]*(deriv ? c*cp1-s2 : s*cp1);
	__Umm[4][2]=__Umm_coef[4][2]*(deriv ? -2.0*cp1*s: std::pow(cp1, 2));
}

void TidalDissipation::calculate_torque_power(const DissipatingBody &body,
		const SpinOrbitLockInfo &lock, short body_index,
		Dissipation::Derivative derivative)
{
	double &unlocked_power=rate_entry(body_index, Dissipation::POWER,
				   derivative, false),
		   &locked_power=rate_entry(body_index, Dissipation::POWER, 
				   derivative, true),
		   &unlocked_torque_x=rate_entry(body_index, Dissipation::TORQUEX,
				   derivative, false),
		   &locked_torque_x=rate_entry(body_index, Dissipation::TORQUEX,
				   derivative, true),
		   &unlocked_torque_z=rate_entry(body_index, Dissipation::TORQUEZ,
				   derivative, false),
		   &locked_torque_z=rate_entry(body_index, Dissipation::TORQUEZ, 
				   derivative, true);
	unlocked_torque_x=unlocked_torque_z=unlocked_power=locked_torque_x=
		locked_torque_z=locked_power=0;
	for(int m=-2; m<=2; ++m) {
		double m_spin_freq=m*body.current_spin();
		for(int mp=-2; mp<=2; mp+=2) {
			bool locked_term=lock(mp, m);
			double &power=(locked_term ? locked_power : unlocked_power),
				   &torque_x=(locked_term ? locked_torque_x :
						   unlocked_torque_x),
				   &torque_z=(locked_term ? locked_torque_z :
						   unlocked_torque_z);
			PhaseLag::Derivative phase_lag_deriv=PhaseLag::NO_DERIV;
			switch(derivative) {
				case Dissipation::NO_DERIV : case Dissipation::RADIUS :
				case Dissipation::INCLINATION :
					phase_lag_deriv=PhaseLag::NO_DERIV; break;
				case Dissipation::AGE : phase_lag_deriv=PhaseLag::AGE; break;
				case Dissipation::MOMENT_OF_INERTIA : 
				case Dissipation::SPIN_ANGMOM : 
				case Dissipation::SEMIMAJOR :
						   phase_lag_deriv=PhaseLag::FORCING_FREQUENCY;
						   break;
				default :
#ifdef DEBUG
						   assert(false)
#endif
						;
			};
			int mp_ind=mp/2+1, m_ind=m+2;
			double Umm_squared=std::pow(__Umm[m_ind][mp_ind], 2),
				   forcing_freq=(locked_term ? 0 :
						   mp*__orbital_frequency-m_spin_freq),
				   mod_phase_lag=body.modified_phase_lag(m, 
						   forcing_freq, lock, phase_lag_deriv);
			if(derivative==Dissipation::MOMENT_OF_INERTIA || 
					derivative==Dissipation::SPIN_ANGMOM)
				mod_phase_lag=
					-m*mod_phase_lag
					+
					body.modified_phase_lag(m, forcing_freq, lock,
						   PhaseLag::SPIN_FREQUENCY);
			else if(derivative==Dissipation::SEMIMAJOR) mod_phase_lag*=mp;
			torque_z+=Umm_squared*m*mod_phase_lag;
			power+=Umm_squared*mp*mod_phase_lag;
			torque_x+=__Umm[m_ind][mp_ind]*(
					(mp>-2 ? 
					 __torque_x_minus_coef[mp_ind]*__Umm[m_ind-1][mp_ind] : 0)+
					(mp<2 ? 
					 __torque_x_plus_coef[mp_ind]*__Umm[m_ind+1][mp_ind] : 0));
		}
	}
}

void TidalDissipation::calculate_semimajor_decay(short body_index)
{
	double coef=__mass_product/(2.0*std::pow(__orbital_energy, 2))*\
				AstroConst::G*AstroConst::solar_mass/AstroConst::AU*
				std::pow(AstroConst::day/AstroConst::solar_radius, 2);
	for(short locked=0; locked<=1; ++locked) {
		for(int fill_deriv=Dissipation::NO_DERIV;
				fill_deriv<Dissipation::NUM_DERIVATIVES; ++fill_deriv)
			rate_entry(body_index, Dissipation::SEMIMAJOR_DECAY, 
					static_cast<Dissipation::Derivative>(fill_deriv),
					locked)=
				coef*rate_entry(body_index, Dissipation::POWER, 
						static_cast<Dissipation::Derivative>(fill_deriv),
						locked);
		rate_entry(body_index, Dissipation::SEMIMAJOR_DECAY,
				Dissipation::SEMIMAJOR, locked)+=rate_entry(body_index,
					Dissipation::POWER, Dissipation::NO_DERIV, locked)/
					__orbital_energy;
	}
}

void TidalDissipation::calculate_inclination_decay(short body_index)
{
	double sin_inclination=std::sin(__inclination[body_index]),
		   cos_inclination=std::cos(__inclination[body_index]),
		   x_coef=1.0/__spin_angular_momentum[body_index]
			   +
			   cos_inclination/__orbital_angular_momentum,
		   z_coef=-sin_inclination/__orbital_angular_momentum;
	for(short locked=0; locked<=1; ++locked) {
		for(int fill_deriv=Dissipation::NO_DERIV;
				fill_deriv<Dissipation::NUM_DERIVATIVES; ++fill_deriv)
			rate_entry(body_index, Dissipation::INCLINATION_DECAY,
					static_cast<Dissipation::Derivative>(fill_deriv),
					locked)=x_coef*rate_entry(body_index,
						Dissipation::TORQUEX,
						static_cast<Dissipation::Derivative>(fill_deriv),
						locked)
				+
				z_coef*rate_entry(body_index,  Dissipation::TORQUEZ,
						static_cast<Dissipation::Derivative>(fill_deriv),
						locked);
		rate_entry(body_index, Dissipation::INCLINATION_DECAY,
				Dissipation::INCLINATION, locked)-=
			(sin_inclination*rate_entry(body_index, Dissipation::TORQUEX,
										Dissipation::NO_DERIV, locked)
			 +
			 cos_inclination*rate_entry(body_index, Dissipation::TORQUEZ,
				 Dissipation::NO_DERIV, locked)
			)/__orbital_angular_momentum;
		rate_entry(body_index, Dissipation::INCLINATION_DECAY,
				Dissipation::SPIN_ANGMOM, locked)-=rate_entry(body_index,
					Dissipation::TORQUEX, Dissipation::NO_DERIV, locked)/
					std::pow(__spin_angular_momentum[body_index], 2);
		rate_entry(body_index, Dissipation::INCLINATION_DECAY,
				Dissipation::SEMIMAJOR, locked)-=0.5/__semimajor*(
					cos_inclination/__orbital_angular_momentum
					+z_coef);
	}
}

void TidalDissipation::calculate_eccentricity_decay(short)
{
	throw Error::NotImplemented("Eccentric orbits");
}

TidalDissipation::TidalDissipation(const DissipatingBody &body1, 
		const DissipatingBody &body2, double semimajor, double eccentricity,
		const SpinOrbitLockInfo &lock1, const SpinOrbitLockInfo &lock2) :
	__orbital_frequency(std::sqrt(
			AstroConst::G*
			(body1.current_mass()+body2.current_mass())*
			AstroConst::solar_mass/
			std::pow(semimajor*AstroConst::AU, 3))*AstroConst::day),
	__mass_product(body1.current_mass()*body2.current_mass()),
	__reduced_mass(
			__mass_product/(body1.current_mass()+body2.current_mass())),
	__orbital_energy(-__mass_product/(2.0*semimajor)*
			AstroConst::G*AstroConst::solar_mass*
			std::pow(AstroConst::day/AstroConst::solar_radius, 2)/
			AstroConst::AU),
	__orbital_angular_momentum(__reduced_mass*__orbital_frequency*
			std::sqrt(1.0-std::pow(eccentricity, 2))*
			std::pow(semimajor/Rsun_AU, 2)),
	__Umm(std::valarray<double>(3), 5),
	__dissipation_rate(NaN, 
			2*Dissipation::NUM_QUANTITIES*Dissipation::NUM_DERIVATIVES),
	__spin_angular_momentum(2), __inclination(2)
{
	__spin_angular_momentum[0]=
		body1.current_moment_of_inertia()*body1.current_spin();
	__spin_angular_momentum[1]=
		body2.current_moment_of_inertia()*body2.current_spin();
	__inclination[0]=body1.current_inclination();
	__inclination[1]=body2.current_inclination();

	double torque_common=AstroConst::G*AstroConst::solar_mass*
			   			 std::pow(AstroConst::solar_radius/
								  std::pow(semimajor*AstroConst::AU, 2), 3)*
						 AstroConst::Gyr*AstroConst::day,
		   other_mass=body2.current_mass();
#ifdef DEBUG
	assert(Dissipation::INCLINATION==Dissipation::NUM_DERIVATIVES-1);
#endif
	for(short body_index=0; body_index<2; ++body_index) {
		const DissipatingBody &body=(body_index==0 ? body1 : body2);
		const SpinOrbitLockInfo &lock=(body_index==0 ? lock1 : lock2);
		double torque_scale=std::pow(body.current_radius(),5)*
			std::pow(other_mass, 2)*torque_common;
		fill_Umm(body.current_inclination(), false);
		for(int deriv_ind=Dissipation::NO_DERIV;
				deriv_ind<Dissipation::NUM_DERIVATIVES; ++deriv_ind) {
			Dissipation::Derivative derivative=
				static_cast<Dissipation::Derivative>(deriv_ind);
			if(derivative==Dissipation::RADIUS)
				torque_scale*=5.0/body.current_radius();
			else if(derivative==Dissipation::SEMIMAJOR)
				torque_scale*=-1.5/semimajor;
			else if(derivative==Dissipation::SPIN_ANGMOM)
				torque_scale/=body.current_moment_of_inertia();
			else if(derivative==Dissipation::MOMENT_OF_INERTIA)
				torque_scale*=-body.current_spin()/
					body.current_moment_of_inertia();
			else if(derivative==Dissipation::INCLINATION)
				fill_Umm(body.current_inclination(), true);
			calculate_torque_power(body, lock, body_index,
					derivative);
			for(short locked=0; locked<=1; ++locked) {
				rate_entry(body_index, Dissipation::POWER, derivative,
						locked)*=torque_scale*__orbital_frequency;
				rate_entry(body_index, Dissipation::TORQUEX, derivative,
						locked)*=torque_scale;
				rate_entry(body_index, Dissipation::TORQUEZ, derivative,
						locked)*=torque_scale;
				if(derivative==Dissipation::SEMIMAJOR) {
					rate_entry(body_index, Dissipation::POWER, derivative,
							locked)-=7.5/semimajor*rate_entry(body_index,
								Dissipation::POWER, Dissipation::NO_DERIV,
								locked);
					rate_entry(body_index, Dissipation::TORQUEX, derivative,
							locked)-=6.0/semimajor*rate_entry(body_index, 
								Dissipation::TORQUEX, Dissipation::NO_DERIV,
								locked);
					rate_entry(body_index, Dissipation::TORQUEZ, derivative,
							locked)-=6.0/semimajor*rate_entry(body_index, 
								Dissipation::TORQUEZ, Dissipation::NO_DERIV,
								locked);
				}
			}
		}
		other_mass=body1.current_mass();
	}
}

double TidalDissipation::operator()(short body_index, 
		Dissipation::Quantity quantity,
		Dissipation::Derivative derivative,
		bool locked)
{
	if(std::isnan(rate_entry(body_index, quantity, derivative, locked)))
		switch(quantity) {
			case Dissipation::SEMIMAJOR_DECAY : 
				calculate_semimajor_decay(body_index);
				break;
			case Dissipation::INCLINATION_DECAY : 
				calculate_inclination_decay(body_index);
				break;
			default :
#ifdef DEBUG
				assert(false)
#endif
					;
		}
	return rate_entry(body_index, quantity, derivative, locked);
}
