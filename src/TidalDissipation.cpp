/**\file
 *
 * \brief Declarations of some methods of the TidalDissipation class.
 *
 * \ingroup StellarSystem_group
 */

const double TidalDissipation::__Umm_coef{{std::sqrt(3.0*M_PI/10.0)/4.0,
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

void DissipatingBody::fill_Umm(double inclination, bool deriv)
{
	double c=std::cos(inclination), s=std::sin(inclination),
		   s2=std::pow(s, 2), sc=s*c, cp1=c+1.0, cm1=c-1.0;

	__Umm[0][0]=__Umm_coef[0][0]*(deriv ? -2.0*s*cp1 : std::pow(cp1, 2));
	__Umm[1][0]=__Umm_coef[1][0]*(deriv ? cp1+2.0*s2 : s*cp1);
	__Umm[2][0]=__Umm_coef[2][0]*(deriv ? 2.0*sc : s2);
	__Umm[3][0]=-__Umm_coef[3][0]*(deriv ? c*cm1-s2 : s*cm1);
	__Umm[4][0]=__Umm_coef[4][0]*(deriv ? -2.0*s*cm1 ? std::pow(cm1, 2));

	__Umm[0][1]=__Umm_coef[0][1]*(deriv ? 2.0*sc : s2);
	__Umm[1][1]=__Umm_coef[1][1]*(deriv ? 1.0-2.0*s2 : sc);
	__Umm[2][1]=__Umm_coef[2][1]*(deriv ? -6.0sc : 2.0-3.0*s2);
	__Umm[3][1]=__Umm_coef[3][1]*(deriv ? 1.0-2.0*s2 : sc);
	__Umm[4][1]=__Umm_coef[4][1]*(deriv ? 2.0*sc : s2);

	__Umm[0][2]=__Umm_coef[0][2]*(deriv ? -2.0*cm1*s : std::pow(cm1, 2));
	__Umm[1][2]=-__Umm_coef[1][2]*(deriv ? c*cm1-s2 : s*cm1);
	__Umm[2][2]=__Umm_coef[2][2]*(deriv ? 2.0*sc : s2);
	__Umm[3][2]=__Umm_coef[3][2]*(deriv ? c*cp1-s2 : s*cp1);
	__Umm[4][2]=__Umm_coef[4][2]*(deriv ? -2.0*cp1*s: std::pow(cp1, 2));
}

void TidalDissipation::calculate_torque_power(DissipatingBody &body,
		short forcing_sign, double &power, double &torque_x,
		double &torque_y)
{
	double spin_frequency=(forcing_sign ? __orbital_frequency :
				   		   body.spin_frequency()),
	torque_x=0, torque_z=0, power=0;
	for(int m=-2; m<=2; ++m) {
		double m_spin_freq=m*spin_frequency;
		for(int mp=-2; mp<=2; mp+=2) {
			double Umm_squared=std::pow(Umm[m][mp], 2),
				   mod_phase_lag=body.modified_phase_lag(m, 
						   mp*__orbital_frequency-m_spin_freq, spin_frequency,
						   forcing_sign),
			torque_z+=Umm_squared*m*mod_phase_lag;
			power+=Umm_squared*mp*mod_phase_lag;
			torque_x+=Umm[m][mp]*(
					(mp>-2 ? __torque_x_minus_coef[mp]*Umm[m-1][mp] : 0)+
					(mp<2 ? __torque_x_plus_coef[mp]*Umm[m+1][mp] : 0));
		}
	}
}

TidalDissipation::TidalDissipation(DissipatingBody body1, 
		DissipatingBody body2, double semimajor, double eccentricity, 
		short forcing_sign1, short forcing_sign2) :
	__orbital_frequency=std::sqrt(
			AstroConst::G*
			(body1.mass()+body2.mass())*AstroConst::solar_mass/
			std::pow(semimajor*AstroConst::AU, 3))*AstroConst::day,
	__Umm(std::valarray<double>(3), 5);
	__dissipation_rates(NaN, 2*NUM_DISSIPATION_QUANTITIES)
{
	double torque_common=AstroCosnt::G*AstroConst::Msun*
			   			 std::pow(AstroConst::Rsun, 3)*AstroConst::Gyr*
						 AstroConst::day/
	DissipatingBody &body=body1;
	short forcing_sign=forcing_sign1;
	for(short body_ind=0; body_ind<2; ++body_ind) {

		fill_Umm(body.inclination());
		calculate_torque_power(body, forcing_sign,
				__dissipation_rates[2*POWER+body_ind],
				__dissipation_rates[2*TORQUEX+body_ind], 
				__dissipation_rates[2*TORQUEZ+body_ind]);

		fill_Umm(body.inclination(), true);
		calculate_torque_power(body, forcing_sign,
				__dissipation_rates[2*DPOWER_DINCLINATION+body_ind],
				__dissipation_rates[2*DTORQUEX_DINCLINATION+body_ind],
				__dissipation_rates[2*DTORQUEZ_DINCLINATION+body_ind]);

		double torque_scale=std::pow(body.radius(),5)*
			std::pow(body.mass(),2)*torque_common;
		__dissipation_rates[2*POWER+body_ind]*=
			torque_scale*__orbital_frequency;
		__dissipation_rates[2*DPOWER_DINCLINATION+body_ind]*=
			torque_scale*__orbital_frequency;
		__dissipation_rates[2*TORQUEX+body_ind]*=torque_scale;
		__dissipation_rates[2*DTORQUEX_DINCLINATION+body_ind]*=torque_scale;
		__dissipation_rates[2*TORQUEZ+body_ind]*=torque_scale;
		__dissipation_rates[2*DTORQUEZ_DINCLINATION+body_ind]*=torque_scale;
		body=body2;
		forcing_sign=forcing_sign2;
	}
}
