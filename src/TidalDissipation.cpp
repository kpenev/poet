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

void DissipatingBody::fill_Umm(double inclination, 
		std::valarray< std::valarray<double> > &Umm)
{
	double c=std::cos(inclination), s=std::sin(inclination),
		   s2=std::pow(s, 2), sc=s*c;

	Umm[0][0]=__Umm_coef[0][0]*std::pow(1.0+c, 2);
	Umm[1][0]=__Umm_coef[1][0]*(s+sc);
	Umm[2][0]=__Umm_coef[2][0]*s2;
	Umm[3][0]=__Umm_coef[3][0]*(s-sc);
	Umm[4][0]=__Umm_coef[4][0]*std::pow(1.0-c, 2);

	Umm[0][1]=__Umm_coef[0][1]*s2;
	Umm[1][1]=__Umm_coef[1][1]*sc;
	Umm[2][1]=__Umm_coef[2][1]*(2.0-3.0*s2);
	Umm[3][1]=__Umm_coef[3][1]*sc;
	Umm[4][1]=__Umm_coef[4][1]*s2;

	Umm[0][2]=__Umm_coef[0][2]*std::pow(1.0-c, 2);
	Umm[1][2]=__Umm_coef[1][2]*(s-sc);
	Umm[2][2]=__Umm_coef[2][2]*s2;
	Umm[3][2]=__Umm_coef[3][2]*(s+sc);
	Umm[4][2]=__Umm_coef[4][2]*std::pow(1.0+c, 2);
}

TidalDissipation::TidalDissipation(DissipatingBody body1, 
		DissipatingBody body2, double semimajor, double eccentricity, 
		short forcing_sign1, short forcing_sign2)
{
	std::valarray< std::valarray<double> > Umm1(std::valarray<double>(3), 5),
										   Umm2(std::valarray<double>(3), 5);
	fill_Umm(body1.inclination(), Umm1);
	fill_Umm(body2.inclination(), Umm2);
	double orbital_frequency=std::sqrt(
			AstroConst::G*
			(body1.mass()+body2.mass())*AstroConst::solar_mass/
			std::pow(semimajor*AstroConst::AU, 3))*AstroConst::day,
		   torque_common=AstroCosnt::G*AstroConst::Msun*
			   			 std::pow(AstroConst::Rsun, 3)*AstroConst::Gyr*
						 AstroConst::day/
						 std::pow(semimajor*AstroConst::AU, 6),
		   torque1_scale=std::pow(body1.radius(),5)*std::pow(body1.mass(),2)*
			   			 torque_common,
		   torque2_scale=std::pow(body2.radius(),5)*std::pow(body2.mass(),2)*
			   			 torque_common,
		   spin_frequency1=(forcing_sign1 ? orbital_frequency :
				   			body1.spin_frequency()),
		   spin_frequency2=(forcing_sign1 ? orbital_frequency :
				   			body1.spin_frequency());
		   __torque_x=0, __torque_z=0, __power=0;
	for(int m=-2; m<=2; ++m) {
		double m_orbit_freq=m*orbital_frequency;
		for(int mp=-2; mp<=2; mp+=2) {
			double Umm1_squared=std::pow(Umm1[m][mp], 2),
				   Umm2_squared=std::pow(Umm2[m][mp], 2),
				   mod_phase_lag1=modified_phase_lag(m, 
						   mp*orbital_frequency-m_orbit_freq, spin_frequency1,
						   forcing_sign1),
				   mod_phase_lag2=modified_phase_lag(m, 
						   mp*orbital_frequency-m_orbit_freq, spin_frequency2,
						   forcing_sign2);
			__torque1_z+=Umm1_squared*m*mod_phase_lag1;
			__torque2_z+=Umm2_squared*m*mod_phase_lag2;
			__power1+=Umm1_squared*mp*mod_phase_lag1;
			__power2+=Umm2_squared*mp*mod_phase_lag2;
			__torque1_x+=Umm1[m][mp]*(
					(mp>-2 ? __torque_x_minus_coef[mp]*Umm1[m-1][mp] : 0)+
					(mp<2 ? __torque_x_plus_coef[mp]*Umm1[m+1][mp] : 0));
			__torque2_x+=Umm2[m][mp]*(
					(mp>-2 ? __torque_x_minus_coef[mp]*Umm2[m-1][mp] : 0)+
					(mp<2 ? __torque_x_plus_coef[mp]*Umm2[m+1][mp] : 0));
		}
	}
	power1*=T0*orbital_frequency;
	power2*=T0*orbital_frequency;
	torque1_x*=T0*AstroConst::Gyr;
	torque2_x*=T0*AstroConst::Gyr;
	torque1_z*=T0*AstroConst::Gyr;
	torque2_z*=T0*AstroConst::Gyr;
}
