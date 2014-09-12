#include "DissipatingZone.h"

EccentricityExpansionCoefficients DissipatingZone::__pms;

const double DissipatingZone::__Umm_coef[][3]={
	{std::sqrt(3.0*M_PI/10.0)/4.0,
	 -std::sqrt(6.0*M_PI/5.0)/4.0,
	 std::sqrt(3.0*M_PI/10.0)/4.0},

	{-std::sqrt(3.0*M_PI/10.0)/2.0,
	 -std::sqrt(6.0*M_PI/5.0)/2.0,
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

const double DissipatingZone::__torque_x_plus_coef[]={1.0, 
													  std::sqrt(1.5),
													  std::sqrt(1.5),
													  1.0,
													  0.0};

const double DissipatingZone::__torque_x_minus_coef[]={0.0,           //m=-2
													   1.0,           //m=-1
													   std::sqrt(1.5),//m=0
													   std::sqrt(1.5),//m=1
													   1.0};          //m=2

void DissipatingZone::fill_Umm()
{
	if(__Ummp_inclination==inclination()) return;
	__Ummp_inclination=inclination();
	double c=std::cos(__Ummp_inclination), s=std::sin(__Ummp_inclination),
		   s2=std::pow(s, 2), sc=s*c, cp1=c+1.0, cm1=c-1.0;

	__Ummp[0][0]=__Umm_coef[0][0]*std::pow(cp1, 2);
	__Ummp_deriv[0][0]=-__Umm_coef[0][0]*2.0*s*cp1;

	__Ummp[1][0]=__Umm_coef[1][0]*s*cp1;
	__Ummp_deriv[1][0]=__Umm_coef[1][0]*(cp1+2.0*s2);

	__Ummp[2][0]=__Umm_coef[2][0]*s2;
	__Ummp_deriv[2][0]=__Umm_coef[2][0]*2.0*sc;

	__Ummp[3][0]=-__Umm_coef[3][0]*s*cm1;
	__Ummp_deriv[3][0]=-__Umm_coef[3][0]*(c*cm1-s2);

	__Ummp[4][0]=__Umm_coef[4][0]*std::pow(cm1, 2);
	__Ummp_deriv[4][0]=-__Umm_coef[4][0]*2.0*s*cm1;



	__Ummp[0][1]=__Umm_coef[0][1]*s2;
	__Ummp_deriv[0][1]=__Umm_coef[0][1]*2.0*sc;

	__Ummp[1][1]=__Umm_coef[1][1]*sc;
	__Ummp_deriv[1][1]=__Umm_coef[1][1]*(1.0-2.0*s2);

	__Ummp[2][1]=__Umm_coef[2][1]*(2.0-3.0*s2);
	__Ummp_deriv[2][1]=-__Umm_coef[2][1]*6.0*sc;

	__Ummp[3][1]=__Umm_coef[3][1]*sc;
	__Ummp_deriv[3][1]=__Umm_coef[3][1]*(1.0-2.0*s2);

	__Ummp[4][1]=__Umm_coef[4][1]*s2;
	__Ummp_deriv[4][1]=__Umm_coef[4][1]*2.0*sc;



	__Ummp[0][2]=__Umm_coef[0][2]*std::pow(cm1, 2);
	__Ummp_deriv[0][2]=-__Umm_coef[0][2]*2.0*cm1*s;

	__Ummp[1][2]=-__Umm_coef[1][2]*s*cm1;
	__Ummp_deriv[1][2]=-__Umm_coef[1][2]*(c*cm1-s2);

	__Ummp[2][2]=__Umm_coef[2][2]*s2;
	__Ummp_deriv[2][2]=__Umm_coef[2][2]*2.0*sc;

	__Ummp[3][2]=__Umm_coef[3][2]*s*cp1;
	__Ummp_deriv[3][2]=__Umm_coef[3][2]*(c*cp1-s2);

	__Ummp[4][2]=__Umm_coef[4][2]*std::pow(cp1, 2);
	__Ummp_deriv[4][2]=-__Umm_coef[4][2]*2.0*cp1*s;
}

void DissipatingZone::potential_term(double e, int m, int mp, 
		unsigned e_order, double &no_deriv, double &inclination_deriv, 
		double &eccentricity_deriv)
{
	no_deriv=inclination_deriv=eccentricity_deriv=0;
	for(int i=0; i<3; ++i) {
		double pms=__pms(2*(i-1), mp, e, e_order, false);
		no_deriv+=pms*__Ummp[m+2][i];
		inclination_deriv+=pms*__Ummp_deriv[m+2][i];
		eccentricity_deriv+=
			__pms(2*(i-1), mp, e, e_order, true)*__Ummp[m+2][i];
	}
}

DissipatingZone::DissipatingZone()
	: __Ummp_inclination(NaN), __Ummp(5), __Ummp_deriv(5),
	__power(2*Dissipation::END_DIMENSIONLESS_DERIV),
	__torque_x(2*Dissipation::END_DIMENSIONLESS_DERIV),
	__torque_y(2*Dissipation::END_DIMENSIONLESS_DERIV),
	__torque_z(2*Dissipation::END_DIMENSIONLESS_DERIV)
{
	for(int i=0; i<5; ++i) {
		__Ummp[i].resize(3);
		__Ummp_deriv[i].resize(3);
	}
}

void DissipatingZone::set_orbit(double orbital_frequency,
		double eccentricity)
{
	fill_Umm();
	__power=0;
	__torque_x=0;
	__torque_y=0;
	__torque_z=0;
	int e_order=eccentricity_order(eccentricity);

	for(int mp=-e_order-2; mp<=e_order+2; ++mp) {
		double U_mm1mp_value=0, U_mm1mp_i_deriv=0, U_mm1mp_e_deriv=0,
			   U_mp1mp_value, U_mp1mp_i_deriv, U_mp1mp_e_deriv;
		potential_term(eccentricity, -2, mp, e_order,  U_mp1mp_value,
				U_mp1mp_i_deriv, U_mp1mp_e_deriv);
		for(int m=-2; m<=2; ++m) {
			int m_ind=m+2;
			bool locked_term=locked(mp, m);
			double U_mmp_value=U_mp1mp_value, 
				   U_mmp_i_deriv=U_mp1mp_i_deriv, 
				   U_mmp_e_deriv=U_mp1mp_e_deriv;
			if(m<2)
				potential_term(eccentricity, m+1, mp, e_order, U_mp1mp_value,
						U_mp1mp_i_deriv, U_mp1mp_e_deriv);
			else U_mp1mp_value=U_mp1mp_i_deriv=U_mp1mp_e_deriv=0;

			for(int deriv=Dissipation::NO_DERIV; 
					(m!=0 || mp!=0) && 
					deriv<Dissipation::END_DIMENSIONLESS_DERIV; ++deriv) {
				Dissipation::Derivative phase_lag_deriv=
					(deriv<Dissipation::END_PHASE_LAG_DERIV
					 ? static_cast<Dissipation::Derivative>(deriv)
					 : Dissipation::NO_DERIV);
				double mod_phase_lag_above,
					   mod_phase_lag_below=modified_phase_lag(mp, m,
							   orbital_frequency, phase_lag_deriv,
							   mod_phase_lag_above),
					   love_coef=love_coefficient(mp, m,
							   (phase_lag_deriv==Derivative::AGE
								? Derivative::age : Derivative::NO_DERIV));
					   U_mmp, U_mp1mp, U_mm1mp; 
				if(deriv<Dissipation::END_PHASE_LAG_DERIV) {
					U_mmp=U_mmp_value;
					U_mp1mp=U_mp1mp_value;
					U_mm1mp=U_mm1mp_value;
				} else if(deriv==Dissipation::INCLINATION) {
					U_mmp=U_mmp_i_deriv;
					U_mp1mp=U_mp1mp_i_deriv;
					U_mm1mp=U_mm1mp_i_deriv;
				} else {
					U_mmp=U_mmp_e_deriv;
					U_mp1mp=U_mp1mp_e_deriv;
					U_mm1mp=U_mm1mp_e_deriv;
				}
				double U_mmp_squared=std::pow(U_mmp, 2),
					   term_power=U_mmp_squared*mp,
					   term_torque_z=U_mmp_squared*m,
					   term_torque_x=U_mmp*(
								__torque_x_minus_coef[m_ind]*U_mm1mp+
								__torque_x_plus_coef[m_ind]*U_mp1mp);
				if(!locked_term) mod_phase_lag_above=mod_phase_lag_below;
				int deriv_ind=2*deriv;
				__power[deriv_ind]+=term_power*mod_phase_lag_below;
				__torque_z[deriv_ind]+=term_torque_z*mod_phase_lag_below;
				__torque_x[deriv_ind]+=term_torque_x*mod_phase_lag_below;
				__torque_y[deriv_ind+1]=
					-(__torque_y[deriv_ind]-=term_torque_x*love_coef);
				__power[deriv_ind+1]+=term_power*mod_phase_lag_above;
				__torque_z[deriv_ind+1]+=term_torque_z*mod_phase_lag_above;
				__torque_x[deriv_ind+1]+=term_torque_x*mod_phase_lag_above;
			}
			U_mm1mp_value=U_mmp_value;
			U_mm1mp_i_deriv=U_mmp_i_deriv;
			U_mm1mp_e_deriv=U_mmp_e_deriv;
		}
	}
}
