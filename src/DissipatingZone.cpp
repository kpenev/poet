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
	if(__Ummp_inclination==__inclination) return;
	__Ummp_inclination=__inclination;
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
		double &no_deriv, double &inclination_deriv,
		double &eccentricity_deriv)
{
	no_deriv=inclination_deriv=eccentricity_deriv=0;
	for(int i=0; i<3; ++i) {
		double pms=__pms(2*(i-1), mp, e, __e_order, false);
		no_deriv+=pms*__Ummp[m+2][i];
		inclination_deriv+=pms*__Ummp_deriv[m+2][i];
		eccentricity_deriv+=
			__pms(2*(i-1), mp, e, __e_order, true)*__Ummp[m+2][i];
	}
}

void DissipatingZone::fix_forcing_frequency(const SpinOrbitLockInfo &limit, 
		int orbital_frequency_multiplier, int spin_frequency_multiplier
		double &forcing_frequency)
{
#ifdef DEBUG
	assert(limit.lock_direction());
#endif
	int expected_sign=limit.spin_frequency_multiplier()*
		(orbital_frequency_multiplier*limit.spin_frequency_multiplier()
		 -
		 limit.orbital_frequency_multiplier()*spin_frequency_multiplier);
	if(expected_sign*forcing_frequency>0) return;
	if(limit.lock_direction()*spin_frequency_multiplier*expected_sign<0)
		forcing_frequency=(expected_sign>0
						   ? std::numeric_limits::epsilon()
						   : -std::numeric_limits::epsilon());
}

double DissipatingZone::forcing_frequency(int orbital_frequency_multiplier,
		int spin_frequency_multiplier, double orbital_frequency)
{
#ifdef DEBUG
	assert(__lock || 
			(__lock.lock_direction()*__other_lock.lock_direction()==-1));
#endif
	if(__lock(orbital_frequency_multiplier, spin_frequency_multiplier))
		return 0;
	double forcing_freq=orbital_frequency_multiplier*orbital_frequency
						-
						spin_frequnecy_multiplier*spin_frequency;
	fix_forcing_frequency(__lock, forcing_freq);
	fix_forcing_frequency(__other_lock, forcing_freq);
}

void update_lock_to_lower_e_order(SpinOrbitLockInfo &lock)
{
#ifdef DEBUG
	assert(lock.lock_direction());
#endif
	if(std::abs(lock.orbital_frequency_multiplier())>__e_order+2
			&& lock.spin_frequency_multiplier==2)
		lock.set_lock((lock.orbital_frequency_multiplier()
					   -
					   lock.lock_direction())/2,
					  1, lock.lock_direction());
	if(lock.orbital_frequency_multiplier()>__e_order+2) {
		if(lock.lock_direction()>0)
			lock.set_lock(__e_order+2, 1, 1);
		else lock.set_lock(1, 0, 1);
	} else if(lock.orbital_frequency_multiplier()<-__e_order-2) {
		if(lock.lock_direction()>0) lock.set_lock(1, 0, 1)
		else lock.set_lock(-__e_order-2, 1, -1);
	}
}

DissipatingZone::DissipatingZone()
	: __e_order(0), __Ummp_inclination(NaN), __Ummp(5), __Ummp_deriv(5),
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

void DissipatingZone::configure(double age, double orbital_frequency,
		double eccentricity, double orbital_angmom, double spin_angmom,
		double inclination, double periapsis)
{
	__orbital_angmom=orbital_angom;
	__orbital_frequnecy=orbital_frequency;
	if(__lock) {
		__spin_frequency=__lock.spin(orbital_frequency);
		__angular_momentum=__spin_frequency*moment_of_inertia();
	} else {
		__angular_momentum=spin_angmom; 
		__spin_frequency=spin_angmom/moment_of_inertia();
	}
	__inclination=inclination;
	__periapsis=periapsis;
	if(std::isnan(orbital_frequency)) return;

	fill_Umm();
	__power=0;
	__torque_x=0;
	__torque_y=0;
	__torque_z=0;

	for(int mp=-__e_order-2; mp<=__e_order+2; ++mp) {
		double U_mm1mp_value=0, U_mm1mp_i_deriv=0, U_mm1mp_e_deriv=0,
			   U_mp1mp_value, U_mp1mp_i_deriv, U_mp1mp_e_deriv;
		potential_term(eccentricity, -2, mp, U_mp1mp_value, U_mp1mp_i_deriv,
				U_mp1mp_e_deriv);
		for(int m=-2; m<=2; ++m) {
			int m_ind=m+2;
			bool locked_term=locked(mp, m);
			double U_mmp_value=U_mp1mp_value, 
				   U_mmp_i_deriv=U_mp1mp_i_deriv, 
				   U_mmp_e_deriv=U_mp1mp_e_deriv;
			if(m<2)
				potential_term(eccentricity, m+1, mp, U_mp1mp_value,
						U_mp1mp_i_deriv, U_mp1mp_e_deriv);
			else U_mp1mp_value=U_mp1mp_i_deriv=U_mp1mp_e_deriv=0;

			for(int deriv=Dissipation::NO_DERIV; 
					(m!=0 || mp!=0) && 
					deriv<Dissipation::END_DIMENSIONLESS_DERIV; ++deriv) {
				Dissipation::Derivative phase_lag_deriv=
					(deriv<Dissipation::END_PHASE_LAG_DERIV
					 ? static_cast<Dissipation::Derivative>(deriv)
					 : Dissipation::NO_DERIV);
				double tidal_frequency=forcing_frequency(mp, m),
					   mod_phase_lag_above,
					   mod_phase_lag_below=modified_phase_lag(mp, m,
							   tidal_frequency, phase_lag_deriv,
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
				if(!locked_term && (tidal_frequency!=0 ||
									__lock.lock_direction()<0))
					mod_phase_lag_above=mod_phase_lag_below;
				else if(!locked_term && tidal_frequency==0
						&& __lock.lock_direction()>0)
					mod_phase_lag_below=mod_phase_lag_above;
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

double DissipatingZone::periapsis_evolution(
		const Eigen::Vector3d &orbit_torque,
		const Eigen::Vector3d &zone_torque,
		Dissipation::Derivative deriv=Dissipation::NO_DERIV,
		const Eigen::Vector3d &orbit_torque_deriv,
		const Eigen::Vector3d &zone_torque_deriv)
{
	double sin_inc=std::sin(__inclination),
		   cos_inc=std::cos(__inclination),
		   zone_y_torque, orbit_y_torque;
	if(deriv==Dissipation::NO_DERIV) {
		orbit_y_torque=orbit_torque[2];
		zone_y_torque=zone_torque[2];
	} else {
		orbit_y_torque=orbit_torque_deriv[2];
		zone_y_torque=zone_torque_deriv[2];
	}
	double result=-orbit_y_torque*cos_inc/(__orbital_angmom*sin_inc)
				  +
				  zone_y_torque/(angular_momentum()*sin_inc);
	if(		deriv==Dissipation::NO_DERIV 
			|| deriv==Dissipation::AGE 
			|| deriv==Dissipation::ECCENTRICITY
			|| deriv==Dissipation::PERIAPSIS
			|| deriv==Dissipation::RADIUS
			|| deriv==Dissipation::MOMENT_OF_INERTIA
			|| deriv==Dissipation::SEMIMAJOR)
		return result;
	else if(deriv==Dissipation::SPIN_FREQUENCY ||
			deriv==Dissipation::SPIN_ANGMOM)
		return result
			   -
			   zone_torque[2]/(std::pow(angular_momentum(), 2)*sin_inc)
			   *(deriv==Dissipation::SPIN_FREQUENCY ? moment_of_inertia():1);
	else if(deriv==INCLINATION) 
		return result
			   -
			   (
					orbit_torque[2]/__orbital_angmom
					+
					zone_torque[2]*cos_inc/angular_momentum()
			   )/std::pow(sin_inc, 2);
#ifdef DEBUG
	else assert(false);
#endif
}

double DissipatingZone::inclination_evolution(
		const Eigen::Vector3d &orbit_torque,
		const Eigen::Vector3d &zone_torque,
		Dissipation::Derivative deriv=Dissipation::NO_DERIV,
		const Eigen::Vector3d &orbit_torque_deriv,
		const Eigen::Vector3d &zone_torque_deriv)
{
	double sin_inc=std::sin(__inclination),
		   cos_inc=std::cos(__inclination),
		   zone_x_torque,
		   orbit_x_torque,
		   orbit_z_torque;
	if(deriv==Dissipation::NO_DERIV) {
		orbit_x_torque=orbit_torque[1];
		orbit_z_torque=orbit_torque[3]
		zone_x_torque=zone_torque[1];
	} else {
		orbit_x_torque=orbit_torque_deriv[1];
		orbit_z_torque=orbit_torque_deriv[3];
		zone_x_torque=zone_torque_deriv[1];
	}
	double result=(orbit_x_torque*cos_inc-orbit_z_torque*sin_inc)
			      /__orbital_angmom
				  -
				  zone_x_torque/angular_momentum();
	if(		deriv==Dissipation::NO_DERIV 
			|| deriv==Dissipation::AGE 
			|| deriv==Dissipation::ECCENTRICITY
			|| deriv==Dissipation::PERIAPSIS
			|| deriv==Dissipation::RADIUS
			|| deriv==Dissipation::MOMENT_OF_INERTIA
			|| deriv==Dissipation::SEMIMAJOR)
		return result;
	else if(deriv==Dissipation::SPIN_FREQUENCY ||
			deriv==Dissipation::SPIN_ANGMOM)
		return result
			   +
			   zone_torque[1]/std::pow(angular_momentum(), 2)
			   *(deriv==Dissipation::SPIN_FREQUENCY ? moment_of_inertia():1);
	else if(deriv==INCLINATION)
		return result
			   +
			   (orbit_torque[3]*cos_inc + orbit_torque[1]*sin_inc)
			   /angular_momentum();
#ifdef DEBUG
	else assert(false);
#endif
}

void DissipatingZone::release_lock()
{
	if(__lock.spin_frequency_multiplier()==2) {
#ifdef DEBUG
		assert(__lock.orbital_frequency_multiplier()%2==1);
#endif
		__other_lock.set_lock((__lock.orbital_frequency_multiplier()+1)/2,
				1, -1);
		__lock.set_lock((__lock.orbital_frequency_multiplier()-1)/2, 1, -1);
	}
	update_lock_to_lower_e_order(__lock);
	update_lock_to_lower_e_order(__other_lock);
	if(__lock.spin_frequency_multiplier()==0) {
		__lock=__other_lock;
		__other_lock.set_lock(1, 0, 1);
	}
}

void DissipatingZone::release_lock(short direction)
{
#ifdef DEBUG
	assert(__lock);
	assert(direction==1 || direction==-1);
#endif
	__lock.lock_direction(direction);
	int orbit_mult=2*__lock.orbital_frequency_multiplier()+direction;
	if(orbit_mult%2) __other_lock.set_lock(orbit_mult, 2, -direction);
	else __other_lock.set_lock(orbit_mult/2, 1, -direction);
	update_lock_to_lower_e_order(__other_lock);
}

Eigen::Vector3D zone_to_zone_transform(const DissipatingZone &from_zone,
		const DissipationgZone &to_zone, const Eigen::Vector3D &vector,
		Dissipation::Derivative deriv, bool with_respect_to_from)
{
#ifdef DEBUG
	assert(deriv==Dissipation::NO_DERIV || deriv==Dissipation::INCLINATION
			|| deriv==Dissipation::PERIAPSIS);
#endif
	double cos_dw=std::cos(to_zone.periapsis()-from_zone.periapsis()),
		   sin_dw=std::sin(to_zone.periapsis()-from_zone.periapsis()),
		   sin_from=std::sin(from_zone.inclination()),
		   cos_from=std::cos(from_zone.inclination()),
		   sin_to=std::sin(to_zone.inclination()),
		   cos_to=std::cos(to_zone.inclination()),
		   sin_sin, cos_cos, sin_cos, cos_sin;
	if(deriv==Dissipation::PERIAPSIS) {
		double buffer=cos_dw;
		if(with_respect_to_from) {cos_dw=sin_dw; sin_dw=-buffer;}
		else {cos_dw=-sin_dw; sin_dw=buffer;}
	} else if(deriv==Dissipation::INCLINATION) {
		if(with_respect_to_from) {
			sin_sin=sin_to*cos_from;
			cos_cos=-cos_to*sin_from;
			sin_cos=-sin_to*sin_from;
			cos_sin=cos_to*cos_from;
		} else {
			sin_sin=cos_to*sin_from;
			cos_cos=-sin_to*cos_from;
			sin_cos=cos_to*cos_from;
			cos_sin=-sin_to*sin_from;
		}
	} else {
		sin_sin=sin_to*sin_from;
		cos_cos=cos_to*cos_from;
		sin_cos=sin_to*cos_from;
		cos_sin=cos_to*sin_from;
	}
	return Eigen::Vector3D(
			(sin_sin+cos_cos*cos_dw)*vector[0] - cos_to*sin_dw*vector[1]
			+ (sin_cos-cos_sin*cos_dw)*vector[2],


			cos_from*sin_dw*vector[0] + cos_dw*vector[1]
			- sin_from*sin_dw*vector[2],

			(cos_sin-sin_cos*cos_dw)*vector[0] + sin_to*sin_dw*vector[1]
			+ (cos_cos+sin_sin*cos_dw)*vector[2]);
}

void change_e_order(unsigned new_e_order)
{
	if(__lock) {
	   if(__lock.orbital_frequency_multiplier()>__e_order+2) release_lock();
	   __e_order=new_e_order;
	   return;
	}
#ifdef DEBUG
	assert(__lock.lock_direction()*__other_lock.lock_direction()<0);
	assert((__lock.spin_frequency_multiplier()==1 ||
				__lock.spin_frequency_multiplier()==2));
	assert((__other_lock.spin_frequency_multiplier()>=0 && 
				__other_lock.spin_frequency_multiplier()<=2));
	assert(__lock.lock_direction()*__lock.spin_frequency_multiplier()
			*spin_frequency()
			>
			__lock.lock_direction()*__lock.orbital_frequency_multiplier()*
			__orbital_frequency);
	if(__other_lock.spin_frequency_multiplier())
		assert(__other_lock.lock_direction()
				*__other_lock.spin_frequency_multiplier()
				*spin_frequency()
				>
				__other_lock.lock_direction()
				*__other_lock.orbital_frequency_multiplier()
				*__orbital_frequency);
	else assert(__lock.lock_direction()
				*__lock.orbital_frequency_multiplier()
				>0);
#endif
	if(new_e_order>__e_order) {
		if(__other_lock.spin_frequency_miltiplier()==0) {
			int orb_mult=std::ceil(std::abs(2*spin_frequency())
								/__orbital_frequency);
			if(orb_mult%2==1 && new_e_order+2>=orb_mult)
				__other_lock.set_lock(__lock.lock_direction()*orb_mult, 2,
									  -__lock.lock_direction());
			else if(new_e_order+2>=orb_mult/2)
				__other_lock.set_lock(__lock.lock_direction()*orb_mult/2, 1,
									  -__lock.lock_direction());
		} else if(__lock.spin_frequency_multiplier()==1
				&& __other_lock.spin_frequency_multiplier()==1
				&& __lock.orbital_frequency_multiplier()
			   +
			   __other_lock.orbital_frequency_multiplier()<=new_e_order+2) {
			int mid_mult=__lock.orbital_frequency_multiplier()
				+
				__other_lock.orbital_frequency_multiplier();
			if(2*spin_frequency()<mid_mult*__orbital_frequency)
				(__lock.lock_direction()>0 ? __other_lock : __lock)
					.set_lock(mid_mult, 2, -1);
			else (__lock.lock_direction()>0 ? __lock : __other_lock)
				.set_lock(mid_mult, 2, 1);
		}
		__e_order=new_e_order;
	} else {
		__e_order=new_e_order;
		update_lock_to_lower_e_order(__lock);
		update_lock_to_lower_e_order(__other_lock);
		if(__lock.spin_frequency_multiplier()==0) {
			__lock=__other_lock;
			__other_lock.set_lock(1, 0, 1);
		}
	}
}

