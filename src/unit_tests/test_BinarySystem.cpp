#include "test_BinarySystem.h"

void test_BinarySystem::test_orbit_diff_eq(
		RandomDiskPlanetSystem &system,
		const std::valarray<double> &expected, bool diff_eq)
{
	std::valarray<double> returned_orbit, *to_compare;
	EvolModeType evol_mode=system().fill_orbit(returned_orbit);
	if(diff_eq) {
		to_compare=new std::valarray<double>(expected.size());
		system().differential_equations(system().age(), &(returned_orbit[0]),
				evol_mode, &((*to_compare)[0]));
	} else to_compare=&returned_orbit;
	std::ostringstream msg;
	msg << "The " << (diff_eq ? "differential equations" : "orbit")
		<< " created by the binary system " << (diff_eq ? "have" : "has")
		<< " a different size (" << to_compare->size() << ") than expected (" 
		<< expected.size() << ")";
	TEST_ASSERT_MSG(to_compare->size()==expected.size(), msg.str().c_str());
	for(unsigned i=0; i<expected.size(); ++i) {
		msg.str("");
		msg << "System:" << std::endl
			<< system << "Expected "
			<< (diff_eq ? "differential equation[" : "orbit[")
			<< i << "]=" << expected[i] << ", got: " << (*to_compare)[i]
			<< ", difference: " << (*to_compare)[i] - expected[i];
		TEST_ASSERT_MSG(check_diff((*to_compare)[i], expected[i],
								   1e-10, 1e-15), msg.str().c_str());
	}
	if(diff_eq) delete to_compare;
}

Lags test_BinarySystem::signed_lags(const RandomDiskPlanetSystem &system,
		unsigned zone_ind) const
{
	using namespace SystemParameters;
	Lags result;
	for(int mp=-2; mp<=2; ++mp) 
		for(int m=-2; m<=2; ++m) {
			if(system.lock(zone_ind)(mp, m)) result(m, mp)=0;
			else {
				double forcing_freq=
					system.orbital_frequency()*mp
					-
					system.quantity(
						static_cast<Quantity>(FIRST_ZONE_ANGVEL+zone_ind))*m;
				result(m, mp)=(forcing_freq<0 ? -1.0 : 1.0)
							  *system.lags(zone_ind)(m, mp);
			}
		}
	return result;
}

Lags test_BinarySystem::locked_lags_below(
		const RandomDiskPlanetSystem &system, unsigned zone_ind) const
{
	using namespace SystemParameters;
	Lags result;
	assert(system.lock(zone_ind));
	for(int mp=-2; mp<=2; ++mp) 
		for(int m=-2; m<=2; ++m)
			if(system.lock(zone_ind)(mp, m)) 
				result(m, mp)=(m<0 ? 1 : -1)*system.lags(zone_ind)(m, mp);
			else result(m, mp)=0;
	return result;
}

void test_BinarySystem::fill_torques_angvel_in_orbit_coord(
		const RandomDiskPlanetSystem &system,
		std::vector<Eigen::Vector2d> &nonlocked_tidal_torques,
		std::vector<Eigen::Vector2d> &angular_velocities,
		Eigen::Vector2d &locked_tidal_torque) const
{
	using namespace SystemParameters;
	for(unsigned zone_ind=0; zone_ind<4; ++zone_ind) {
		Lags lags=signed_lags(system, zone_ind);
		double inclination=system.quantity(
						static_cast<Quantity>(FIRST_INCLINATION+zone_ind)),
			   x_torque=dimensionless_torque_x_Lai(inclination, lags),
			   z_torque=dimensionless_torque_z_Lai(inclination, lags),
			   torque_norm;
		if(zone_ind<2)
			torque_norm=torque_norm_Lai(
					system.quantity(SECONDARY_MASS), 
					system.quantity(PRIMARY_RADIUS),
					system.quantity(SEMIMAJOR));
		else 
			torque_norm=torque_norm_Lai(
					system.quantity(PRIMARY_MASS), 
					system.quantity(SECONDARY_RADIUS),
					system.quantity(SEMIMAJOR));
		double sin_inc=std::sin(inclination),
			   cos_inc=std::cos(inclination);
		x_torque*=torque_norm;
		z_torque*=torque_norm;
		nonlocked_tidal_torques[zone_ind](0)=
			-x_torque*cos_inc+z_torque*sin_inc;
		nonlocked_tidal_torques[zone_ind](1)=
			x_torque*sin_inc+z_torque*cos_inc;
		if(system.lock(zone_ind)) {
			lags=locked_lags_below(system, zone_ind);
			x_torque=
				torque_norm*dimensionless_torque_x_Lai(inclination, lags);
			z_torque=
				torque_norm*dimensionless_torque_z_Lai(inclination, lags);
			locked_tidal_torque(0)=-x_torque*cos_inc+z_torque*sin_inc;
			locked_tidal_torque(1)=x_torque*sin_inc+z_torque*cos_inc;
		}
		angular_velocities[zone_ind]=
			system.quantity(
					static_cast<Quantity>(FIRST_ZONE_ANGVEL+zone_ind))
			*
			Eigen::Vector2d(sin_inc, cos_inc);
	}
}

void test_BinarySystem::fill_nontidal_torques_in_orbit_coord(
		const RandomDiskPlanetSystem &system,
		const std::vector<Eigen::Vector2d> &angular_velocities,
		std::vector<Eigen::Vector2d> &nontidal_torques) const
{
	using namespace SystemParameters;
	for(unsigned body_ind=0; body_ind<2; ++body_ind) {
		//Wind
		nontidal_torques[2*body_ind]=
			-system.quantity(
					static_cast<Quantity>(FIRST_WIND_STRENGTH+body_ind))
			*angular_velocities[2*body_ind]
			*std::pow(std::min(
						system.quantity(
							static_cast<Quantity>(FIRST_ZONE_ANGVEL
												  +2*body_ind)),
						system.quantity(
							static_cast<Quantity>(FIRST_WIND_SAT_FREQ
												  +body_ind))
						),
					2)
			*std::sqrt(system.quantity(
						static_cast<Quantity>(FIRST_ZONE_RADIUS+2*body_ind))
					   /
					   system.quantity(
						static_cast<Quantity>(FIRST_ZONE_MASS+2*body_ind)));

		//Core growth
		double core_growth=system.quantity(
				static_cast<Quantity>(FIRST_CORE_MASS_DERIV+body_ind));
		nontidal_torques[2*body_ind+1]=
			2.0/3.0*core_growth
			*std::pow(system.quantity(
				static_cast<Quantity>(FIRST_ZONE_RADIUS+2*body_ind+1)), 2)
			*angular_velocities[2*body_ind+(core_growth>0 ? 0 : 1)];

		//Differential rotation torque
		double 
			inertia1=system.quantity(
					static_cast<Quantity>(FIRST_ZONE_INERTIA+2*body_ind)),
			inertia2=system.quantity(
					static_cast<Quantity>(FIRST_ZONE_INERTIA+2*body_ind+1));
		nontidal_torques[2*body_ind+1]+=
			(inertia1*inertia2)/(inertia1+inertia2)
			*(
					angular_velocities[2*body_ind]
					-
					angular_velocities[2*body_ind+1]
			 )/system.quantity(
				 static_cast<Quantity>(FIRST_COUPLING_TIMESCALE+body_ind));

		//Add core growth and differential rotation torque to envelopes
		nontidal_torques[2*body_ind]-=nontidal_torques[2*body_ind+1];
	}
}

void test_BinarySystem::fill_diff_eq(const RandomDiskPlanetSystem &system,
		const Eigen::Vector2d &orbit_torque,
		const std::vector<Eigen::Vector2d> &zone_torques,
		std::valarray<double> &expected_diff_eq)
{
	using namespace SystemParameters;
	double m1=system.quantity(PRIMARY_MASS),
		   m2=system.quantity(SECONDARY_MASS),
		   a=system.quantity(SEMIMAJOR),
		   worb=system.orbital_frequency(),
		   orbital_angmom=(m1*m2)/(m1+m2)*std::pow(a,2)*worb,
		   orbit_rotation=orbit_torque(0)/orbital_angmom;
	expected_diff_eq[0]=2.0*orbit_torque(1)*a/orbital_angmom;
	unsigned locked_found=0;
	for(unsigned zone_ind=0; zone_ind<4; ++zone_ind) {
		double sin_inc=std::sin(system.quantity(
					static_cast<Quantity>(FIRST_INCLINATION+zone_ind))),
			   cos_inc=std::cos(system.quantity(
					static_cast<Quantity>(FIRST_INCLINATION+zone_ind)));
		expected_diff_eq[2+zone_ind]=
			(
			 zone_torques[zone_ind](0)*cos_inc
			 -
			 zone_torques[zone_ind](1)*sin_inc
			)
			/system.quantity(static_cast<Quantity>(FIRST_ZONE_ANGVEL
					 							   +zone_ind))
			/system.quantity(static_cast<Quantity>(FIRST_ZONE_INERTIA
												   +zone_ind))
			-
			orbit_rotation;
		if(system.lock(zone_ind)) ++locked_found;
		else {
			expected_diff_eq[9+zone_ind-locked_found]=
				zone_torques[zone_ind](0)*sin_inc
				+
				zone_torques[zone_ind](1)*cos_inc;
		}
	}
	if(locked_found==0) expected_diff_eq[0]*=6.5*std::pow(a, 5.5);
}

void test_BinarySystem::test_fill_orbit_locked_surface()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(LOCKED_SURFACE_SPIN);
		std::valarray<double> expected_orbit(1);
		expected_orbit[0]=system_maker.quantity(PRIMARY_CORE_INERTIA)
						  *system_maker.quantity(PRIMARY_ANGVEL_CORE);
		test_orbit_diff_eq(system_maker, expected_orbit);
	}
}

void test_BinarySystem::test_fill_orbit_single()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(SINGLE);
		std::valarray<double> expected_orbit(4);
		expected_orbit[0]=system_maker.quantity(PRIMARY_INCLINATION_CORE);
		expected_orbit[1]=system_maker.quantity(PRIMARY_PERIAPSIS_CORE);
		expected_orbit[2]=system_maker.quantity(PRIMARY_ENV_INERTIA)
						  *system_maker.quantity(PRIMARY_ANGVEL_ENV);
		expected_orbit[3]=system_maker.quantity(PRIMARY_CORE_INERTIA)
						  *system_maker.quantity(PRIMARY_ANGVEL_CORE);
		test_orbit_diff_eq(system_maker, expected_orbit);
	}
}

void test_BinarySystem::test_fill_orbit_binary_no_locks()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(BINARY, 0, 0);
		std::valarray<double> expected_orbit(13);
		expected_orbit[0]=std::pow(system_maker.quantity(SEMIMAJOR), 6.5);
		expected_orbit[1]=system_maker.quantity(ECCENTRICITY);
		expected_orbit[2]=system_maker.quantity(PRIMARY_INCLINATION_ENV);
		expected_orbit[3]=system_maker.quantity(PRIMARY_INCLINATION_CORE);
		expected_orbit[4]=system_maker.quantity(SECONDARY_INCLINATION_ENV);
		expected_orbit[5]=system_maker.quantity(SECONDARY_INCLINATION_CORE);
		expected_orbit[6]=system_maker.quantity(PRIMARY_PERIAPSIS_CORE);
		expected_orbit[7]=system_maker.quantity(SECONDARY_PERIAPSIS_ENV);
		expected_orbit[8]=system_maker.quantity(SECONDARY_PERIAPSIS_CORE);
		expected_orbit[9]=system_maker.quantity(PRIMARY_ENV_INERTIA)
						  *system_maker.quantity(PRIMARY_ANGVEL_ENV);
		expected_orbit[10]=system_maker.quantity(PRIMARY_CORE_INERTIA)
						   *system_maker.quantity(PRIMARY_ANGVEL_CORE);
		expected_orbit[11]=system_maker.quantity(SECONDARY_ENV_INERTIA)
						  *system_maker.quantity(SECONDARY_ANGVEL_ENV);
		expected_orbit[12]=system_maker.quantity(SECONDARY_CORE_INERTIA)
						   *system_maker.quantity(SECONDARY_ANGVEL_CORE);
		test_orbit_diff_eq(system_maker, expected_orbit);
	}
}

void test_BinarySystem::test_fill_orbit_binary_locks()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		std::cerr << "\rStarting test#" << i;
		std::cerr.flush();
		RandomDiskPlanetSystem system_maker(BINARY, 1, 1);
		std::valarray<double> 
			expected_orbit(13-system_maker.num_locked_zones());
		unsigned ind=0;
		expected_orbit[ind++]=system_maker.quantity(SEMIMAJOR);
		expected_orbit[ind++]=system_maker.quantity(ECCENTRICITY);
		expected_orbit[ind++]=system_maker.quantity(PRIMARY_INCLINATION_ENV);
		expected_orbit[ind++]=system_maker.quantity(
				PRIMARY_INCLINATION_CORE);
		expected_orbit[ind++]=
			system_maker.quantity(SECONDARY_INCLINATION_ENV);
		expected_orbit[ind++]=
			system_maker.quantity(SECONDARY_INCLINATION_CORE);
		expected_orbit[ind++]=system_maker.quantity(PRIMARY_PERIAPSIS_CORE);
		expected_orbit[ind++]=system_maker.quantity(SECONDARY_PERIAPSIS_ENV);
		expected_orbit[ind++]=
			system_maker.quantity(SECONDARY_PERIAPSIS_CORE);
		if(!system_maker.lock(0))
			expected_orbit[ind++]=system_maker.quantity(PRIMARY_ENV_INERTIA)
								  *system_maker.quantity(PRIMARY_ANGVEL_ENV);
		if(!system_maker.lock(1))
			expected_orbit[ind++]=system_maker.quantity(PRIMARY_CORE_INERTIA)
								  *
								  system_maker.quantity(PRIMARY_ANGVEL_CORE);
		if(!system_maker.lock(2))
			expected_orbit[ind++]=
				system_maker.quantity(SECONDARY_ENV_INERTIA)
				*system_maker.quantity(SECONDARY_ANGVEL_ENV);
		if(!system_maker.lock(3))
			expected_orbit[ind++]=
				system_maker.quantity(SECONDARY_CORE_INERTIA)
				*system_maker.quantity(SECONDARY_ANGVEL_CORE);
		test_orbit_diff_eq(system_maker, expected_orbit);
	}
	std::cerr << std::endl;
}

void test_BinarySystem::test_locked_surface_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(LOCKED_SURFACE_SPIN);
		std::valarray<double> expected_diff_eq(1);
		double Icore=system_maker.quantity(PRIMARY_CORE_INERTIA),
			   Ienv=system_maker.quantity(PRIMARY_ENV_INERTIA),
			   core_growth=system_maker.quantity(PRIMARY_CORE_MASS_DERIV);
		expected_diff_eq[0]=
			Icore*Ienv/(Icore+Ienv)
			/system_maker.quantity(PRIMARY_COUPLING_TIMESCALE)
			*(system_maker.quantity(DISK_LOCK_FREQ)
			  -
			  system_maker.quantity(PRIMARY_ANGVEL_CORE))
			+
			2.0/3.0*core_growth*
			std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2)*
			(core_growth>0 ? system_maker.quantity(DISK_LOCK_FREQ)
			 			   : system_maker.quantity(PRIMARY_ANGVEL_CORE));
		test_orbit_diff_eq(system_maker, expected_diff_eq, true);
	}
}

void test_BinarySystem::test_single_aligned_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(SINGLE, 0, 0, true, true, true,
				true, true, true, true, true);
		std::valarray<double> expected_diff_eq(0.0, 4);
		double Icore=system_maker.quantity(PRIMARY_CORE_INERTIA),
			   Ienv=system_maker.quantity(PRIMARY_ENV_INERTIA),
			   angmom_loss=system_maker.quantity(PRIMARY_WIND_STRENGTH)
						   *system_maker.quantity(PRIMARY_ANGVEL_ENV)
						   *std::pow(std::min(
								system_maker.quantity(PRIMARY_ANGVEL_ENV),
								system_maker.quantity(PRIMARY_WIND_SAT_FREQ))
								,
								2)
						   *std::sqrt(system_maker.quantity(PRIMARY_RADIUS)
									  /system_maker.quantity(PRIMARY_MASS));
		expected_diff_eq[3]=
			Icore*Ienv/(Icore+Ienv)
			/system_maker.quantity(PRIMARY_COUPLING_TIMESCALE)
			*(system_maker.quantity(PRIMARY_ANGVEL_ENV)
			  -
			  system_maker.quantity(PRIMARY_ANGVEL_CORE))
			+
			2.0/3.0*system_maker.quantity(PRIMARY_CORE_MASS_DERIV)
			*std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2)
			*(system_maker.quantity(PRIMARY_CORE_MASS_DERIV)>0
			  ? system_maker.quantity(PRIMARY_ANGVEL_ENV)
			  : system_maker.quantity(PRIMARY_ANGVEL_CORE));
		expected_diff_eq[2]=-expected_diff_eq[3]-angmom_loss;

		test_orbit_diff_eq(system_maker, expected_diff_eq, true);
	}
}

void test_BinarySystem::test_single_zero_periapsis_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(SINGLE, 0, 0, true, false, true);
		std::valarray<double> expected_diff_eq(0.0, 4);
		double Icore=system_maker.quantity(PRIMARY_CORE_INERTIA),
			   Ienv=system_maker.quantity(PRIMARY_ENV_INERTIA),
			   angmom_loss=system_maker.quantity(PRIMARY_WIND_STRENGTH)
						   *system_maker.quantity(PRIMARY_ANGVEL_ENV)
						   *std::pow(std::min(
								system_maker.quantity(PRIMARY_ANGVEL_ENV),
								system_maker.quantity(PRIMARY_WIND_SAT_FREQ))
								,
								2)
						   *std::sqrt(system_maker.quantity(PRIMARY_RADIUS)
									  /system_maker.quantity(PRIMARY_MASS)),
			   sin_inc=
				   std::sin(system_maker.quantity(PRIMARY_INCLINATION_CORE)),
			   cos_inc=
				   std::cos(system_maker.quantity(PRIMARY_INCLINATION_CORE)),
			   abs_angmom_env=Ienv*system_maker.quantity(PRIMARY_ANGVEL_ENV),
			   abs_angmom_core=
				   Icore*system_maker.quantity(PRIMARY_ANGVEL_CORE);
		Eigen::Vector2d 
			angmom_env(0, abs_angmom_env),
			angmom_core(sin_inc*abs_angmom_core, cos_inc*abs_angmom_core),
			coup_rot=(Ienv*angmom_core-Icore*angmom_env)
					 /system_maker.quantity(PRIMARY_COUPLING_TIMESCALE)
					 /(Icore+Ienv)
					 -
					 2.0/3.0*system_maker.quantity(PRIMARY_CORE_MASS_DERIV)
					 *std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2)
					 *(system_maker.quantity(PRIMARY_CORE_MASS_DERIV)>0
					   ? angmom_env/Ienv
					   : angmom_core/Icore);
		expected_diff_eq[0]=-coup_rot(0)*(1.0/abs_angmom_env
										  +
										  cos_inc/abs_angmom_core)
							+coup_rot(1)*sin_inc/abs_angmom_core;
		expected_diff_eq[2]=-angmom_loss+coup_rot[1];
		expected_diff_eq[3]=
			-coup_rot.dot(Eigen::RowVector2d(sin_inc, cos_inc));

		test_orbit_diff_eq(system_maker, expected_diff_eq, true);
	}
}

void test_BinarySystem::test_binary_no_locks_circular_aligned_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(BINARY, 0, 0, true, true, true,
											true, true, true, true, true);
		std::valarray<double> expected_diff_eq(0.0, 13);
		double orbit_power=0,
			   m1=system_maker.quantity(PRIMARY_MASS),
			   m2=system_maker.quantity(SECONDARY_MASS),
			   r1=system_maker.quantity(PRIMARY_RADIUS),
			   r2=system_maker.quantity(SECONDARY_RADIUS),
			   a=system_maker.quantity(SEMIMAJOR);
		for(unsigned zone_ind=0; zone_ind<4; ++zone_ind) {
			double power_norm, torque_norm;
			if(zone_ind<2) {
				power_norm=power_norm_Lai(m1, m2, r1, a);
				torque_norm=torque_norm_Lai(m2, r1, a);
			} else {
				power_norm=power_norm_Lai(m2, m1, r2, a);
				torque_norm=torque_norm_Lai(m1, r2, a);
			}
			Lags zone_lags=signed_lags(system_maker, zone_ind);
			orbit_power-=power_norm*dimensionless_power_Lai(0, zone_lags);
			expected_diff_eq[9+zone_ind]=
				torque_norm*dimensionless_torque_z_Lai(0, zone_lags);
		}
		double w1=system_maker.quantity(PRIMARY_ANGVEL_ENV),
			   w2=system_maker.quantity(PRIMARY_ANGVEL_CORE),
			   w3=system_maker.quantity(SECONDARY_ANGVEL_ENV),
			   w4=system_maker.quantity(SECONDARY_ANGVEL_CORE),
			   wsat1=system_maker.quantity(PRIMARY_WIND_SAT_FREQ),
			   wsat2=system_maker.quantity(SECONDARY_WIND_SAT_FREQ),
			   primary_wind_loss=system_maker.quantity(PRIMARY_WIND_STRENGTH)
				   				 *w1*std::pow(std::min(w1, wsat1), 2)
								 *std::sqrt(
									system_maker.quantity(PRIMARY_RADIUS)
									/system_maker.quantity(PRIMARY_MASS)),
			   secondary_wind_loss=
				   system_maker.quantity(SECONDARY_WIND_STRENGTH)
				   *w3*std::pow(std::min(w3, wsat2), 2)
				   *std::sqrt(system_maker.quantity(SECONDARY_RADIUS)
						   	  /system_maker.quantity(SECONDARY_MASS)),
			   i1=system_maker.quantity(PRIMARY_ENV_INERTIA),
			   i2=system_maker.quantity(PRIMARY_CORE_INERTIA),
			   i3=system_maker.quantity(SECONDARY_ENV_INERTIA),
			   i4=system_maker.quantity(SECONDARY_CORE_INERTIA),
			   primary_coupling_torque=
				   (i1*i2)/(i1+i2)*(w2-w1)
				   /system_maker.quantity(PRIMARY_COUPLING_TIMESCALE),
			   secondary_coupling_torque=
				   (i3*i4)/(i3+i4)*(w4-w3)
				   /system_maker.quantity(SECONDARY_COUPLING_TIMESCALE);
		if(system_maker.quantity(PRIMARY_CORE_MASS_DERIV)>0)
			primary_coupling_torque-=
				2.0/3.0*w1*system_maker.quantity(PRIMARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2);
		else 
			primary_coupling_torque-=
				2.0/3.0*w2*system_maker.quantity(PRIMARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2);
		if(system_maker.quantity(SECONDARY_CORE_MASS_DERIV)>0)
			secondary_coupling_torque-=
				2.0/3.0*w3*system_maker.quantity(SECONDARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(SECONDARY_CORE_RADIUS), 2);
		else
			secondary_coupling_torque-=
				2.0/3.0*w4*system_maker.quantity(SECONDARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(SECONDARY_CORE_RADIUS), 2);
		expected_diff_eq[0]=13.0*std::pow(a, 7.5)*orbit_power/(m1*m2)
							*std::pow(AstroConst::solar_radius, 3)
							/AstroConst::G/AstroConst::solar_mass
							/std::pow(AstroConst::day, 2);
		expected_diff_eq[9]+=primary_coupling_torque-primary_wind_loss;
		expected_diff_eq[10]-=primary_coupling_torque;
		expected_diff_eq[11]+=secondary_coupling_torque-secondary_wind_loss;
		expected_diff_eq[12]-=secondary_coupling_torque;
		test_orbit_diff_eq(system_maker, expected_diff_eq, true);
	}
}

void test_BinarySystem::test_binary_no_locks_circular_inclined_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(BINARY, 0, 0, true, false, false,
											true, false, false, true, true);
		std::valarray<double> expected_diff_eq(0.0, 13);
		double orbit_power=0,
			   m1=system_maker.quantity(PRIMARY_MASS),
			   m2=system_maker.quantity(SECONDARY_MASS),
			   r1=system_maker.quantity(PRIMARY_RADIUS),
			   r2=system_maker.quantity(SECONDARY_RADIUS),
			   a=system_maker.quantity(SEMIMAJOR),
			   orbital_energy=-(m1*m2)/(2.0*a)
							  *AstroConst::G*AstroConst::solar_mass
							  *std::pow(AstroConst::day, 2)
				   			  /std::pow(AstroConst::solar_radius, 3),
			   worb=std::sqrt(AstroConst::G*(m1+m2)*AstroConst::solar_mass
					   		  *std::pow(AstroConst::day, 2)
							  /std::pow(a*AstroConst::solar_radius, 3)),
			   orbital_angmom=(m1*m2)/(m1+m2)*std::pow(a,2)*worb,
			   orbit_rotation=0;
		for(unsigned zone_ind=0; zone_ind<4; ++zone_ind) {
			const Lags &lag=signed_lags(system_maker, zone_ind);
			double power_norm, torque_norm,
				   inclination=system_maker.quantity(
						  static_cast<Quantity>(FIRST_INCLINATION+zone_ind)),
				   sin_inc=std::sin(inclination),
				   cos_inc=std::cos(inclination);
			if(zone_ind<2) {
				power_norm=power_norm_Lai(m1, m2, r1, a);
				torque_norm=torque_norm_Lai(m2, r1, a);
			} else {
				power_norm=power_norm_Lai(m2, m1, r2, a);
				torque_norm=torque_norm_Lai(m1, r2, a);
			}
			double x_torque=torque_norm
							*dimensionless_torque_x_Lai(inclination, lag),
				   z_torque=torque_norm
							*dimensionless_torque_z_Lai(inclination, lag);
			orbit_power-=power_norm*dimensionless_power_Lai(inclination,lag);
			orbit_rotation+=z_torque*sin_inc-x_torque*cos_inc;
			expected_diff_eq[9+zone_ind]=
				torque_norm*dimensionless_torque_z_Lai(inclination, lag);
			double zone_rotation=
				-x_torque/system_maker.quantity(
						static_cast<Quantity>(FIRST_ZONE_ANGVEL+zone_ind))
				/system_maker.quantity(
						static_cast<Quantity>(FIRST_ZONE_INERTIA+zone_ind));
			expected_diff_eq[2+zone_ind]=zone_rotation;
		}
		orbit_rotation/=orbital_angmom;
		double i1=system_maker.quantity(PRIMARY_ENV_INERTIA),
			   i2=system_maker.quantity(PRIMARY_CORE_INERTIA),
			   i3=system_maker.quantity(SECONDARY_ENV_INERTIA),
			   i4=system_maker.quantity(SECONDARY_CORE_INERTIA),
			   wsat1=system_maker.quantity(PRIMARY_WIND_SAT_FREQ),
			   wsat2=system_maker.quantity(SECONDARY_WIND_SAT_FREQ),
			   w1=system_maker.quantity(PRIMARY_ANGVEL_ENV),
			   w2=system_maker.quantity(PRIMARY_ANGVEL_CORE),
			   w3=system_maker.quantity(SECONDARY_ANGVEL_ENV),
			   w4=system_maker.quantity(SECONDARY_ANGVEL_CORE),
			   core_growth1=system_maker.quantity(PRIMARY_CORE_MASS_DERIV),
			   core_growth2=system_maker.quantity(SECONDARY_CORE_MASS_DERIV);
		expected_diff_eq[9]-=system_maker.quantity(PRIMARY_WIND_STRENGTH)
							 *w1*std::pow(std::min(w1, wsat1), 2)
							 *std::sqrt(r1/m1);
		expected_diff_eq[11]-=system_maker.quantity(SECONDARY_WIND_STRENGTH)
							  *w3*std::pow(std::min(w3, wsat2), 2)
							  *std::sqrt(r2/m2);
		Eigen::Vector2d 
			w1_vec=w1*Eigen::Vector2d(
				std::sin(system_maker.quantity(PRIMARY_INCLINATION_ENV)),
				std::cos(system_maker.quantity(PRIMARY_INCLINATION_ENV))),
			w2_vec=w2*Eigen::Vector2d(
				std::sin(system_maker.quantity(PRIMARY_INCLINATION_CORE)),
				std::cos(system_maker.quantity(PRIMARY_INCLINATION_CORE))),
			w3_vec=w3*Eigen::Vector2d(
				std::sin(system_maker.quantity(SECONDARY_INCLINATION_ENV)),
				std::cos(system_maker.quantity(SECONDARY_INCLINATION_ENV))),
			w4_vec=w4*Eigen::Vector2d(
				std::sin(system_maker.quantity(SECONDARY_INCLINATION_CORE)),
				std::cos(system_maker.quantity(SECONDARY_INCLINATION_CORE))),
			primary_coupling_torque=
				(i1*i2)/(i1+i2)*(w2_vec-w1_vec)
				/system_maker.quantity(PRIMARY_COUPLING_TIMESCALE),
			secondary_coupling_torque=
				(i3*i4)/(i3+i4)*(w4_vec-w3_vec)
				/system_maker.quantity(SECONDARY_COUPLING_TIMESCALE);
		primary_coupling_torque-=2.0/3.0*core_growth1
			*std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2)
			*(core_growth1>0 ? w1_vec : w2_vec);
		secondary_coupling_torque-=2.0/3.0*core_growth2
			*std::pow(system_maker.quantity(SECONDARY_CORE_RADIUS), 2)
			*(core_growth2>0 ? w3_vec : w4_vec);
		expected_diff_eq[9]+=primary_coupling_torque.dot(w1_vec)/w1;
		expected_diff_eq[10]-=primary_coupling_torque.dot(w2_vec)/w2;
		expected_diff_eq[11]+=secondary_coupling_torque.dot(w3_vec)/w3;
		expected_diff_eq[12]-=secondary_coupling_torque.dot(w4_vec)/w4;
		expected_diff_eq[2]+=
			(primary_coupling_torque(0)*std::cos(
					system_maker.quantity(PRIMARY_INCLINATION_ENV))
			 -
			 primary_coupling_torque(1)*std::sin(
					system_maker.quantity(PRIMARY_INCLINATION_ENV)))
			/system_maker.quantity(static_cast<Quantity>(PRIMARY_ANGVEL_ENV))
			/system_maker.quantity(
					static_cast<Quantity>(PRIMARY_ENV_INERTIA))
			+
			orbit_rotation;
		expected_diff_eq[3]+=
			(-primary_coupling_torque(0)*std::cos(
					system_maker.quantity(PRIMARY_INCLINATION_CORE))
			 +
			 primary_coupling_torque(1)*std::sin(
					system_maker.quantity(PRIMARY_INCLINATION_CORE)))
			/system_maker.quantity(
					static_cast<Quantity>(PRIMARY_ANGVEL_CORE))
			/system_maker.quantity(
					static_cast<Quantity>(PRIMARY_CORE_INERTIA))
			+
			orbit_rotation;
		expected_diff_eq[4]+=
			(secondary_coupling_torque(0)*std::cos(
					system_maker.quantity(SECONDARY_INCLINATION_ENV))
			 -
			 secondary_coupling_torque(1)*std::sin(
					system_maker.quantity(SECONDARY_INCLINATION_ENV)))
			/system_maker.quantity(
					static_cast<Quantity>(SECONDARY_ANGVEL_ENV))
			/system_maker.quantity(
					static_cast<Quantity>(SECONDARY_ENV_INERTIA))
			+
			orbit_rotation;
		expected_diff_eq[5]+=
			(-secondary_coupling_torque(0)*std::cos(
					system_maker.quantity(SECONDARY_INCLINATION_CORE))
			 +
			 secondary_coupling_torque(1)*std::sin(
					system_maker.quantity(SECONDARY_INCLINATION_CORE)))
			/system_maker.quantity(
					static_cast<Quantity>(SECONDARY_ANGVEL_CORE))
			/system_maker.quantity(
					static_cast<Quantity>(SECONDARY_CORE_INERTIA))
			+
			orbit_rotation;
		expected_diff_eq[0]=-6.5*orbit_power*std::pow(a, 6.5)/orbital_energy;
		test_orbit_diff_eq(system_maker, expected_diff_eq, true);
	}
}

void test_BinarySystem::test_binary_1lock_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(BINARY, 1, 1, true, false, false,
											true, false, false, true, true);
		std::valarray<double> expected_diff_eq(0.0, 13);
		std::vector<Eigen::Vector2d> nonlocked_tidal_torques(4),
									 nontidal_torques(4),
									 angular_velocities(4);
		Eigen::Vector2d locked_tidal_torque;
		fill_torques_angvel_in_orbit_coord(system_maker,
										   nonlocked_tidal_torques,
										   angular_velocities,
										   locked_tidal_torque);
		fill_nontidal_torques_in_orbit_coord(system_maker,angular_velocities,
											 nontidal_torques);
		Eigen::Vector2d orbit_torque(0, 0);
		int locked_zone_ind=-1;
		for(unsigned zone_ind=0; zone_ind<4; ++zone_ind) {
			if(system_maker.lock(zone_ind)) locked_zone_ind=zone_ind;
			orbit_torque-=nonlocked_tidal_torques[zone_ind];
		}
		Eigen::Vector2d locked_spin_dir;
		if(locked_zone_ind>=0)
			locked_spin_dir=angular_velocities[locked_zone_ind]
							/angular_velocities[locked_zone_ind].norm();
		double m1=system_maker.quantity(PRIMARY_MASS),
			   m2=system_maker.quantity(SECONDARY_MASS),
			   a=system_maker.quantity(SEMIMAJOR),
			   locked_inertia=(locked_zone_ind>=0
					   		   ? system_maker.quantity(
								   static_cast<Quantity>(FIRST_ZONE_INERTIA
									   					 +locked_zone_ind))
							   : NaN),
			   orbit_coef=3.0*(m1+m2)/(m1*m2*std::pow(a, 2));
		std::vector<Eigen::Vector2d> zone_torques(4);
		for(unsigned zone_ind=0; zone_ind<4; ++zone_ind)
			zone_torques[zone_ind]=nonlocked_tidal_torques[zone_ind]
								   +
								   nontidal_torques[zone_ind];
		if(locked_zone_ind>=0) {
			double lock_coef=
				(
				 zone_torques[locked_zone_ind].dot(locked_spin_dir)
				 /locked_inertia
				 -
				 system_maker.quantity(
					 static_cast<Quantity>(FIRST_ZONE_ANGVEL
						 				   +locked_zone_ind))
				 *system_maker.quantity(
					 static_cast<Quantity>(FIRST_ZONE_INERTIA_DERIV
						 				   +locked_zone_ind))
				 /locked_inertia
				 +
				 orbit_coef*orbit_torque(1)
				)
				/
				(
				 orbit_coef*locked_tidal_torque(1)
				 -
				 locked_tidal_torque.dot(locked_spin_dir)/locked_inertia
				);
			zone_torques[locked_zone_ind]+=lock_coef*locked_tidal_torque;
			orbit_torque-=lock_coef*locked_tidal_torque;
		}
		fill_diff_eq(system_maker, orbit_torque, zone_torques,
					 expected_diff_eq);
		test_orbit_diff_eq(system_maker, expected_diff_eq, true);
	}
}

test_BinarySystem::test_BinarySystem(unsigned ntests,
			const std::string &eccentricity_expansion) : __ntests(ntests)
{
	DissipatingZone::read_eccentricity_expansion(eccentricity_expansion);
/*	TEST_ADD(test_BinarySystem::test_fill_orbit_locked_surface);
	TEST_ADD(test_BinarySystem::test_fill_orbit_single);
	TEST_ADD(test_BinarySystem::test_fill_orbit_binary_no_locks);
	TEST_ADD(test_BinarySystem::test_fill_orbit_binary_locks);
	TEST_ADD(test_BinarySystem::test_locked_surface_diff_eq);
	TEST_ADD(test_BinarySystem::test_single_aligned_diff_eq);
	TEST_ADD(test_BinarySystem::test_single_zero_periapsis_diff_eq);
	TEST_ADD(
		test_BinarySystem::test_binary_no_locks_circular_aligned_diff_eq);
	TEST_ADD(
		test_BinarySystem::test_binary_no_locks_circular_inclined_diff_eq);*/
	TEST_ADD(test_BinarySystem::test_binary_1lock_diff_eq);
}

#ifdef STANDALONE
int main()
{
	std::srand(7);
	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_BinarySystem tests(1);
	return (tests.run(output, true) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
