#include "test_BinarySystem.h"

void test_BinarySystem::test_orbit_diff_eq(RandomDiskPlanetSystem &system,
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
								   1e-10,1e-15), msg.str().c_str());
	}
	if(diff_eq) delete to_compare;
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
		RandomDiskPlanetSystem system_maker(BINARY, 1, 4);
		std::valarray<double> 
			expected_orbit(13-system_maker.num_locked_zones());
		unsigned ind=0;
		expected_orbit[ind++]=system_maker.quantity(SEMIMAJOR);
		expected_orbit[ind++]=system_maker.quantity(ECCENTRICITY);
		expected_orbit[ind++]=system_maker.quantity(PRIMARY_INCLINATION_ENV);
		expected_orbit[ind++]=system_maker.quantity(PRIMARY_INCLINATION_CORE);
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
}

void test_BinarySystem::test_locked_surface_diff_eq()
{
	using namespace SystemParameters;
	for(unsigned i=0; i<__ntests; ++i) {
		RandomDiskPlanetSystem system_maker(LOCKED_SURFACE_SPIN);
		std::valarray<double> expected_diff_eq(1);
		double Icore=system_maker.quantity(PRIMARY_CORE_INERTIA),
			   Ienv=system_maker.quantity(PRIMARY_ENV_INERTIA);
		expected_diff_eq[0]=
			Icore*Ienv/(Icore+Ienv)
			/system_maker.quantity(PRIMARY_COUPLING_TIMESCALE)
			*(system_maker.quantity(DISK_LOCK_FREQ)
			  -
			  system_maker.quantity(PRIMARY_ANGVEL_CORE))
			+
			2.0/3.0*system_maker.quantity(PRIMARY_CORE_MASS_DERIV)*
			std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2)*
			system_maker.quantity(DISK_LOCK_FREQ);
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
			*system_maker.quantity(PRIMARY_ANGVEL_ENV);
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
					 *angmom_env/Ienv;
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
											true, true, true);
		std::valarray<double> expected_diff_eq(0.0, 13);
		double orbit_power=0,
			   m1=system_maker.quantity(PRIMARY_MASS),
			   m2=system_maker.quantity(SECONDARY_MASS),
			   r1=system_maker.quantity(PRIMARY_RADIUS),
			   r2=system_maker.quantity(SECONDARY_RADIUS),
			   a=system_maker.quantity(SEMIMAJOR);
		for(unsigned zone_ind=0; zone_ind<4; ++zone_ind) {
			const Lags &lag=system_maker.lags(zone_ind);
			double power_norm, torque_norm;
			if(zone_ind<2) {
				power_norm=power_norm_Lai(m1, m2, r1, a);
				torque_norm=torque_norm_Lai(m2, r1, a);
			} else {
				power_norm=power_norm_Lai(m2, m1, r2, a);
				torque_norm=torque_norm_Lai(m1, r2, a);
			}
			orbit_power-=power_norm*dimensionless_power_Lai(0, lag);
			expected_diff_eq[9+zone_ind]=
				torque_norm*dimensionless_torque_z_Lai(0, lag);
		}
		double orbit_energy=m1*m2/(2.0*a)*AstroConst::G
							*AstroConst::solar_mass
							*std::pow(AstroConst::day, 2) 
							/std::pow(AstroConst::solar_radius,3),
			   a_rate=-a*orbit_power/orbit_energy,
			   w1=system_maker.quantity(PRIMARY_ANGVEL_ENV),
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
			primary_coupling_torque+=
				2.0/3.0*w2*system_maker.quantity(PRIMARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(PRIMARY_CORE_RADIUS), 2);
		if(system_maker.quantity(SECONDARY_CORE_MASS_DERIV)>0)
			secondary_coupling_torque-=
				2.0/3.0*w3*system_maker.quantity(SECONDARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(SECONDARY_CORE_RADIUS), 2);
		else
			secondary_coupling_torque+=
				2.0/3.0*w4*system_maker.quantity(SECONDARY_CORE_MASS_DERIV)*
				std::pow(system_maker.quantity(SECONDARY_CORE_RADIUS), 2);
		expected_diff_eq[0]=-13.0*std::pow(a, 7.5)*orbit_power/(m1*m2)
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
}

void test_BinarySystem::test_binary_locks_diff_eq()
{
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
	TEST_ADD(test_BinarySystem::test_single_zero_periapsis_diff_eq);*/
	TEST_ADD(
		test_BinarySystem::test_binary_no_locks_circular_aligned_diff_eq);
	TEST_ADD(
		test_BinarySystem::test_binary_no_locks_circular_inclined_diff_eq);
	TEST_ADD(test_BinarySystem::test_binary_locks_diff_eq);
}

#ifdef STANDALONE
int main()
{
	std::srand(std::time(NULL));
	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_BinarySystem tests(10000);
	return (tests.run(output, false) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
