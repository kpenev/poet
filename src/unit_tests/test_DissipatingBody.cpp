#include "test_DissipatingBody.h"

TwoZoneBody *test_DissipatingBody::random_body()
{
	double m_env=std::pow(10.0, uniform_rand(-2.0, 1.0)),
		   m_core=std::pow(10.0, uniform_rand(-5.0, 0.0))*m_env,
		   m_other=std::pow(10.0, uniform_rand(-2.0, 1.0)),
		   r_env=std::pow(10.0, uniform_rand(-1, 1)),
		   r_core=std::pow(10.0, uniform_rand(-1, 0))*r_env,
		   a=uniform_rand(3, 30),
		   orbit_freq=orbital_angular_velocity(m_env, m_other, a),
		   periapsis_core=uniform_rand(0, 2.0*M_PI),
		   inertia_env=std::pow(10.0, uniform_rand(-5, 0))
			   		   *m_env*std::pow(r_env, 2),
		   inertia_core=std::pow(10.0, uniform_rand(-5, 0))
			   			*m_core*std::pow(r_core, 2),
		   spin_freq_env=std::pow(10.0, uniform_rand(-1, 1))*orbit_freq,
		   spin_freq_core=std::pow(10.0, uniform_rand(-1, 1))*orbit_freq,
		   coupling_timescale=uniform_rand(0, 100),
		   wind_strength=uniform_rand(0, 100),
		   wind_sat_freq=std::pow(10.0, uniform_rand(-2, 2))*orbit_freq,
		   age=uniform_rand(0, 100);
	Lags lags_env, lags_core;
	for(int m=-2; m<=2; ++m)
		for(int mp=-2; mp<=0; ++mp) {
			lags_env(m, mp)=(uniform_rand(0, 1)<0.2 ? 0
													: uniform_rand(0, 10));
			lags_env(-m, -mp)=lags_env(m, mp);
			lags_core(m, mp)=(uniform_rand(0, 1)<0.2 ? 0
													 : uniform_rand(0, 10));
			lags_core(-m, -mp)=lags_core(m, mp);
		}
	ConstPhaseLagDissipatingZone envelope(lags_env, inertia_env, r_env,
										  m_env),
								 core(lags_core, inertia_core, r_core,
									  m_core);
	TwoZoneBody *result=new TwoZoneBody(envelope, core, coupling_timescale,
										wind_strength, wind_sat_freq);
	std::valarray<double> angmom(2), inclination(2);
	angmom[0]=spin_freq_env*inertia_env;
	angmom[1]=spin_freq_core*inertia_core;
	inclination[0]=uniform_rand(0, M_PI);
	inclination[1]=uniform_rand(0, M_PI);

	result->configure(age, m_other, a, 0, &(angmom[0]), &(inclination[0]),
					  &periapsis_core, false, false, true);
	return result;
}

test_DissipatingBody::test_DissipatingBody(unsigned ntests,
			const std::string &eccentricity_expansion) : __ntests(ntests)
{
	DissipatingZone::read_eccentricity_expansion(eccentricity_expansion);
	TEST_ADD(test_DissipatingBody::test_Lai_torque_power);
}

void test_DissipatingBody::test_Lai_torque_power()
{
	for(unsigned test_ind=0; test_ind<__ntests; ++test_ind) {
		TwoZoneBody *body=random_body();
		delete body;
	}
}

void test_DissipatingBody::test_orbit_rates_two_zones()
{
}

#ifdef STANDALONE
int main()
{
	std::srand(std::time(NULL));
	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_DissipatingBody tests(10000);
	return (tests.run(output, false) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
