#include "test_DissipatingZone.h"

void test_DissipatingZone::single_test(TestingDissipatingZone &zone,
	   	double orbital_frequency, double eccentricity, double m1, double m2,
		double spin_angmom, double inclination, double periapsis,
		double torque_z, double torque_x, double power)
{
	__current_failed=true;
	std::ostringstream msg_start;
	zone.describe(msg_start);
	msg_start << ", orbital freq=" << orbital_frequency
		<< ", e=" << eccentricity;
	std::ostringstream msg;
	zone.configure(1.0, orbital_frequency, eccentricity,
				   orbital_angmom_from_freq(m1, m2, orbital_frequency,
					   						eccentricity),
				   spin_angmom, inclination, periapsis);
	double above_value=zone.tidal_torque_z(true),
		   below_value=zone.tidal_torque_z(false);
	msg << msg_start.str()
		<< ", torque_z=(" << above_value << ", " << below_value
		<< "), expected=" << torque_z
		<< ", differences=(" << above_value-torque_z << ", "
		<< below_value-torque_z << ").";
	TEST_ASSERT_MSG(check_diff(above_value, torque_z, 1e-10, 1e-15), 
			msg.str().c_str());
	TEST_ASSERT_MSG(check_diff(below_value, torque_z, 1e-10, 1e-15), 
			msg.str().c_str());
	above_value=zone.tidal_torque_x(true);
	below_value=zone.tidal_torque_x(false);
	msg.str("");
	msg << msg_start.str()
		<< ", torque_x=(" << above_value << ", " << below_value
		<< "), expected=" << torque_x
		<< ", difference=(" << above_value-torque_x << ", "
		<< below_value-torque_x << ").";
	TEST_ASSERT_MSG(check_diff(above_value, torque_x, 1e-10, 1e-15), 
			msg.str().c_str());
	TEST_ASSERT_MSG(check_diff(below_value, torque_x, 1e-10, 1e-15), 
			msg.str().c_str());
	above_value=zone.tidal_power(true);
	below_value=zone.tidal_power(false);
	msg.str("");
	msg << msg_start.str()
		<< ", power=(" << above_value << ", " << below_value
		<< "), expected=" << power
		<< ", difference=(" << above_value-power << ", "
		<< below_value-power << ").";
	TEST_ASSERT_MSG(check_diff(above_value, power, 1e-10, 1e-15), 
			msg.str().c_str());
	TEST_ASSERT_MSG(check_diff(below_value, power, 1e-10, 1e-15), 
			msg.str().c_str());
	__current_failed=false;
}

test_DissipatingZone::test_DissipatingZone(
		const std::string &eccentricity_expansion)
{
	DissipatingZone::read_eccentricity_expansion(eccentricity_expansion);
	TEST_ADD(test_DissipatingZone::test_Lai);
	TEST_ADD(test_DissipatingZone::test_inclination_periapsis_evol_deriv);
}

///Tests the torque expression agains Eq. (27) in Lai 2012.
void test_DissipatingZone::test_Lai()
{
	std::cerr << std::endl << std::endl;
	Lags lags;
	for(int e_order=0; e_order<3; ++e_order) {
		for(double inclination=0; inclination<3*M_PI/4; inclination+=M_PI/8){
			for(int i=0; i<std::pow(2, 15); ++i) {
				int enabled_terms_mask=i;
				for(int m=-2; m<=2; ++m) {
					for(int mp=-2-e_order; mp<=0; ++mp) {
						lags(m, mp)=enabled_terms_mask%2;
						lags(-m, -mp)=-enabled_terms_mask%2;
						enabled_terms_mask/=2;
					}
				}
				ConstPhaseLagDissipatingZone zone(lags, 1.0, 1.0, 1.0);
				zone.change_e_order(e_order);
				if(i%1000==0) 
					std::cerr << "Lai test " << i << "\r";
				single_test(zone, 1.0, 0, 1.0, 1e-3, 1.0, inclination, 0.0, 
						dimensionless_torque_z_Lai(inclination, lags), 
						dimensionless_torque_x_Lai(inclination, lags), 
						dimensionless_power_Lai(inclination, lags));
				if(e_order==0) {
					single_test(zone, 1.0, 0.5, 1.0, 1e-3, 0.5, inclination,
						   		M_PI/3,
								dimensionless_torque_z_Lai(inclination,lags),
								dimensionless_torque_x_Lai(inclination,lags),
								dimensionless_power_Lai(inclination, lags));
					single_test(zone, 10.0, 0.5, 1.0, 1e-3, 0.5, inclination,
						   		M_PI/3,
								dimensionless_torque_z_Lai(inclination,lags),
								dimensionless_torque_x_Lai(inclination,lags),
								dimensionless_power_Lai(inclination, lags));
				}
				TEST_ASSERT_MSG(!__current_failed, "Test failure.");
			}
		}
	}
}

void test_DissipatingZone::test_inclination_periapsis_evol_deriv()
{
}

#ifdef STANDALONE
int main()
{
	std::cout.setf(std::ios_base::scientific);
	std::cout.precision(16);
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_DissipatingZone tests;
	return (tests.run(output, false) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
