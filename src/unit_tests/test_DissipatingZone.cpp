#include "test_DissipatingZone.h"

double &Lags::operator()(int m, int mp)
{
	return (*this)[std::pair<int, int>(m, mp)];
}

double Lags::operator()(int m, int mp) const
{
	Lags::const_iterator it=this->find(std::pair<int, int>(m, mp));
	assert(it!=this->end());
	return it->second;
}

void ConstPhaseLagDissipatingZone::describe(std::ostream &os) const
{
	os << "I=" << inclination() << ", O(e)=" << eccentricity_order();
	for(Lags::const_iterator i=__lags.begin(); i!=__lags.end(); ++i) 
		os << ", D'_{" << i->first.first << "," << i->first.second << "}="
			<< i->second;
}

double test_DissipatingZone::torque_z_Lai(double inclination,
		const Lags &lags)
{
	double c=std::cos(inclination), s=std::sin(inclination);
	return 0.15*M_PI*(
			std::pow(1+c, 4)/2*lags(2, 2)
			+ std::pow(s*(1+c), 2)*lags(1, 2)
			- std::pow(s*(1-c), 2)*lags(-1, 2)
			- std::pow(1-c, 4)/2*lags(-2, 2))
		+ 3.0*M_PI/5.0*(std::pow(s, 4)/2*lags(2, 0)
				+ std::pow(s*c, 2)*lags(1, 0));
}

double test_DissipatingZone::torque_x_Lai(double inclination,
		const Lags &lags)
{
	double c=std::cos(inclination), s=std::sin(inclination);
	return 0.15*M_PI*(s*std::pow(1+c, 3)*lags(2,2)/2
			+ s*std::pow(1+c, 2)*(2.0-c)*lags(1,2)
			+ 3.0*std::pow(s, 3)*lags(0,2)
			+ s*std::pow(1-c, 2)*(2.0+c)*lags(-1, 2)
			+ s*std::pow(1-c, 3)*lags(-2, 2)/2)
		- 3.0*M_PI/5.0*(std::pow(s, 3)*c*lags(2,0)/2
				+ s*std::pow(c, 3)*lags(1,0));
}

double test_DissipatingZone::power_Lai(double inclination,
		const Lags &lags)
{
	double c=std::cos(inclination), s=std::sin(inclination);
	return 0.15*M_PI*(
			std::pow(1+c, 4)*lags(2,2)/2
			+ 2.0*std::pow(s*(1+c), 2)*lags(1,2)
			+ 3.0*std::pow(s, 4)*lags(0, 2)
			+ 2.0*std::pow(s*(1-c), 2)*lags(-1, 2)
			+ std::pow(1-c, 4)*lags(-2, 2)/2);
}

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
						torque_z_Lai(inclination, lags), 
						torque_x_Lai(inclination, lags), 
						power_Lai(inclination, lags));
				if(e_order==0) {
					single_test(zone, 1.0, 0.5, 1.0, 1e-3, 0.5, inclination,
						   		M_PI/3, torque_z_Lai(inclination, lags), 
								torque_x_Lai(inclination, lags), 
								power_Lai(inclination, lags));
					single_test(zone, 10.0, 0.5, 1.0, 1e-3, 0.5, inclination,
						   		M_PI/3, torque_z_Lai(inclination, lags), 
								torque_x_Lai(inclination, lags), 
								power_Lai(inclination, lags));
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
