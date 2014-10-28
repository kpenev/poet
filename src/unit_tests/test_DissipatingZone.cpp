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
	os << "I=" << __inclination << ", O(e)=" << __e_order;
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
	   	double orbital_frequency, double eccentricity, double torque_z,
		double torque_x, double power)
{
	__current_failed=true;
	std::ostringstream msg_start;
	zone.describe(msg_start);
	msg_start << ", orbital freq=" << orbital_frequency
		<< ", e=" << eccentricity;
	std::ostringstream msg;
	zone.set_orbit(orbital_frequency, eccentricity);
	double value=zone.tidal_torque_z();
	msg << msg_start.str()
		<< ", torque_z=" << value
		<< ", expected=" << torque_z
		<< ", difference=" << value-torque_z;
	TEST_ASSERT_MSG(check_diff(value, torque_z, 1e-10, 1e-15), 
			msg.str().c_str());
	value=zone.tidal_torque_x();
	msg.str("");
	msg << msg_start.str()
		<< ", torque_x=" << value
		<< ", expected=" << torque_x
		<< ", difference=" << value-torque_x;
	TEST_ASSERT_MSG(check_diff(value, torque_x, 1e-10, 1e-15), 
			msg.str().c_str());
	value=zone.tidal_power();
	msg.str("");
	msg << msg_start.str()
		<< ", power=" << value
		<< ", expected=" << power
		<< ", difference=" << value-power;
	TEST_ASSERT_MSG(check_diff(value, power, 1e-10, 1e-15), 
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
				ConstPhaseLagDissipatingZone zone(inclination, e_order,
						lags);
				if(i%1000==0) 
					std::cerr << "Lai test " << i << "\r";
				single_test(zone, 0, 0, torque_z_Lai(inclination, lags), 
						torque_x_Lai(inclination, lags), 
						power_Lai(inclination, lags));
				if(e_order==0) {
					single_test(zone, 0, 0.5, torque_z_Lai(inclination, lags), 
							torque_x_Lai(inclination, lags), 
							power_Lai(inclination, lags));
					single_test(zone, 1.0, 0.5, torque_z_Lai(inclination, lags), 
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
