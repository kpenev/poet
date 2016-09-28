#include "test_DissipatingBody.h"

TwoZoneBody *test_DissipatingBody::random_body(double &other_mass, double &a,
        Lags &lags_env, Lags &lags_core, bool no_periapsis,
        bool same_inclination) const
{
    other_mass=std::pow(10.0, uniform_rand(-2.0, 1.0));
    a=uniform_rand(3, 30);
    double m_env=std::pow(10.0, uniform_rand(-2.0, 1.0)),
           m_core=std::pow(10.0, uniform_rand(-5.0, 0.0))*m_env,
           r_env=std::pow(10.0, uniform_rand(-1, 1)),
           r_core=std::pow(10.0, uniform_rand(-1, 0))*r_env,
           orbit_freq=orbital_angular_velocity(m_env, other_mass, a),
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
    if(no_periapsis) periapsis_core=0;
    for(int m=-2; m<=2; ++m)
        for(int mp=-2; mp<=0; ++mp) {
            lags_env(m, mp)=(uniform_rand(0, 1)<0.2 ? 0
                                                    : uniform_rand(0, 10));
            lags_env(-m, -mp)=-lags_env(m, mp);
            lags_core(m, mp)=(uniform_rand(0, 1)<0.2 ? 0
                                                     : uniform_rand(0, 10));
            lags_core(-m, -mp)=-lags_core(m, mp);
        }
    ConstPhaseLagDissipatingZone 
        *envelope=new ConstPhaseLagDissipatingZone(lags_env, inertia_env,
                                                   r_env, m_env),
        *core=new ConstPhaseLagDissipatingZone(lags_core, inertia_core,
                                               r_core, m_core);
    TwoZoneBody *result=new TwoZoneBody(*envelope, *core, coupling_timescale,
                                        wind_strength, wind_sat_freq);
    std::valarray<double> angmom(2), inclination(2);
    angmom[0]=spin_freq_env*inertia_env;
    angmom[1]=spin_freq_core*inertia_core;
    inclination[0]=uniform_rand(0, M_PI);
    inclination[1]=(same_inclination ? inclination[0]
                                     : uniform_rand(0, M_PI));

    result->configure(age, other_mass, a, 0, &(angmom[0]), &(inclination[0]),
                      &periapsis_core, false, false, true);
    return result;
}

void test_DissipatingBody::test_Lai_torque_power()
{
    std::vector<std::string> quantity_names(4);
    quantity_names[0]="x-torque";
    quantity_names[1]="y-torque";
    quantity_names[2]="z-torque";
    quantity_names[3]="power";
    for(unsigned test_ind=0; test_ind<__ntests; ++test_ind) {
        double other_mass, a;
        Lags lags_env, lags_core;
        TwoZoneBody *body=random_body(other_mass, a, lags_env, lags_core);
        double power_norm=power_norm_Lai(body->mass(), other_mass,
                body->radius(), a),
               torque_norm=torque_norm_Lai(other_mass, body->radius(), a);
        for(unsigned zone_ind=0; zone_ind<2; ++zone_ind) {
            const Lags &lags=(zone_ind==0 ? lags_env : lags_core);
            double inclination=body->zone(zone_ind).inclination(),
                   expected_power=power_norm
                                  *
                                  dimensionless_power_Lai(inclination, lags),
                   power_above=body->tidal_power(zone_ind, true),
                   power_below=body->tidal_power(zone_ind, false);
            Eigen::Vector3d 
                torque_above=body->tidal_torque(zone_ind, true),
                torque_below=body->tidal_torque(zone_ind, false),
                expected_torque=torque_norm*Eigen::Vector3d(
                        dimensionless_torque_x_Lai(inclination, lags),
                        0,
                        dimensionless_torque_z_Lai(inclination, lags));
                for(int i=0; i<4; ++i) {
                    double expected=(i==3 ? expected_power
                                          : expected_torque(i)),
                           above=(i==3 ? power_above : torque_above(i)),
                           below=(i==3 ? power_below : torque_below(i));
                    std::ostringstream msg;
                    msg << "Expected " << (zone_ind==0 ? "envelope" : "core")
                        << " tidal " << quantity_names[i] 
                        << "=" << expected << " got (" << above << ", " 
                        << below << "), differences (" << above-expected 
                        << ", " << below-expected << ")";
                    TEST_ASSERT_MSG(check_diff(above, expected, 1e-10,1e-15),
                                    msg.str().c_str());
                    TEST_ASSERT_MSG(check_diff(below, expected, 1e-10,1e-15),
                                    msg.str().c_str());
                }
        }
        delete &(body->zone(0));
        delete &(body->zone(1));
        delete body;
    }
}

void test_DissipatingBody::test_orbit_rates_same_periapsis()
{
    for(unsigned test_ind=0; test_ind<__ntests; ++test_ind) {
        double other_mass, a;
        Lags lags_env, lags_core;
        TwoZoneBody *body=random_body(other_mass, a, lags_env, lags_core,
                                      true);
        double power_norm=power_norm_Lai(body->mass(), other_mass,
                body->radius(), a),
               torque_norm=torque_norm_Lai(other_mass, body->radius(), a);
        double energy_gain=0;
        Eigen::Vector3d orbit_torque_env(0, 0, 0),
                        orbit_torque_core(0, 0, 0);
        for(unsigned zone_ind=0; zone_ind<2; ++zone_ind) {
            const Lags &lags=(zone_ind==0 ? lags_env : lags_core);
            double this_inclination=body->zone(zone_ind).inclination(),
                   other_inclination=body->zone(1-zone_ind).inclination(),
                   sin_inc_diff=std::sin(this_inclination-other_inclination),
                   cos_inc_diff=std::cos(this_inclination-other_inclination);
            Eigen::Vector3d &this_orbit_torque(zone_ind==0
                                               ? orbit_torque_env
                                               : orbit_torque_core),
                            &other_orbit_torque(zone_ind==0
                                                ? orbit_torque_core
                                                : orbit_torque_env);
            energy_gain-=power_norm*dimensionless_power_Lai(this_inclination,
                                                            lags);
            Eigen::Vector3d 
                zone_torque(dimensionless_torque_x_Lai(this_inclination,
                                                       lags),
                            0,
                            dimensionless_torque_z_Lai(this_inclination,
                                                       lags));
            zone_torque*=torque_norm;
            this_orbit_torque-=zone_torque;
            other_orbit_torque[0]-=zone_torque[0]*cos_inc_diff
                                   -
                                   zone_torque[2]*sin_inc_diff;
            other_orbit_torque[2]-=zone_torque[0]*sin_inc_diff
                                   +
                                   zone_torque[2]*cos_inc_diff;
        }
        double value=body->tidal_orbit_energy_gain();
        std::ostringstream msg;
        msg << "Expected orbit energy gain: " 
            << energy_gain << " got: " << value << ", difference: " 
            << energy_gain-value;
        TEST_ASSERT_MSG(check_diff(value, energy_gain, 1e-10, 1e-15),
                        msg.str().c_str());
        for(unsigned zone_ind=0; zone_ind<2; ++zone_ind) {
            const Eigen::Vector3d &expected_torque=(zone_ind==0
                                                    ? orbit_torque_env
                                                    : orbit_torque_core);
            for(unsigned method=0; method<2-zone_ind; ++method) {
                const Eigen::Vector3d &got_torque=
                    (method==0
                     ? body->tidal_orbit_torque(body->zone(zone_ind))
                     :body->tidal_orbit_torque());
                for(int i=0; i<3; ++i) {
                    msg.str("");
                    msg << "Expected orbit torque(" << i << ") for zone"
                        << zone_ind << ": " << expected_torque(i) << ", got: "
                        << got_torque(i) << ", difference: "
                        << got_torque(i)-expected_torque(i);
                    TEST_ASSERT_MSG(check_diff(got_torque(i),
                                               expected_torque(i), 1e-10,
                                               1e-15),
                                    msg.str().c_str());
                }
            }
        }
    }
}

test_DissipatingBody::test_DissipatingBody(unsigned ntests,
            const std::string &eccentricity_expansion) : __ntests(ntests)
{
    DissipatingZone::read_eccentricity_expansion(eccentricity_expansion);
    TEST_ADD(test_DissipatingBody::test_Lai_torque_power);
    TEST_ADD(test_DissipatingBody::test_orbit_rates_two_zones);
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
