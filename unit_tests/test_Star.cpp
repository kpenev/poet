/**\file
 *
 * \brief Definitions of some of the methods  exercising the Star class.
 *
 * Lots of cleanup is necessary.
 *
 * \ingroup UnitTests_group
 */

#include "test_Star.h"
#include <assert.h>

test_Star::test_Star() : d()
{
	using namespace std;
	d.init_random_star(&star);
	TEST_ADD(test_Star::test_basic);
	TEST_ADD(test_Star::test_evolution);
	TEST_ADD(test_Star::test_tidal);
	TEST_ADD(test_Star::test_ang_momentum);
	TEST_ADD(test_Star::test_spin);
	TEST_ADD(test_Star::test_wind_torque);
	TEST_ADD(test_Star::test_diff_rot);

}

void test_Star::test_basic() {
	using namespace std;
	cout<<"Age: "<<star->current_age()<<endl;
	cout<<"Lrad: "<<star->get_current_angular_momentum(radiative)<<endl;
	cout<<"Lconv: "<<star->get_current_angular_momentum(convective)<<endl;
	cout<<"Lifetime: "<<star->get_lifetime()<<endl;
	cout<<"Mass: "<<star->get_mass()<<endl;
	cout<<"Name: "<<star->get_name()<<endl;
	cout<<"Radius: "<<star->get_radius()<<endl;
}

void test_Star::test_known() {

}

void test_Star::test_ang_momentum() {
	using namespace std;
	for (double age=d.evolution_ages[0] + 0.2; age < 10; age+=0.5) {
		double rad_L = star->get_angular_momentum(age, radiative);
		double conv_L = star->get_angular_momentum(age, convective);
		double tot_L = star->get_angular_momentum(age, total);
		TEST_ASSERT(isEqual(rad_L+conv_L, tot_L));
		//test for expected ang momenta
	}

}

void test_Star::test_tidal() {
	double step = 1e-6;
	for (double freq=1e-3; freq < 0.5; freq+=step) {
		double testDeriv = star->get_tidal_Q(freq + step/2)-
			star->get_tidal_Q(freq - step/2);
		testDeriv /= step;
		double deriv = star->get_tidal_Q_deriv(freq);
		if (deriv < 1e-20) break;
		TEST_ASSERT(approxEqual(testDeriv, deriv));
	}
}

void test_Star::test_evolution() {
	using namespace std;
	double start = 1;
	double step = 0.01;
	for (double age=start; age < 10; age+=step) {
		//core inertia test
		double testDeriv = star->core_inertia_gain(age + step/2) -
			star->core_inertia_gain(age - step/2);
		testDeriv /= step;
		double deriv = star->core_inertia_gain_deriv(age);
		TEST_ASSERT(approxEqual(testDeriv, deriv));
		//moment of inertia test
		double testRadIDeriv = star->moment_of_inertia(age+step/2, radiative) -
			star->moment_of_inertia(age-step/2, radiative);
		testRadIDeriv /= step;
		double radIDeriv = star->moment_of_inertia_deriv(age, radiative);
		TEST_ASSERT(approxEqual(testRadIDeriv, radIDeriv));

		double testConvIDeriv = star->moment_of_inertia(age+step/2, convective) -
			star->moment_of_inertia(age-step/2, convective);
		testConvIDeriv /= step;
		double convIDeriv = star->moment_of_inertia_deriv(age, convective);
		TEST_ASSERT(approxEqual(testConvIDeriv, convIDeriv));

		double totIDeriv = star->moment_of_inertia_deriv(age, total);
		TEST_ASSERT(isEqual(radIDeriv + convIDeriv, totIDeriv));
		double radI = star->moment_of_inertia(age, radiative);
		double convI = star->moment_of_inertia(age, convective);
		double totI = star->moment_of_inertia(age, total);
		TEST_ASSERT(isEqual(radI+convI, totI));

		//ang momentum test
		double radL = star->get_angular_momentum(age, radiative);
		double convL = star->get_angular_momentum(age, convective);
		double totL = star->get_angular_momentum(age, total);
		double radSpin = star->spin_frequency(age, radiative);
		double convSpin = star->spin_frequency(age, convective);
		TEST_ASSERT(isEqual(radSpin*radI, radL));
		TEST_ASSERT(isEqual(convSpin*convI, convL));
		TEST_ASSERT(isEqual(radL+convL, totL));
		double radSpinFromL = star->spin_frequency(age, radiative, radL);
		double convSpinFromL = star->spin_frequency(age, convective, convL);
		TEST_ASSERT(isEqual(radSpinFromL, radSpin));
		TEST_ASSERT(isEqual(convSpinFromL, convSpin));
		//diff. rotation test
		double diffRot = star->differential_rotation(age);
		double testDiffRot = star->differential_rotation(age, convL, radL);
		TEST_ASSERT(isEqual(diffRot, testDiffRot));
		//wind torque test
		double windTorque = star->wind_torque(age);
		double testWindTorque = star->wind_torque(age, convSpin);
		TEST_ASSERT(isEqual(windTorque, testWindTorque));
		//test radius
		double testAgeDeriv = star->get_radius(age + step/2) -
			star->get_radius(age - step/2);
		testAgeDeriv /= step;
		double ageDeriv = star->get_logradius_deriv(age);
		ageDeriv *= star->get_radius(age);
		TEST_ASSERT(approxEqual(testAgeDeriv, ageDeriv));
		TEST_ASSERT(isEqual(star->differential_rotation(age),
					star->differential_rotation(age, convL, radL)));
		double radMass = star->get_zone_mass(age, radiative);
		double convMass = star->get_zone_mass(age, convective);
		double totMass = star->get_zone_mass(age, total);
		TEST_ASSERT(radMass > 0 && convMass > 0 && totMass > 0);
		TEST_ASSERT(isEqual(radMass + convMass, totMass));
	}
}

void test_Star::test_spin() {
	using namespace std;
	double step = 0.1;
	for (double age = 1; age < 10; age+=step) {
		for (double L = 0.1; L < 5; L += 0.5) {
			double radSpin = star->spin_frequency(age, radiative, L);
			double convSpin = star->spin_frequency(age, convective, L);
			double radI = star->moment_of_inertia(age, radiative);
			double convI = star->moment_of_inertia(age, convective);
			TEST_ASSERT(isEqual(radSpin*radI, L));
			TEST_ASSERT(isEqual(convSpin*convI, L));

			double testRadSpinDeriv =
				star->spin_frequency(age+step/2, radiative, L) -
				star->spin_frequency(age-step/2, radiative, L);
			testRadSpinDeriv /= step;
			double radSpinDeriv =
				star->spin_frequency_age_deriv(age, radiative, L);
			TEST_ASSERT(approxEqual(testRadSpinDeriv, radSpinDeriv));

			double testConvSpinDeriv =
				star->spin_frequency(age+step/2, convective, L) -
				star->spin_frequency(age-step/2, convective, L);
			testConvSpinDeriv /= step;
			double convSpinDeriv =
				star->spin_frequency_age_deriv(age, convective, L);
			TEST_ASSERT(approxEqual(testConvSpinDeriv, convSpinDeriv));
		}
	}
}

void test_Star::test_wind_torque() {
	double ageStep = 0.05;
	double freqStep = 1e-5;
	for (double age=0.5; age < 10; age += ageStep) {
		std::cout << "at age " << age << std::endl;
		for (double freq=1e-4; freq < 0.1; freq += freqStep) {
			double testAgeDeriv = star->wind_torque(age+ageStep/2, freq) -
				star->wind_torque(age-ageStep/2, freq);
			testAgeDeriv /= ageStep;
			double ageDeriv = star->wind_torque_age_deriv(age, freq);
			//TEST_ASSERT(approxEqual(testAgeDeriv, ageDeriv));
			std::cout<<testAgeDeriv<<" "<<ageDeriv<<std::endl;
			assert(approxEqual(testAgeDeriv, ageDeriv));
			double testFreqDeriv = star->wind_torque(age, freq+freqStep/2) -
				star->wind_torque(age, freq-freqStep/2);
			testFreqDeriv /= freqStep;
			double freqDeriv = star->wind_torque_freq_deriv(age, freq);
			//TEST_ASSERT(approxEqual(testFreqDeriv, freqDeriv));
			assert(approxEqual(testFreqDeriv, freqDeriv));
		}
		double convSpin = star->spin_frequency(age, convective);
		//TEST_ASSERT(isEqual(star->wind_torque(age),
			//		star->wind_torque(age, convSpin)));
		assert(isEqual(star->wind_torque(age),
							star->wind_torque(age, convSpin)));
	}
}

void test_Star::test_diff_rot() {
	using namespace std;
	double ageStep = 0.1;
	double radLStep = 0.02;
	double convLStep = 0.02;
	for (double age=0.5; age < 10; age+=ageStep) {
		double convL = star->get_angular_momentum(age, convective);
		double radL = star->get_angular_momentum(age, radiative);
		TEST_ASSERT(isEqual(star->differential_rotation(age),
					star->differential_rotation(age, convL, radL)));
		double diffRot = star->differential_rotation(age);
		double convSpin = star->spin_frequency(age, convective);
		TEST_ASSERT(isEqual(
					star->differential_rotation_torque(age, diffRot, convSpin),
					star->differential_rotation_torque(age)));
		for (double convL = 0.2; convL < 5; convL += convLStep) {
			for (double radL = 0.1; radL < 4; radL += radLStep) {
				if (radL >= convL) continue;
				double testAgeDeriv =
					star->differential_rotation(age+ageStep/2, convL, radL)-
					star->differential_rotation(age-ageStep/2, convL, radL);
				testAgeDeriv /= ageStep;
				double ageDeriv = star->differential_rotation_deriv(age, convL, radL);
				TEST_ASSERT(approxEqual(testAgeDeriv, ageDeriv));

				double testConvDeriv =
					star->differential_rotation(age,convL+convLStep/2, radL)-
					star->differential_rotation(age, convL-convLStep/2, radL);
				testConvDeriv /= convLStep;
				double convDeriv =
					star->differential_rotation_deriv(
							age, convL, radL, convective);
				TEST_ASSERT(approxEqual(testConvDeriv, convDeriv));

				double testRadDeriv =
					star->differential_rotation(age,convL, radL+radLStep/2)-
					star->differential_rotation(age, convL, radL-radLStep/2);
				testRadDeriv /= radLStep;
				double radDeriv =
					star->differential_rotation_deriv(
							age, convL, radL, radiative);
				TEST_ASSERT(approxEqual(testRadDeriv, radDeriv));
			}
		}

	}
}

#ifdef STANDALONE
int main()
{
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_Star tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif

