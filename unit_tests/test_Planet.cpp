/**\file
 *
 * \brief The entire declaration/definition of the unit tests that exercise
 * the Planet class.
 *
 * Much cleanup is necessary.
 *
 * \ingroup UnitTests_group
 */

#include "test_Star.h"
#include "Star.h"
#include "StellarSystem.h"
#include "Planet.h"
#include "Common.h"
#include <assert.h>
#include <limits>

class test_Planet : public Test::Suite {
	Star* star;
	PlanetData d;
	Planet* planet;
	double start_age, init_rad_L, init_conv_L, curr_age;
	std::valarray<double> ages, rad_Ls, conv_Ls;
public:
  test_Planet(): d(), start_age(1), init_rad_L(3), init_conv_L(6) {
	  using namespace std;
    d.init_random_planet(&planet);
    star = d.star;
    cout << "finished init"<<endl;
    TEST_ADD(test_Planet::test_basic);
    //TEST_ADD(test_Planet::test_known);
    TEST_ADD(test_Planet::test_angmom);
    TEST_ADD(test_Planet::test_evolution);
    //TEST_ADD(test_Planet::test_tidal);
  }
  void test_basic() {
    using namespace std;
    cout << "Mass: " << planet->get_mass() << endl;
    cout << "Lifetime: " << planet->get_lifetime() << endl;
    cout << "Current semimajor: " <<
      planet->get_current_semimajor() << endl;
    cout << "Radius: " << planet->get_radius() << endl;
  }
  
  void test_known() {
	  //this test depends on the exact data passed to Star and Planet in the
	  //constructor
	  using namespace std;
	  using namespace AstroConst;
	  double step = 0.02;
	  double pmass = planet->get_mass()*jupiter_mass;
	  double smass = star->get_mass()*solar_mass;
	  for (double age=0.1; age < 1; age += step) {
		  double semimajor = planet->get_semimajor(age);
		  double semiDeriv = planet->get_semimajor_derivative(age);
		  assert(approxEqual(semimajor, 0.1*(1-age*age)));
		  assert(approxEqual(semiDeriv, -0.2*age));
		  double a = planet->get_semimajor(age)*AU;
		  double physSemiDeriv = semiDeriv*AU/Gyr;
		  double angmomDeriv = 0.5*pmass*smass*sqrt(G/a/(pmass+smass))*physSemiDeriv;
		  angmomDeriv *= (day*Gyr/solar_mass/solar_radius/solar_radius);
		  assert(approxEqual(planet->orbital_angmom_deriv(age), angmomDeriv));
		  assert(approxEqual(
				  planet->orbital_angmom_deriv(semimajor, semiDeriv),
				  angmomDeriv));
		//  double assumedStellarFreq = 0.17;
	  }
  }

  void test_evolution() {
    using namespace std;
    double step = 0.02;
    for (double age=0.1; age < planet->get_lifetime(); age+=step) {
      double testSemiDeriv = 
	(planet->get_semimajor(age+step/2) -
	 planet->get_semimajor(age-step/2))/step;
      double semiDeriv = planet->get_semimajor_derivative(age);
      TEST_ASSERT(approxEqual(testSemiDeriv, semiDeriv));
      TEST_ASSERT(isEqual(semiDeriv, planet->tidal_decay(age)));
      double minSemi = planet->minimum_semimajor(age);
      TEST_ASSERT(minSemi >= 0);
      double angmomDeriv = planet->orbital_angmom_deriv(age);
      if (angmomDeriv > 0)
    	  cout << "Warning: angular momentum is increasing" << endl;
      double orbPeriod = planet->orbital_period_age(age);
      double semi = planet->get_semimajor(age);
      double testOrbPeriod = planet->orbital_period_semimajor(semi);
      TEST_ASSERT(isEqual(testOrbPeriod, orbPeriod));
      
    }
  }

  double try_calc(bool use_a6p5, double semimajor) {
	  double rstar = 4.173e8;
	  double mplanet = 3.7974e28;
	  double Q = 127453;
	  double mstar = 1.98892e30;

	if(use_a6p5)
		return -4.5*6.5*std::sqrt(AstroConst::G/mstar)*
			std::pow(rstar, 5.0)*
			mplanet/Q*AstroConst::Gyr/
			std::pow(AstroConst::AU, 6.5);
	else {
		double a=semimajor*AstroConst::AU;
		return -4.5*std::sqrt(AstroConst::G/(a*mstar))*
			std::pow(rstar/a, 5.0)*
			mplanet/Q*AstroConst::Gyr/AstroConst::AU;
	}
  }

  void test_tidal() {
	  using namespace std;
	  double ageStep = 0.05;
	  double semiStep = 5e-4;
	  double freqStep = 1e-3;
	  for (double age=0.1; age < 2; age+=ageStep) {
		  cout << age << endl;
		  for (double semi=1e-2; semi < 0.1; semi+=semiStep) {
			  for (double freq=1e-2; freq < 0.1; freq+=freqStep) {
				  double testAgeDeriv =
						  (planet->tidal_decay(age+ageStep/2, semi, freq, false)
								  -planet->tidal_decay(age-ageStep/2, semi, freq, false))
								  /ageStep;
				  double ageDeriv =
						  planet->tidal_decay_age_deriv(age, semi, freq, false);
				  assert(approxEqual(testAgeDeriv, ageDeriv));
	  
				  double testAgeDeriv6p5 =
						  (planet->tidal_decay(age+ageStep/2, semi, freq, true)
								  -planet->tidal_decay(age-ageStep/2, semi, freq, true))
								  /ageStep;
				  double ageDeriv6p5 = planet->tidal_decay_age_deriv(age, semi, freq, true);
				  assert(approxEqual(testAgeDeriv6p5,
						  6.5*pow(semi, 5.5)*testAgeDeriv));
				  assert(approxEqual(ageDeriv6p5, 6.5*pow(semi, 5.5)*ageDeriv));

				  double testSemiDeriv =
						  (planet->tidal_decay(age,semi+semiStep/2,freq,false)
								  -planet->tidal_decay(age,semi-semiStep/2,freq,false))
								  /semiStep;
				  double semiDeriv =
						  planet->tidal_decay_semimajor_deriv(age, semi, freq, false);
				  assert(approxEqual(testSemiDeriv, semiDeriv));
				  double testSemiDeriv6p5 =
						  (planet->tidal_decay(age,semi+semiStep/2,freq,true)
								  -planet->tidal_decay(age,semi-semiStep/2,freq,true))
								  /semiStep;
				  double semiDeriv6p5 =
						  planet->tidal_decay_semimajor_deriv(age, semi, freq, true);
		//		  cout << "f " << testSemiDeriv6p5 << " " << semiDeriv6p5 << endl;
				  if (!isEqual(testSemiDeriv6p5, semiDeriv6p5))
					  assert(approxEqual(testSemiDeriv6p5, semiDeriv6p5));

				  double testFreqDeriv =
						  (planet->tidal_decay(age, semi, freq+freqStep/2, false)
								  -planet->tidal_decay(age, semi, freq-freqStep/2, false))
								  /freqStep;
				  double freqDeriv =
						  planet->tidal_decay_star_spin_deriv(age, semi, freq, false);

				  if (!isEqual(testFreqDeriv, freqDeriv)) //in case both are close to 0
					  assert(approxEqual(testFreqDeriv, freqDeriv));
				  double testFreqDeriv6p5 =
						  (planet->tidal_decay(age, semi, freq+freqStep/2, true)
								  -planet->tidal_decay(age, semi, freq-freqStep/2, true))
								  /freqStep;
				  double freqDeriv6p5 =
						  planet->tidal_decay_star_spin_deriv(age, semi, freq, true);

	//			  cout << "d"<<testFreqDeriv6p5 << " "<< 6.5*pow(semi, 5.5)*testFreqDeriv << endl;
				  assert(approxEqual(testFreqDeriv6p5, 6.5*pow(semi, 5.5)*testFreqDeriv));
				  assert(isEqual(freqDeriv6p5, 6.5*pow(semi, 5.5)*freqDeriv));
		}
      }
    }
  }

  void test_angmom() {
	  using namespace std;
	  double semiStep=1e-4, derivStep=1e-3;
	  for (double semi=1e-2; semi < 0.1; semi+=semiStep) {
		  cout << semi << endl;
		  for (double deriv=-0.1; deriv > -1; deriv -= derivStep) {
			  double testAngmomDerivDeriv =
					  (planet->orbital_angmom_deriv(semi+semiStep/2, deriv)
							  -planet->orbital_angmom_deriv(semi-semiStep/2, deriv))
							  /semiStep;
			  double angmomDerivDeriv =
					  planet->orbital_angmom_deriv_semimajor_deriv(semi, deriv);
			  TEST_ASSERT(approxEqual(testAngmomDerivDeriv, angmomDerivDeriv));
      }
    }
  }
};

#ifdef STANDALONE
int main()
{
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_Planet tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);
	
}
#endif

