/**\file
 * 
 * \brief Defines some of the methods of the test suite that exercised the
 * EvolvingStellarQuantity class.
 *
 * \ingroup UnitTests_group
 */

#include "test_EvolvingStellarQuantity.h"

const OneArgumentDiffFunction *test_EvolvingStellarQuantity::closest_track(
		double m, const std::valarray<double> &track_masses,
		const std::list<const OneArgumentDiffFunction *> &tracks)
{
	if(m<track_masses[0]) return tracks.front();
	std::list<const OneArgumentDiffFunction *>::const_iterator 
		track_iter=tracks.begin();
	for(size_t i=1; i<track_masses.size(); i++) {
		if(track_masses[i]>m) {
			if(std::abs(track_masses[i]-m)<std::abs(track_masses[i-1]-m)) {
				return *(++track_iter);
			} else return *track_iter;
		}
		track_iter++;
	}
	return tracks.back();

}

void test_EvolvingStellarQuantity::test_low_mass_no_death_interp(bool log_age)
{
	try {
		double dm=0.01;
		int num_tracks=static_cast<int>((max_low_mass-min_stellar_mass)/dm)+1;
		std::valarray<double> track_masses(num_tracks);
		std::valarray< std::valarray<double> > poly_coef(std::valarray<double>(3), 3);
		poly_coef[0][0]=0.1; poly_coef[0][1]=0.2; poly_coef[0][2]=0.3;
		poly_coef[1][0]=0.2; poly_coef[1][1]=0.3; poly_coef[1][2]=0.4;
		poly_coef[2][0]=0.3; poly_coef[2][1]=0.4; poly_coef[2][2]=0.5;
		std::list<const OneArgumentDiffFunction *> tracks;
		int mass_index=0; 
		for(double m=min_stellar_mass; m<max_low_mass; m+=dm) {
			track_masses[mass_index]=m;
			tracks.push_back(exact_track(poly_coef, m));
			mass_index++;
		}
		for(unsigned rand_mass_ind=0; rand_mass_ind<100; rand_mass_ind++) {
			double rand_mass=rand_value(min_stellar_mass, max_low_mass);
			EvolvingStellarQuantity q(rand_mass, track_masses, tracks,
					log_age, max_low_mass, 0.0, 0.0);
			PolynomialEvolutionTrack *exact=exact_track(poly_coef, rand_mass);
			for(unsigned rand_age_ind=0; rand_age_ind<100; rand_age_ind++) {
				double rand_age=rand_value(min_age, max_age),
					   expected=(*exact)(rand_age);
				if(log_age) rand_age=std::exp(rand_age);
				double got=q(rand_age);
				std::ostringstream test_msg;
				test_msg << "Evaluating interpolation over polynomial tracks"
				   "(rand_mass #" << rand_mass_ind << " = " << rand_mass
				   << ", rand_age #" << rand_age_ind << " = " << rand_age
				   << "): expected " << expected << ", got " << got;
				TEST_ASSERT_MSG(check_diff(expected, got, 1e-12, 1e-12),
						test_msg.str().c_str()) 
				for(unsigned deriv_order=0; deriv_order<poly_coef.size(); 
						deriv_order++) {
					const FunctionDerivatives *q_deriv=q.deriv(rand_age);
					expected=exact->order(deriv_order);
					got=q_deriv->order(deriv_order);
					if(log_age && deriv_order>=1) {
						expected/=rand_age;
						if(deriv_order==2)
							expected-=exact->order(1)/(rand_age*rand_age);
					}
					test_msg.str("");
					test_msg << "Evaluating interpolation derivative "
						"order " << deriv_order << " over polynomial "
						"tracks (rand_mass #" << rand_mass_ind << " = "
						<< rand_mass << ", rand_age #" << rand_age_ind 
						<< " = " << rand_age << "): expected " 
						<< expected << ", got " << got;
					TEST_ASSERT_MSG(check_diff(expected, got, 1e-12, 1e-12), 
							test_msg.str().c_str()) 
				}
			}
		}
	} catch(Error::General &ex) {
		TEST_ASSERT_MSG(false, (std::string("test_low_mass_no_death_interp: "
						"unexpected exception thrown:")+ex.what()+":"+
					ex.get_message()).c_str())
	} catch(std::exception &ex) {
		TEST_ASSERT_MSG(false, 
			(std::string("test_low_mass_no_death_interp: unexpected exception"
						 "thrown:")+ex.what()).c_str())
	} catch(...) {
		TEST_ASSERT_MSG(false, "test_low_mass_no_death_interp: unexpected "
				"unknown exception thrown:")
	}
}

void test_EvolvingStellarQuantity::test_high_mass_no_death_interp(bool log_age)
{
	try {
		double dm=0.01;
		int num_tracks=static_cast<int>((max_stellar_mass-max_low_mass)/dm)+1;
		std::valarray<double> track_masses(num_tracks);
		std::valarray< std::valarray<double> > poly_coef(std::valarray<double>(3), 3);
		poly_coef[0][0]=0.1; poly_coef[0][1]=0.2; poly_coef[0][2]=0.3;
		poly_coef[1][0]=0.2; poly_coef[1][1]=0.3; poly_coef[1][2]=0.4;
		poly_coef[2][0]=0.3; poly_coef[2][1]=0.4; poly_coef[2][2]=0.5;
		std::list<const OneArgumentDiffFunction *> tracks;
		int mass_index=0; 
		for(double m=max_low_mass; m<max_stellar_mass; m+=dm) {
			track_masses[mass_index]=m;
			tracks.push_back(exact_track(poly_coef, m));
			mass_index++;
		}
		for(unsigned rand_mass_ind=0; rand_mass_ind<100; rand_mass_ind++) {
			double rand_mass=rand_value(max_low_mass, max_stellar_mass);
			EvolvingStellarQuantity q(rand_mass, track_masses, tracks, log_age,
					max_low_mass, 0.0, 0.0);
			PolynomialEvolutionTrack *exact=exact_track(poly_coef, rand_mass);
			const OneArgumentDiffFunction *closest=closest_track(rand_mass, 
					track_masses, tracks);
			for(unsigned rand_age_ind=0; rand_age_ind<100; rand_age_ind++) {
				double rand_age=rand_value(min_age, max_age),
					   expected=(*exact)(rand_age),
					   max_error=std::abs((*closest)(rand_age)-expected),
					   exp_rand_age=std::exp(rand_age),
					   got=q((log_age ? exp_rand_age : rand_age));
				std::ostringstream test_msg;
				test_msg << "Evaluating interpolation over polynomial tracks"
				   "(rand_mass #" << rand_mass_ind << " = " << rand_mass
				   << ", rand_age #" << rand_age_ind << " = " << exp_rand_age
				   << "): expected " << expected << "+-" << max_error 
				   << ", got " << got;
				TEST_ASSERT_MSG(check_diff(expected, got, 1e-12, max_error),
						test_msg.str().c_str()) 
				for(unsigned deriv_order=0; deriv_order<poly_coef.size(); 
						deriv_order++) {
					const FunctionDerivatives *q_deriv=q.deriv(
							(log_age ? exp_rand_age : rand_age));
					expected=exact->order(deriv_order);
					if(log_age && deriv_order>=1) {
						expected/=exp_rand_age;
						if(deriv_order==2)
							expected-=exact->order(1)/(exp_rand_age*exp_rand_age);
					}
					got=q_deriv->order(deriv_order);
					max_error=std::abs(
							closest->deriv(rand_age)->order(deriv_order)-expected);
					test_msg.str("");
					test_msg << "Evaluating interpolation derivative "
						"order " << deriv_order << " over polynomial "
						"tracks (rand_mass #" << rand_mass_ind << " = "
						<< rand_mass << ", rand_age #" << rand_age_ind 
						<< " = " << rand_age << "): expected " 
						<< expected << "+-" << max_error << ", got " << got;
					TEST_ASSERT_MSG(check_diff(expected, got, 1e-12, max_error), 
							test_msg.str().c_str()) 
				}
			}
		}
	} catch(std::exception &ex) {
		TEST_ASSERT_MSG(false, 
			(std::string("test_high_mass_no_death_interp: unexpected "
						 "exception thrown:")+ex.what()).c_str())
	} catch(...) {
		TEST_ASSERT_MSG(false, "test_high_mass_no_death_interp: unexpected "
				"unknown exception thrown:")
	}
}

test_EvolvingStellarQuantity::test_EvolvingStellarQuantity()
{
	TEST_ADD(test_EvolvingStellarQuantity::test_low_mass_lin_age_interp)
	TEST_ADD(test_EvolvingStellarQuantity::test_low_mass_log_age_interp)
	TEST_ADD(test_EvolvingStellarQuantity::test_high_mass_lin_age_interp)
	TEST_ADD(test_EvolvingStellarQuantity::test_high_mass_log_age_interp)
}

#ifdef STANDALONE
int main()
{
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_EvolvingStellarQuantity tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);  
}
#endif
