/**\file
 *
 * \brief Defines the test methods of the test suite that exercises the
 * StellarEvolution class.
 *
 * \ingroup UnitTests_group
 */

#include "test_StellarEvolution.h"

void test_StellarEvolution::check_evolution(
		const PolynomialEvolutionTrack *track, 
		const EvolvingStellarQuantity *quantity,
		const std::string &quant_name, double test_mass,
		std::valarray< std::valarray<double> > poly_coef,
		PolynomialEvolutionTrack *closest_track)
{
	const unsigned num_tests=100;
	for(unsigned test=0; test<num_tests; test++) {
		double rand_age=rand_value(quantity->range_low(),
				quantity->range_high());
		std::ostringstream test_message;
		double expected=(*track)(rand_age);
		double got=(*quantity)(rand_age);
		double max_diff=std::abs(closest_track ?
					(*closest_track)(rand_age)-expected : 1e-8);
		test_message << "Testing interpolation of " << quant_name 
			<< " to M=" << test_mass << " for age=" << rand_age
			<< " exact track: " << *track << " derived from: "
			<< poly_coef << ", expected " << expected  << "+-" << max_diff
			<< ", got " << got;
		TEST_ASSERT_MSG(check_diff(expected, got, 1e-3, max_diff),
				test_message.str().c_str());
		for(unsigned deriv_order=0; deriv_order<2; deriv_order++) {
			const FunctionDerivatives *q_deriv=quantity->deriv(rand_age);
			expected=track->order(deriv_order);
			got=q_deriv->order(deriv_order);
			if(closest_track) 
				max_diff=std::abs(closest_track->order(deriv_order)-expected);
			else max_diff=1e-8;
			test_message.str("");
			test_message << "Testing order " << deriv_order 
				<< " derivative of " << quant_name 
				<< " to M=" << test_mass << " for age=" << rand_age
				<< " exact track: " << *track << " derived from: "
				<< poly_coef << ", expected " << expected  << "+-" << max_diff
				<< ", got " << got;
			TEST_ASSERT_MSG(check_diff(expected, got, 1e-3, max_diff), 
					test_message.str().c_str()) 
		}
	}
	delete track;
	delete quantity;
	if(closest_track) delete closest_track;
}

void test_StellarEvolution::polynomial_evolution_check(
		double low_mass_age_scaling, double high_mass_age_scaling)
{
	const unsigned num_tests=1;
	try {
		std::valarray<double> masses(100), common_ages(100);
		std::valarray< std::valarray<double> >
			r_poly_coef=rand_poly_coef(), Iconv_poly_coef=rand_poly_coef(),
			Irad_poly_coef=rand_poly_coef(), Mrad_poly_coef=rand_poly_coef(),
			Rcore_poly_coef, Mconv_poly_coef=-Mrad_poly_coef,
			Itot_poly_coef=Iconv_poly_coef+Irad_poly_coef;
		Mconv_poly_coef[0][1]+=1;
		for(unsigned i=0; i<masses.size(); i++) 
			masses[i]=min_stellar_mass + 
				((max_stellar_mass-min_stellar_mass)/masses.size())*i;
		for(unsigned i=0; i<common_ages.size(); i++)
			common_ages[i]=min_age + ((max_age-min_age)/common_ages.size())*i;
		PolynomialStellarEvolution evolution(masses, common_ages, r_poly_coef,
				Iconv_poly_coef, Itot_poly_coef, Mrad_poly_coef,
				Rcore_poly_coef, low_mass_age_scaling,
				high_mass_age_scaling);
		for(unsigned test=0; test<num_tests; test++) {
			double rand_mass=rand_value(masses[0], masses[masses.size()-1]);
			bool high_mass=rand_mass>max_low_mass;
			double closest_mass = -1;
			if(high_mass) {
				double *below_mass=&masses[masses.size()-1];
				while(*below_mass>rand_mass) below_mass--;
				if(rand_mass-*below_mass<*(below_mass+1)-rand_mass)
					closest_mass=*below_mass;
				else closest_mass=*(below_mass+1);
			}
			TEST_ASSERT(closest_mass > 0);
			check_evolution(exact_track(r_poly_coef, rand_mass,
						low_mass_age_scaling, high_mass_age_scaling),
					evolution.interpolate_radius(rand_mass), "R", rand_mass,
					r_poly_coef,
					(high_mass ? exact_track(r_poly_coef, closest_mass,
											 low_mass_age_scaling, 
											 high_mass_age_scaling,
											 rand_mass) : NULL));
			check_evolution(exact_track(Iconv_poly_coef, rand_mass,
						low_mass_age_scaling, high_mass_age_scaling),
					evolution.interpolate_moment_of_inertia(rand_mass, 
						convective), "Iconv", rand_mass, Iconv_poly_coef,
					(high_mass ? exact_track(Iconv_poly_coef, closest_mass,
											 low_mass_age_scaling,
											 high_mass_age_scaling,
											 rand_mass) : NULL));
			check_evolution(exact_track(Itot_poly_coef, rand_mass,
						low_mass_age_scaling, high_mass_age_scaling),
					evolution.interpolate_moment_of_inertia(rand_mass, 
						total), "Itot", rand_mass, Itot_poly_coef,
					(high_mass ? exact_track(Itot_poly_coef, closest_mass,
											 low_mass_age_scaling,
											 high_mass_age_scaling,
											 rand_mass) : NULL));
			check_evolution(exact_track(Mrad_poly_coef, rand_mass,
						low_mass_age_scaling, high_mass_age_scaling),
					evolution.interpolate_zone_mass(rand_mass, 
						radiative), "Mrad", rand_mass, Mrad_poly_coef,
					(high_mass ? exact_track(Mrad_poly_coef, closest_mass,
											 low_mass_age_scaling,
											 high_mass_age_scaling,
											 rand_mass) : NULL));
			check_evolution(exact_track(Rcore_poly_coef, rand_mass,
						low_mass_age_scaling, high_mass_age_scaling),
					evolution.interpolate_core_boundary(rand_mass), "Rcore",
					rand_mass, Rcore_poly_coef,
					(high_mass ? exact_track(Rcore_poly_coef, closest_mass,
											 low_mass_age_scaling,
											 high_mass_age_scaling,
											 rand_mass) : NULL));
		}
	} catch(Error::General &ex) {
		TEST_ASSERT_MSG(false, 
			(std::string("test_polynomial_evolution: unexpected exception"
						 "thrown: ")+ex.what()+": "+ex.get_message()).c_str())
	} catch(std::exception &ex) {
		TEST_ASSERT_MSG(false, 
			(std::string("test_polynomial_evolution: unexpected exception"
						 "thrown: ")+ex.what()).c_str())
	} catch(...) {
		TEST_ASSERT_MSG(false, "test_polynomial_evolution: unexpected "
				"unknown exception thrown:")
	}
}

test_StellarEvolution::test_StellarEvolution()
{
	TEST_ADD(test_StellarEvolution::test_polynomial_evolution)
	TEST_ADD(test_StellarEvolution::test_scaled_polynomial_evolution)
}

#ifdef STANDALONE
int main()
{
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_StellarEvolution tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);  
}
#endif
