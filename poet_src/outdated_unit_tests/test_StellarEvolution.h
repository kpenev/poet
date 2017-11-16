/**\file
 *
 * \brief Defines the test suite that exercises the StellarEvolution class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_STELLAR_EVOLUTION
#define __TEST_STELLAR_EVOLUTION

#include "Common.h"
#include "../StellarEvolution.h"
#include "../Error.h"
#include "test_EvolvingStellarQuantity.h"
#include <sstream>
#include <valarray>
#include <vector>

///\brief The test suite that exercises the StellarEvolution class.
///
///\ingroup UnitTests_group
class test_StellarEvolution : public Test::Suite {
private:

protected:
	///No fixtures at this time.
	void setup() {};

	///No fixtures at this time.
	void tear_down() {};

	///\brief Compares an exact polynomial track with an interpolated
	///quantity.
	///
	///Deletes both the quantity and the track when done.
	///
	///The maximum allowed error is either 1e-12 (if closest_tracks is NULL)
	///or the difference  between the exact and closest track
	void check_evolution(const PolynomialEvolutionTrack *track,
				const EvolvingStellarQuantity *quantity,
				const std::string &quant_name, double test_mass,
				std::valarray< std::valarray<double> > poly_coef,
				PolynomialEvolutionTrack *closest_track=NULL);

	///\brief Checks the age scaling by using polynomial evolutions.
	///
	///Creates a stellar evolution based on polynomial expressions of all 
	///quantities imposing the given age scalings for low and high mass
	///tracks, no smoothing is applied to any stellar quantity and tests that
	///correct results are returned.
	void polynomial_evolution_check(double low_mass_age_scaling=0,
			double high_mass_age_scaling=0);

	///\brief Tests the stellar evolution of polynomial quantities with no
	///smoothing and no mass scaling of age.
	void test_polynomial_evolution()
	{polynomial_evolution_check();}

	///\brief Tests the stellar evolution of polynomial quantities with no
	///smoothing, but low mass tracks age scaling of 1 and high mass tracks
	///age scaling of 2.
	void test_scaled_polynomial_evolution()
	{polynomial_evolution_check(9.5, -4.3);}
public:
	///Add all tests to the test suite.
	test_StellarEvolution();

	///Do nothing.
	~test_StellarEvolution() {}
};

#endif
