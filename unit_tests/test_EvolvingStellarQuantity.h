#ifndef __TEST_EVOLVING_STELLAR_QUANTITY
#define __TEST_EVOLVING_STELLAR_QUANTITY

#include "Common.h"
#include "../StellarEvolution.h"
#include "../Error.h"
#include <sstream>

class test_EvolvingStellarQuantity : public Test::Suite {
private:

	///Returns the track with mass closest to the given mass.
	const OneArgumentDiffFunction *closest_track(double m,
			const std::valarray<double> &track_masses,
			const std::list<const OneArgumentDiffFunction *> &tracks);


	///Tests interpolation over low stellar masses with coiniciding 
	///lifetimes, never asking for ages outside the lifetimes. 
	void test_low_mass_no_death_interp(bool log_age=false);

	///Tests interpolation over low stellar masses with coiniciding 
	///lifetimes, never asking for ages outside the lifetimes using
	///linear age as the independent argument. 
	void test_low_mass_lin_age_interp()
	{test_low_mass_no_death_interp(false);}

	///Tests interpolation over low stellar masses with coiniciding 
	///lifetimes, never asking for ages outside the lifetimes using
	///log(age) as the independent argument. 
	void test_low_mass_log_age_interp()
	{test_low_mass_no_death_interp(true);}

	///Tests the interpolation over high stellar masses with coinciding 
	///lifetimes, never asking for ages outside the lifetimes.
	void test_high_mass_no_death_interp(bool log_age=false);

	///Tests interpolation over high stellar masses with coiniciding 
	///lifetimes, never asking for ages outside the lifetimes using
	///linear age as the independent argument. 
	void test_high_mass_lin_age_interp()
	{test_high_mass_no_death_interp(false);}

	///Tests interpolation over high stellar masses with coiniciding
	///lifetimes, never asking for ages outside the lifetimes using
	///log(age) as the independent argument. 
	void test_high_mass_log_age_interp()
	{test_high_mass_no_death_interp(true);}
protected:
	///No fixtures at this time.
	void setup() {};

	///No fixtures at this time.
	void tear_down() {};
public:
	///Add all tests to the test suite.
	test_EvolvingStellarQuantity();

	///Do nothing.
	~test_EvolvingStellarQuantity() {}
};

#endif
