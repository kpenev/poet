/*
 * test_Star.h
 *
 *  Created on: Oct 29, 2012
 *      Author: stanley
 */

#ifndef TEST_STAR_H_
#define TEST_STAR_H_

#include "Star.h"
#include "Common.h"
#include <assert.h>



class test_Star : public Test::Suite {
private:
	Star* star;
	StarData d;
public:
	test_Star();
	void test_basic();
	void test_known();
	void test_ang_momentum();
	void test_tidal();
	void test_evolution();
	void test_spin();
	void test_wind_torque();
	void test_diff_rot();

};
#endif /* TEST_STAR_H_ */
