#include "Common.h"
#include "SimAssumptions.hpp"
#include "Star.h"
#include "Planet.h"
#include "StellarSystem.h"
#include "OrbitSolver.h"
#include "YRECIO.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <omp.h>
#include <time.h>

void simulateOnce(double Q, RotationScenario rot, double star_mass,
		double planet_mass, double P0)
{
	using namespace AstroConst;
	try {
		YRECEvolution evol;
		evol.load_state("interp_state_data_phs4");
		double semi = G*star_mass*solar_mass*std::pow(P0*day, 2)/4/M_PI/M_PI;
		semi = std::pow(semi, 1.0/3.0)/AU;
		Star star(star_mass, Q, rot.K, WIND_SAT_FREQ, rot.Tc,
				Q_TRANS_WIDTH, rot.initSpin, rot.Tdisk,
				evol, 4, 0, 0);
		Planet planet(&star, planet_mass, 0.714, semi);
		StellarSystem system(&star, &planet);
		double end_age = std::min((const double)MAX_END_AGE,
				star.get_lifetime());
		OrbitSolver solver(MIN_AGE, end_age, 1e-5,
				SPIN_THRES, MAIN_SEQ_START);
		time_t start_time, end_time;
		time(&start_time);
		solver(system, Inf, PLANET_FORM_AGE, semi);
		time(&end_time);
		std::cout << "Solver took " << difftime(end_time, start_time)
			<< " seconds to finish." << std::endl;
		std::cout << "Finished!" << std::endl;
	} catch (Error::General err) {
		std::cout << err.what() << ": " << err.get_message() << std::endl;
	}
}


int main(int argc, char** argv) {
	using namespace AstroConst;
	std::vector<RotationScenario> all_rots;
	all_rots.push_back(RotationScenario(2*M_PI/1.4,0.155,0.012, 0.0025));
	all_rots.push_back(RotationScenario(2*M_PI/7,0.17,0.028, 0.005));
	all_rots.push_back(RotationScenario(2*M_PI/10,0.17,0.030, 0.005));

	//extremely slow evol
	double smass=0.5, pmass=25, Q=1e6, P0=5.84141414;
/*	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	//bad stop hostory interval
	smass=1.0333333; pmass=0.631816; Q=1e6; P0=0.890909;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	//Misses planet death
	smass=0.5; pmass=0.3804226; Q=1e6; P0=0.94848;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	//The following were reported by Michael to throw BadFunctionArguments
	//but they appear to work
	smass=0.5; pmass=17.0882; Q=1e6; P0=4.78687;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.5; pmass=22.0221; Q=1e6; P0=4.84545;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.566667; pmass=25; Q=1e6; P0=2.03333;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.566667; pmass=25; Q=1e6; P0=3.20505;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.566667; pmass=25; Q=1e6; P0=3.26364;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.633333; pmass=15.0528; Q=1e6; P0=5.37273;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.633333; pmass=19.3989; Q=1e6; P0=3.96667;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.633333; pmass=19.3989; Q=1e6; P0=4.08384;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.633333; pmass=19.3989; Q=1e6; P0=4.37677;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.633333; pmass=25; Q=1e6; P0=2.97071;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);
	//End of BadFunctionArgument throwers
	
	smass=0.5; pmass=17.088224742562282; Q=1e6; P0=4.7868686868686874;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	//planet starts inside Roche radius - should fail!
	smass=1.0333333; pmass=0.631816; Q=1e6; P0=0.2;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.5;		pmass=19.3989;		P0=4.08384;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.566667;	pmass= 19.3989;		P0=5.0798;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.566667;	pmass=25;			P0= 3.26364;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.7;		pmass= 15.0528;		P0=4.6697;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.7;		pmass= 17.0882;		P0=0.744444;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.7;		pmass= 17.0882;		P0=2.44343;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.766667;	pmass= 0.431864;	P0= 3.26364;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.766667;	pmass= 19.3989;		P0= 1.21313;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.766667;	pmass= 22.0221;		P0= 3.14646;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.766667; pmass=25;			P0=5.4899;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.5000000000000000; pmass=19.3989385656805; P0=4.083838383838384;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.5666666666666667; pmass=19.3989385656805; P0=5.079797979797980;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.5666666666666667; pmass=25.0000000000000; P0=3.263636363636364;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	//Slow as of 20131004
	smass=0.90000000000000013323; pmass=25; P0=5.9000000000000003553;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);*/

	//Slow as of 20131004
	smass=0.90000000000000010; pmass=10.000000000000010; P0=4.611111111111112;
	simulateOnce(Q, all_rots[0], smass, pmass, P0);

	smass=0.96666666666666679; pmass=7.2297055951858527; P0=2.5020202020202027;
//	simulateOnce(Q, all_rots[0], smass, pmass, P0);
	return 0;
}
