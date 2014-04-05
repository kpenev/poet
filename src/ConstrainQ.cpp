/*
 * ConstrainQ.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: stanley
 */
#include "StellarSystem.h"
#include "StellarQ.h"
#include "Planet.h"
#include "YRECIO.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <initializer_list>
#include <vector>
#include <omp.h>

struct RotationScenario {
	double initSpin;
	double K;
	double Tc;
	RotationScenario(double initSpin, double K, double Tc):
		initSpin(initSpin), K(K), Tc(Tc) {}
};

struct SystemData {
	std::string name;
	double mp, rp, a, ms, rs, Prot;
	double min_age, nom_age, max_age;
};
enum AgeChoice {MIN_AGE, NOM_AGE, MAX_AGE};
const AgeChoice age_choice = NOM_AGE;
std::vector<SystemData> read_file(std::string filename) {
	std::vector<SystemData> allSystems;
	std::ifstream ifs(filename.c_str());
	SystemData temp;
	while (ifs >> temp.name >> temp.mp >> temp.rp >> temp.a >> temp.ms
			>> temp.rs >> temp.Prot >>
			temp.min_age >> temp.nom_age >> temp.max_age) {
		allSystems.push_back(temp);
	}
	ifs.close();
	return allSystems;
}
const double wind_sat_freq = 2.454;
const double start_age = 0.005;
const double min_a0 = 0.01, max_a0 = 0.1,
		min_w0 = 0.63, max_w0 = 6.3;
const double max_err = 0.01;
const double Q_transition_width = 1e-3;
int main(int argc, char** argv) {
	/* For each system in argv[1], and for each Q_*, calculates whether any
	 * initial condition could reproduce current observations. File must be in
	 * following format (one line per system):
	 * Name mp rp a ms rs Prot min_age nom_age max_age*/
	std::vector<double> all_Q;
	all_Q.push_back(1e5);
	all_Q.push_back(5e5);
	all_Q.push_back(1e6);
	all_Q.push_back(5e6);
	all_Q.push_back(1e7);
	all_Q.push_back(5e7);
	all_Q.push_back(1e8);

	std::vector<RotationScenario> all_rots;
	all_rots.push_back(RotationScenario(2*M_PI/0.889,0.155,0.012));
	all_rots.push_back(RotationScenario(2*M_PI/7,0.17,0.028));
	all_rots.push_back(RotationScenario(2*M_PI/10,0.17,0.030));

	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " filename" << std::endl;
		exit(-1);
	}
	std::vector<SystemData> allSystems = read_file(argv[1]);
	YRECEvolution evol;
	evol.load_state("interp_state_data");
	for (size_t Q_index=0; Q_index < all_Q.size(); Q_index++) {
		double Q = all_Q[Q_index];
		std::cout << Q << std::endl;
		std::vector<bool> explained(allSystems.size(), false);
		for (size_t rot_i=0; rot_i < all_rots.size(); rot_i++) {
			double Tc = all_rots[rot_i].Tc;
			double K = all_rots[rot_i].K;
			#pragma omp parallel for
			for (size_t sysIndex=0; sysIndex < allSystems.size(); sysIndex++) {
				SystemData sys = allSystems[sysIndex];
				double age;
				if (age_choice == MIN_AGE) age = sys.min_age;
				else if (age_choice == NOM_AGE) age = sys.nom_age;
				else age = sys.max_age;
				if (age > 9*std::pow(sys.ms, -3)) age = 9*std::pow(sys.ms, -3);
				bool success = false;
				try {
				double conv_spin = 2*M_PI/sys.Prot;
				Star star(sys.name, age, sys.ms, sys.rs, conv_spin,
						conv_spin, Q, K, wind_sat_freq,
						Tc, Q_transition_width, evol);
				Planet planet(&star, sys.mp, sys.rp, sys.a);
				StellarSystem system(&star, &planet);
				success = system.anneal_solve(false, min_a0, max_a0,
						min_w0, max_w0, max_err);
				}
				catch(Error::BadFunctionArguments &e) {
					std::cout << e.get_message() << std::endl;
				}
				if (success) explained[sysIndex] = true;
			}
		}
		int numExplained = 0;
		for (size_t sysIndex=0; sysIndex < allSystems.size(); sysIndex++) {
			std::cout << explained[sysIndex] << " ";
			if (explained[sysIndex]) numExplained++;
		}
		std::cout << std::endl;
		std::cout << numExplained << std::endl;

	}
}
