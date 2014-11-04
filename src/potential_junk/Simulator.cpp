/*
 * Simulator.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: stanley
 */
#include "StellarDissipation.h"
#include "Planet.h"
#include "StellarSystem.h"
#include "OrbitSolver.h"
#include "YRECIO.h"
#include <iostream>
#include <sstream>
#include <omp.h>
class Simulator {
private:
	double wconv_sun, wind_sat_freq;
	static const double Q_trans_width=1e-3;
	double K, Tc;
	static const double min_age=0.005, precision=0.001;
	double age, star_mass, star_radius, conv_spin, tidal_Q;
	double planet_mass, planet_radius, semimajor;
	Star* star;
	Planet* planet;
	StellarSystem* system;
	OrbitSolver* solver;
public:
	Simulator(double age, double star_mass, double star_radius,
			double conv_spin, double K, double Tc,
			double tidal_Q, double planet_mass, double planet_radius,
			double semimajor):
		age(age), star_mass(star_mass), star_radius(star_radius),
		conv_spin(conv_spin), tidal_Q(tidal_Q), planet_mass(planet_mass),
		planet_radius(planet_radius), semimajor(semimajor){

		wconv_sun = 2*M_PI/25.6;
		wind_sat_freq = 10*wconv_sun;

		YRECEvolution evol;
		evol.load_state("interp_state_data");
		star = new Star("Test", age, star_mass, star_radius, conv_spin,
				conv_spin, tidal_Q, K, wind_sat_freq,
				Tc, Q_trans_width, evol);
		planet = new Planet(star, planet_mass, planet_radius, semimajor);
		system = new StellarSystem(star, planet);
		solver = new OrbitSolver(min_age, age, precision,
				stellar_system_diff_eq, stellar_system_jacobian);
	}

	void sim_planet(double init_spin, double init_a) {
		/*Simulate the orbit of a planet from min_age to present day, and
		 * print out (age, semimajor axis (AU), wconv (rad/day),
		 * wrad (rad/day) */
		double start_Lc = star->moment_of_inertia(min_age, convective)*init_spin;
		double start_Lr = star->moment_of_inertia(min_age, radiative)*init_spin;
		double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;
		double present_orbit[3]={std::pow(AU_Rsun*init_a, 6.5), start_Lc,
								start_Lr};
		(*solver)(min_age, present_orbit, stop_evolution, system);
		std::list<double>::const_iterator ages =
				(solver->get_tabulated_indep_var())->begin();
		std::list<double>::const_iterator a =
				(solver->get_tabulated_dep_var(0))->begin();
		std::list<double>::const_iterator Lc =
				(solver->get_tabulated_dep_var(1))->begin();
		std::list<double>::const_iterator Lr =
				(solver->get_tabulated_dep_var(2))->begin();
		for(; ages != solver->get_tabulated_indep_var()->end(); ages++, a++, Lc++, Lr++) {
			std::cout<< *ages << " "<<std::pow(*a,1.0/6.5)/AU_Rsun << " " <<
					*Lc/star->moment_of_inertia(*ages, convective) << " "<<
					*Lr/star->moment_of_inertia(*ages, radiative) <<
					std::endl;
		}
	}

	void write_matrix(double min_spin, double max_spin, double min_semi,
			double max_semi, int semi_steps, int spin_steps,
			std::string filename) {
		/*Outputs a matrix, semi_steps x spin_steps in size, that shows how
		 * good each initial condition (a0, w0) is in reproducing current
		 * conditions. Specifically, outputs L2 norm of fractional error, or
		 * 100 if the planet dies before the current age can be reached.
		 * Output file is in following format.
		 * First line: a0 min_semi max_semi semi_steps
		 * Second line: w0 min_spin max_spin spin_steps
		 * (matrix itself, one semimajor axis per line) */
		using namespace std;
		using namespace AstroConst;
		ofstream ofs(filename.c_str());
		ofs <<"a0 "<<min_semi<<" "<<max_semi<<" "<<semi_steps<<endl;
		ofs <<"w0 "<<min_spin<<" "<<max_spin<<" "<<spin_steps<<endl;
		double final_Lc = star->moment_of_inertia(age, convective)*conv_spin;
		double init_Ic = star->moment_of_inertia(min_age, convective);
		double init_Ir = star->moment_of_inertia(min_age, radiative);
		for (int semi_j=0; semi_j < semi_steps; semi_j++) {
			double* results = new double[semi_steps];
			#pragma omp parallel for
			for (int spin_i=0; spin_i < spin_steps; spin_i++) {
				double spin = min_spin +
						spin_i*(max_spin - min_spin)/spin_steps;
				double semi = min_semi +
						semi_j*(max_semi - min_semi)/semi_steps;
				double start_Lc = init_Ic*spin;
				double start_Lr = init_Ir*spin;
				double present_orbit[3]={std::pow(semi, 6.5), start_Lc,
						start_Lr};
				OrbitSolver parallelSolver(min_age, age, precision,
						stellar_system_diff_eq, stellar_system_jacobian);
				present_orbit[0] *=
						std::pow(AstroConst::AU/AstroConst::solar_radius,6.5);
				parallelSolver(min_age, present_orbit, stop_evolution, system);
				double percent_error =
						parallelSolver.last_error(age, semimajor, final_Lc);
				results[spin_i] = percent_error;
			}
			for (int spin_i=0; spin_i < spin_steps; spin_i++) {
				if (std::isnan(results[spin_i])) ofs<<"100 ";
				else ofs<<std::log10(results[spin_i])<<" ";
			}
			delete[] results;
			ofs<<endl;
		}
		ofs.close();
	}

	~Simulator() {
		delete star;
		delete planet;
		delete system;
		delete solver;
	}
};

int main(int argc, char** argv) {
	double age = 4.2;
	double star_mass=1.13,  star_radius=1.21,
			conv_spin=1.11, tidal_Q=1e8, planet_mass=7.69,
			planet_radius=1.09, semimajor=0.02691;
	//double K=0.155, Tc=0.012;
	//double K=0.17, Tc=0.028;
	double K=0.17, Tc=0.03;
	Simulator s(age,star_mass,star_radius,conv_spin,K,Tc,tidal_Q,
			planet_mass,planet_radius,semimajor);
	std::stringstream ss;
	ss << "corot_14zoom_Q" << tidal_Q << "age" << age << "Tc" << Tc;
	s.write_matrix(0.63,6.3,0.03,0.05,500,500, ss.str());
	return 0;

}
