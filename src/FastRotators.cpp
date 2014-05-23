/*
 * FastRotators.cpp
 *
 *  Created on: Dec 20, 2012
 *      Author: stanley
 */
#include "Common.h"
#include "SimAssumptions.hpp"
#include "AstronomicalConstants.h"
#include "StellarSystem.h"
#include "YRECIO.h"
#include "sampleRotation.h"
#include "sampleHist.h"
#include "OrbitSolver.h"
#include <math.h>
#include <omp.h>
#include <sstream>
#include <exception>
#include <gsl/gsl_rng.h>
class FastRotatorSimulator {
  RotationSampler rotSampler;
  HistSampler smassSampler;
  HistSampler pmassSampler;
  HistSampler pRadiusSampler;
  HistSampler pPeriodSampler;
  YRECEvolution evol;
public:
  FastRotatorSimulator(std::string rot_dist_file,
		       std::string smass_dist_file, std::string pmass_dist_file,
		       std::string pRadius_dist_file, std::string pPeriod_dist_file):
    rotSampler(rot_dist_file), smassSampler(smass_dist_file),
    pmassSampler(pmass_dist_file), pRadiusSampler(pRadius_dist_file),
    pPeriodSampler(pPeriod_dist_file)//, evol("YREC")
  {
    evol.load_state("interp_state_data");
  }
  double iterateOnce(double Q, std::vector<RotationScenario> all_rot,
		     std::ostream* fdata = NULL) {
    /* Randomly create a system, evolve it, and return the fraction of time
     * after MAIN_SEQ_START for which the star's rotation is above
     * SPIN_THRES. Writes orbit data and system data to fdata, if
     * it isn't NULL*/
    using namespace AstroConst;
    double star_mass, planet_mass, planet_radius, P;
#pragma omp critical
    {
      star_mass = smassSampler.randVal();
      planet_mass = pmassSampler.randVal();
      planet_radius = pRadiusSampler.randVal();
      P = pPeriodSampler.randVal();
    }

    double semi = G*star_mass*solar_mass*std::pow(P*day, 2)/4/M_PI/M_PI;
    semi = std::pow(semi, 1.0/3.0)/AU;
    double rot_period = rotSampler.sampleFromDist(star_mass);
    double spin = 2*M_PI/rot_period;

    //select closest rotation scenario
    double min_diff = std::numeric_limits<double>::infinity();
    int best_index = -1;
    for (size_t i=0; i < all_rot.size(); i++) {
      double diff = std::abs(all_rot[i].initSpin - spin);
      if (diff < min_diff) {
	min_diff = diff;
	best_index = i;
      }
    }
	assert(best_index != -1);
    RotationScenario bestModel = all_rot[best_index];
    Star star(star_mass, Q, bestModel.K, WIND_SAT_FREQ, bestModel.Tc,
	      Q_TRANS_WIDTH, bestModel.initSpin, bestModel.Tdisk,
	      evol, 4, 0, 0);
    Planet planet(&star, planet_mass, planet_radius, semi);
    StellarSystem system(&star, &planet);

    double end_age = std::min((const double)MAX_END_AGE, star.lifetime());
    OrbitSolver solver(MIN_AGE, end_age, PRECISION,
		       SPIN_THRES, MAIN_SEQ_START);
    if (fdata != NULL) {
      (*fdata) << "Ms Mp Rp P a0 Prot\n";
      (*fdata) << star_mass<<" "<<planet_mass<<" "<<planet_radius<<
	" "<<P<<" "<<semi<<" "<<rot_period<<"\n";
    }
    double frac = solver.fast_time(system, Inf, PLANET_FORM_AGE,
				   semi)/end_age;
    /*double frac = solver.frac_interesting(
      MIN_AGE, orbit, stop_evolution,
      &system, SPIN_THRES, MAIN_SEQ_START, fdata);*/
    return frac;
  }
  std::vector<double> linspace(double min, double max, int steps) {
	  assert(steps > 1);
    std::vector<double> vals;
    for (int i=0; i < steps; i++) {
      vals.push_back(min + i*(max - min)/(steps - 1));
    }
    return vals;
  }

  void simulateOnce(double Q, RotationScenario rot, double star_mass,
		    double planet_mass, double P0, std::ofstream *ofs,
		    std::string prefix) {
    std::cout.unsetf(std::ios::floatfield);
    std::cout.precision(20);
    //    std::cout << star_mass << " " << planet_mass << " " << P0 << std::endl;

    using namespace AstroConst;
    double semi = G*star_mass*solar_mass*std::pow(P0*day, 2)/4/M_PI/M_PI;
    semi = std::pow(semi, 1.0/3.0)/AU;
    Star star(star_mass, Q, rot.K, WIND_SAT_FREQ, rot.Tc,
	      Q_TRANS_WIDTH, rot.initSpin, rot.Tdisk,
	      evol, 4, 0, 0);
    Planet planet(&star, planet_mass, 0.714, semi);
    StellarSystem system(&star, &planet);
    double end_age = std::min((const double)MAX_END_AGE,
			      star.lifetime());
    OrbitSolver solver(MIN_AGE, end_age, 1e-5,
		       SPIN_THRES, MAIN_SEQ_START);
    try {
      solver(system, Inf, PLANET_FORM_AGE, semi);
    }
    catch (Error::General e) {
      std::cout.unsetf(std::ios::floatfield);
      std::cout.precision(16);
      std::cout << star_mass << " " << planet_mass << " " << P0 << std::endl;
      std::cout << e.what() << std::endl;
    }
    std::valarray<double> ages =
      list_to_valarray(*(solver.get_tabulated_var(AGE)));
    std::valarray<double> semis =
      list_to_valarray(*(solver.get_tabulated_var(SEMIMAJOR)));
    std::valarray<double> semi_derivs =
      list_to_valarray(*(solver.get_tabulated_var_deriv(SEMIMAJOR)));
#pragma omp critical
    {
      (*ofs) << prefix << " Desc" << std::endl;
      double prev_age = -1;
      for (int i=0; i < ages.size(); i++) {
	if (std::isnan(semis[i])) continue;
	if (ages[i] - prev_age > 1e-4)
	  (*ofs) << ages[i] << " " << semis[i] << " "
		 << semi_derivs[i] << std::endl;
	prev_age = ages[i];
      }
    }
  }

  void initDist(std::vector<double> all_Q,
		std::vector<RotationScenario> all_rots,
		std::string filename="all_sims") {
    std::vector<double> smasses = linspace(0.5, 1.1, 10);
    std::vector<double> ln_pmasses = linspace(log(0.05), log(25), 50);
    std::vector<double> P0 = linspace(0.1, 5.9, 100);
    int arrSize = all_Q.size()*all_rots.size()*smasses.size()*P0.size();
    std::ofstream ofs(filename.c_str());

    for (size_t Q_i=0; Q_i < all_Q.size(); Q_i++) {
      for (size_t rot_i=0; rot_i < all_rots.size(); rot_i++) {
	for (size_t smass_i=4; smass_i < smasses.size(); smass_i++) {
	  for (size_t pmass_i=49; pmass_i < ln_pmasses.size();
	       pmass_i++) {
#pragma omp parallel for
	    for (size_t P0_i=92; P0_i < P0.size(); P0_i++) {
	      double Q = all_Q[Q_i];
	      RotationScenario rot = all_rots[rot_i];
	      double star_mass = smasses[smass_i];
	      double planet_mass = exp(ln_pmasses[pmass_i]);
	      double Porb = P0[P0_i];
	      std::stringstream ss;
	      ss << Q_i << " " << rot_i << " " << smass_i <<
		" " << pmass_i << " " << P0_i;
	      simulateOnce(Q, rot, star_mass, planet_mass,
			   Porb, &ofs, ss.str());
	    }
	  }
	}
      }

    }

  }
};

int main() {
  const int numIter = 10000;
  FastRotatorSimulator sim("star_rot_dist", "star_mass_dist",
			   "planet_mass_dist", "planet_radius_dist", "planet_period_dist");
  std::vector<RotationScenario> all_rots;
  all_rots.push_back(RotationScenario(2*M_PI/1.4,0.155,0.012, 0.0025));
  all_rots.push_back(RotationScenario(2*M_PI/7,0.17,0.028, 0.005));
  all_rots.push_back(RotationScenario(2*M_PI/10,0.17,0.030, 0.005));
  std::vector<double> all_Q;
  all_Q.push_back(Q1);
  all_Q.push_back(Q2);
  all_Q.push_back(Q3);
  sim.initDist(all_Q, all_rots);
  return 0;

  for (size_t Q_index=0; Q_index < all_Q.size(); Q_index++) {
    double Q = all_Q[Q_index];
    double counter = 0;
    int numValid = 0;
    int numInvalid = 0;
#pragma omp parallel for
    for (int i=0; i < numIter; i++) {
      bool isValid = true;
      std::stringstream filename, filedata;
      filename << "fastrot/Q" << Q <<"i"<<i;
#pragma omp critical
      if (i%100 == 0)
	std::cout << "count numValid "
		  << counter << " " << numValid << " "
		  << numInvalid << std::endl;
      double frac = -1;
      try {
	frac = sim.iterateOnce(Q, all_rots, &filedata);
      }
      catch(Error::Runtime &e) {
	isValid = false;
#pragma omp critical
	std::cout<<e.get_message()<<std::endl;
      }
      if (frac != 0 && isValid && numValid < 100) {
	std::ofstream ofs(filename.str().c_str());
	ofs << filedata.str();
	ofs << "Frac " << frac;
	ofs.close();
      }
      if (isValid)
	{
#pragma omp atomic
	  counter += frac;
#pragma omp atomic
	  numValid++;
	}
      else {
#pragma omp atomic
	numInvalid++;
      }
    }
    std::cout<<"Q, count, numValid " << Q<<" "<<counter
	     << " " << numValid<< std::endl;
  }
}
