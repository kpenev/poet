#include "sampleRotation.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>

RotationSampler::RotationSampler(std::string filename) {
	/*Initializes object based on rotation-mass relationship in filename.
	 * The file must start with an integer, indicating the number of mass
	 * bins. Then, on another line, it must list the nbins + 1 bin edges.
	 * Starting with the third line, every line must list the rotation periods
	 * (in days) of stars in the corresponding mass bin.*/
	r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, time(NULL));
	setDist(filename);
}
void RotationSampler::setDist(std::string filename) {
	/* Initializes distribution (called by constructor)*/
	std::ifstream ifs(filename.c_str());
	ifs >> nbins;
	assert(nbins > 0);
	bin_edges.resize(nbins + 1, 0);
	rotations.resize(nbins);
	for (int i=0; i < nbins + 1; i++)
		ifs >> bin_edges[i];
	std::string str;
	getline(ifs, str); //get rid of newline
	for (int i=0; i < nbins; i++) {
		getline (ifs,str);
		double period;
		std::stringstream linestream(str);
		while (linestream >> period) {
			rotations[i].push_back(period);
		}
	}
}

double RotationSampler::sampleFromDist(double mass) {
	/*get random rotation period (in days). This is done by locating the mass
	 * bin corresponding to mass, finding a random star within the bin, and
	 * returning its rotation period.*/
	int bin = 0;
	for (int bin=0; bin < nbins; bin++) {
		if (bin == nbins-1 && mass <= bin_edges[bin + 1]) break;
		if (bin != nbins-1 && mass < bin_edges[bin + 1]) break;
	}
	assert(bin != nbins);
	int starIndex = gsl_rng_uniform_int(r, rotations[bin].size());
	return rotations[bin][starIndex];
}

RotationSampler::~RotationSampler() {
	gsl_rng_free(r);
}

