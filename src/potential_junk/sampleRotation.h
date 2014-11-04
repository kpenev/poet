/*
 * sampleRotation.h
 *
 *  Created on: Dec 20, 2012
 *      Author: stanley
 */

#ifndef SAMPLEROTATION_H_
#define SAMPLEROTATION_H_
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
class RotationSampler {
	int nbins;
	std::vector<double> bin_edges;
	std::vector<std::vector<double> > rotations;
	gsl_rng* r;
	void setDist(std::string filename);
public:
	RotationSampler(std::string filename);
	double sampleFromDist(double mass);
	~RotationSampler();
};

#endif /* SAMPLEROTATION_H_ */
