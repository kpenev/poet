/*
 * sampleMass.h
 *
 *  Created on: Dec 20, 2012
 *      Author: stanley
 */

#ifndef SAMPLEHIST_H_
#define SAMPLEHIST_H_
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
class HistSampler {
	std::vector<double> bin_edges;
	std::vector<double> dist;
	gsl_rng * r;
	void setDist(std::string filename);
	int sampleFromDist(std::vector<double> dist);
public:

	HistSampler(std::string filename);

	double randVal();

	~HistSampler();
};

#endif /* SAMPLEHIST_H_ */
