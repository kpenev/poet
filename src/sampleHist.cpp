#include "sampleHist.h"
#include <iostream>
#include <fstream>
#include <assert.h>

HistSampler::HistSampler(std::string filename) {
	/*Initializes object based on histogram in filename. The file must be of
	 * the following format.  It must start with an integer, indicating the
	 * number of bins in the histogram. The next line must list all the bin
	 * edges. For example, if there are 3 bins, and the bin edges were 0, 1,
	 * 3, 4, the bins are [0, 1), [1, 3), and [3, 4]. The line after that must
	 * list the values in all the bins, and they must add up to 1.*/
	r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, 1512);
	setDist(filename);
}
void HistSampler::setDist(std::string filename) {
	//in constructor: initialize distribution based on file
	std::ifstream ifs(filename.c_str());
	int nbins;
	ifs >> nbins;
	assert(nbins > 0);
	dist.resize(nbins, 0);
	bin_edges.resize(nbins + 1, 0);
	for (int i=0; i < nbins + 1; i++)
		ifs >> bin_edges[i];
	for (int i=0; i < nbins; i++)
		ifs >> dist[i];
}

int HistSampler::sampleFromDist(std::vector<double> dist) {
	/*Takes probability distribution dist, and randomly chooses an index from
	 * it. dist must be normalized.*/
	double rand = gsl_rng_uniform(r);
	double cumulProb = 0;
	for (int i=0; i < dist.size(); i++) {
		cumulProb += dist[i];
		if (rand < cumulProb) return i;
	}
	return dist.size() - 1;
}

double HistSampler::randVal() {
	/* Chooses a random value from the histogram. First randomly chooses a
	 * bin, then chooses a value in between the two bin edges uniformly at
	 * random. */
	int index = sampleFromDist(dist);
	double leftEdge = bin_edges[index];
	double rightEdge = bin_edges[index + 1];
	return leftEdge + (rightEdge - leftEdge)*gsl_rng_uniform(r);
}

HistSampler::~HistSampler() {
	gsl_rng_free(r);
}
