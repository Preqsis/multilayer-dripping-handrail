#ifndef OBSERVATION_HPP
#define OBSERVATION_HPP

#include <vector>

class ArgumentParser;

namespace Observation {

void terminate(std::vector<size_t> dim, int n_workers);

double filter_gauss(double**** data, std::vector<size_t> dim, double mu, double fwhm);

// Function for "sim" MPI slave processes
void slave(std::vector<size_t> dim_spec, std::vector<size_t> dim_obs, ArgumentParser* p);

// Function for "sim" MPI master process
void master(std::vector<size_t> dim_spec, std::vector<size_t> dim_obs, int n_workers, ArgumentParser* p);

}

#endif

