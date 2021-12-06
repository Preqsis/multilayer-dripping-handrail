#ifndef OBSERVATION_HPP
#define OBSERVATION_HPP

namespace Observation {

// Function for "sim" MPI slave processes
void slave(std::vector<size_t> dim, ArgumentParser* p);

// Function for "sim" MPI master process
void master(std::vector<size_t> cdim, int n_workers, ArgumentParser* p);

}

#endif

