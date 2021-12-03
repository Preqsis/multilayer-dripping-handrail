#ifndef SIMULATION_HPP
#define SIMULATION_HPP

namespace Simulation {

// dataset write check function
bool writable(bool wf, bool wl, uint s, uint sf, uint sl);

// 'Empty' simulation grid initalization
void grid_init(double*** data, std::vector<size_t> dim, ArgumentParser* p);

// Function for "sim" MPI slave processes
void slave(std::vector<size_t> dim, ArgumentParser* p);

// Function for "sim" MPI master process
void master(std::vector<size_t> comm_dim, int n_workers, ArgumentParser* p);

}

#endif

