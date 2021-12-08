#ifndef RADIATION_HPP
#define RADIATION_HPP

#include <vector>
#include "argparse-cpp/ArgumentParser.h"

namespace Radiation {

void terminate(std::vector<size_t> dim, int n_workers);

double planck(double wl, double T);

void slave(std::vector<size_t> dim_mass, std::vector<size_t> dim_spec, ArgumentParser* p);

void master(std::vector<size_t> dim_mass, std::vector<size_t> dim_spec, int n_workers, ArgumentParser* p);

}

#endif

