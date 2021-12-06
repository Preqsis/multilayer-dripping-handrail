#include <iostream>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include "argparse-cpp/ArgumentParser.h"

#include "Functions.hpp"
namespace fn = Functions;

#include "Constants.hpp"
namespace cs = Constants;

#include "Observation.hpp"

// Function for "sim" MPI slave processes
void Observation::slave(std::vector<size_t> dim, ArgumentParser* p) {}

// Function for "sim" MPI master process
void Observation::master(std::vector<size_t> cdim, int n_workers, ArgumentParser* p) {}

