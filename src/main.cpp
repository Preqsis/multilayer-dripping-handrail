#include <iostream>
#include <mpi.h>
#include <math.h>

#include "argparse-cpp/Argument.h"
#include "argparse-cpp/ArgumentParser.h"

#include "Functions.hpp"
namespace fn = Functions;

#include "Constants.hpp"
namespace cs = Constants;

#include "Simulation.hpp"
#include "Radiation.hpp"
#include "Observation.hpp"

// simulation task "branch"
void sim(int rank, int n_workers, ArgumentParser* p) {
    // Comms dimensions 
    size_t n_jobs               = p->i("idim") * p->i("jdim"); // number of jobs (cells)
    std::vector<size_t> cdim    = {(size_t)std::ceil((double)n_jobs / (double)n_workers), 13};

    // Process specific call
    if (rank == cs::mpi::MASTER) {
        Simulation::master(cdim, n_workers, p);
    } else {
        Simulation::slave(cdim, p);
    }
}

// radiation task "branch"
void rad(int rank, int n_workers, ArgumentParser* p) {
    // Comms dimensions
    std::vector<size_t> dim_mass = {(size_t)p->i("idim"), (size_t)p->i("jdim"), 13};
    std::vector<size_t> dim_spec = {dim_mass[0], dim_mass[1], (size_t)((p->d("wl_high") - p->d("wl_low")) / p->d("wl_step") + 1), 3};

    // Process specific call
    if (rank == cs::mpi::MASTER) {
        Radiation::master(dim_mass, dim_spec, n_workers, p);
    } else {
        Radiation::slave(dim_mass, dim_spec, p);
    }
}

// observation task "branch"
void obs(int rank, int n_workers, ArgumentParser* p) {
    std::vector<size_t> dim_mass = {(size_t)p->i("idim"), (size_t)p->i("jdim"), 13};

    std::vector<size_t> dim_spec = {dim_mass[0], dim_mass[1], (size_t)((p->d("wl_high") - p->d("wl_low")) / p->d("wl_step") + 1), 3};
    std::vector<size_t> dim_obs = {5};

    // Process specific call
    if (rank == cs::mpi::MASTER) {
        Observation::master(dim_spec, dim_obs, n_workers, p);
    } else {
        Observation::slave(dim_spec, dim_obs, p);
    }
}

void ArgumentParserInit(ArgumentParser* p) {
    // Enable verbosity
    Argument<bool>* verbose = new Argument<bool>("verbose", false);
    verbose->setShorthand("v");
    verbose->setHelp("Enables verbose mode.");
    p->addArgument(verbose);

    // Run mass distribution
    Argument<bool>* asim = new Argument<bool>("sim", false);
    asim->setHelp("Run mass distribution task.");
    p->addArgument(asim); // run mass distribution sim

    // Run radiation
    Argument<bool>* arad = new Argument<bool>("rad", false);
    arad->setHelp("Run radiation task.");
    p->addArgument(arad); // run radiation output computation

    // Run obseravation
    Argument<bool>* aobs = new Argument<bool>("obs", false);
    aobs->setHelp("Run mass observation task.");
    p->addArgument(aobs); // run observation and filtering

    // Number of sim steps
    p->addArgument(new Argument<int>("n", 5e5, "Number of simulation steps.")); // number of simulation steps
    p->addArgument(new Argument<int>("step_first", 0)); // first step in range
    p->addArgument(new Argument<int>("step_last", 99)); // last step in range

    // Specify steps range to save
    Argument<int>* save_start = new Argument<int>("save_start");
    save_start->setHelp("First step to save.");
    p->addArgument(save_start);
    Argument<int>* save_end = new Argument<int>("save_end");
    save_end->setHelp("Last step to save.");
    p->addArgument(save_end);

    // Disk dimensions
    Argument<int>* idim = new Argument<int>("idim"); // number of layers
    idim->setRequired(true);
    idim->setHelp("Number of layers(aka. rings).");
    p->addArgument(idim);
    Argument<int>* jdim = new Argument<int>("jdim"); // number of cells in each layer
    jdim->setRequired(true);
    jdim->setHelp("Number of cells in each ring.");
    p->addArgument(jdim);

    // Input / output defs
    Argument<std::string>* outdir = new Argument<std::string>("outdir"); // data output directory
    outdir->setRequired(true);
    outdir->setShorthand("o");
    outdir->setHelp("Output data directory.");
    p->addArgument(outdir);                                     
    Argument<std::string>* sim_file = new Argument<std::string>("sim_file");  // input sim_file
    sim_file->setRequired(false);
    sim_file->setHelp("Input HDF5 mass data file.");
    p->addArgument(sim_file);                                     
    Argument<std::string>* sim_dkey = new Argument<std::string>("sim_dkey"); // initial sim dkey
    sim_dkey->setHelp("Initial data key in HDF5 input mass data file.");
    p->addArgument(sim_dkey);
    Argument<std::string>* blob_file = new Argument<std::string>("blob_file"); // input blob json file
    blob_file->setHelp("Blobs json file.");
    p->addArgument(blob_file);
    Argument<std::string>* rad_file = new Argument<std::string>("rad_file");  // input rad_file
    rad_file->setRequired(false);
    rad_file->setHelp("Input HDF5 spectrum data file.");
    p->addArgument(rad_file);                                     

    // Simulated system parameters
    p->addArgument(new Argument<double>("m_primary", 0.6));
    p->addArgument(new Argument<double>("r_in", 0.01));
    p->addArgument(new Argument<double>("r_out", 2.0));

    // inner / outer mass influx
    p->addArgument(new Argument<double>("Q", 1e17));        // global disc mass influx
    p->addArgument(new Argument<double>("q", 0.5));         // local model mass influx

    // Radiation wavelength specification (range, step)
    p->addArgument(new Argument<double>("wl_low", 2e-5));
    p->addArgument(new Argument<double>("wl_high", 9e-5));
    p->addArgument(new Argument<double>("wl_step", 1e-7));

    Argument<double>* temp_in = new Argument<double>("temp_in", 1e5); // central object temperature
    temp_in->setHelp("Central object surface temperature.");
    p->addArgument(temp_in);

    Argument<double>* temp_atm = new Argument<double>("temp_atm", 1e5); // accretion disk atmosphere temperature
    temp_in->setHelp("Acreetion disk atmosphere temperature.");
    p->addArgument(temp_atm);

    // Simulation model ode stepper
    Argument<std::string>* stepper = new Argument<std::string>("stepper", "fehlberg78");
    stepper->setHelp("ODE stepper method.");
    p->addArgument(stepper);
}

int main(int argc, char **argv) {
    // MPI init
    int rank, size, n_workers; // mpi process rank, num of mpi processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    n_workers = size-1;

    // create and init argParser
    ArgumentParser* p = new ArgumentParser();
    ArgumentParserInit(p);

    // parse cli passed arguments
    if (!p->parse(argc, argv)) {
        if (rank == cs::mpi::MASTER) {
            std::cout << *p; // prints out help msg.
        }
        return 0;
    }

    // Always start on new line
    if (rank == cs::mpi::MASTER) {
        std::cout << std::endl;
    }

    // Mass distribution
    if (p->b("sim")) {
        sim(rank, n_workers, p);
    }

    // Radiative output
    if (p->b("rad")) {
        rad(rank, n_workers, p);
    }

    // 'Observation' and filtering
    if (p->b("obs")) {
        obs(rank, n_workers, p);
    }
    
    // terminate mpi execution env
    MPI_Finalize();

    // Destroy Arg. parser object
    delete p;

    return 0;
}
