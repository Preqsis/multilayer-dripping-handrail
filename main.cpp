#include <iostream>
#include <mpi.h>
#include <math.h>

#include "argparse-cpp/Argument.h"
#include "argparse-cpp/ArgumentParser.h"

#include "functions.cpp"
namespace fn = Functions;

#include "constants.cpp"
namespace cs = Constants;

#include "Simulation.cpp"
#include "Radiation.cpp"

// simulation task
void sim(int rank, int n_workers, ArgumentParser* p) {
    // Comms dimensions 
    size_t n_jobs               = p->i("idim") * p->i("jdim"); // number of jobs (cells)
    std::vector<size_t> cdim    = {(size_t)std::ceil((double)n_jobs / (double)n_workers), 13};
    
    // Process specific call
    if (rank == MASTER) {
        Simulation::master(cdim, n_workers, p);
    } else {
        Simulation::slave(cdim, p);
    }
}

// radiation task
void rad(int rank, int n_workers, ArgumentParser* p) {
    // Comms dimensions
    std::vector<size_t> dim_mass    = {(size_t)p->i("idim"), (size_t)p->i("jdim"), 13};
    std::vector<size_t> dim_spec    = {dim_mass[0], dim_mass[1], (size_t)((p->d("lam_high") - p->d("lam_low")) / p->d("lam_step") + 1), 2};

    // Process specific call
    if (rank == MASTER) {
        Radiation::master(dim_mass, dim_spec, n_workers, p);
    } else {
        Radiation::slave(dim_mass, dim_spec, p);
    }
}

// observation task
void obs(int rank, int n_workers, ArgumentParser* p) {

}

int main(int argc, char **argv) {
    // MPI init
    int rank, size, n_workers; // mpi process rank, num of mpi processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    n_workers = size-1;

    // create argParser
    ArgumentParser* p = new ArgumentParser();

    // taks selection and verbosity 
    p->addArgument(new Argument<bool>("v", false));     // verbosity
    p->addArgument(new Argument<bool>("sim", false));   // run mass distribution sim
    p->addArgument(new Argument<bool>("rad", false));   // run radiation output computation
    p->addArgument(new Argument<bool>("obs", false));   // run observation and filtering

    // Steps (number of steps, range, etc.)
    p->addArgument( new Argument<int>("step_n", 5e5));  // number of simulation steps
    p->addArgument( new Argument<int>("step_first"));   // first step in range
    p->addArgument( new Argument<int>("step_last"));    // last step in range
    p->addArgument( new Argument<int>("step_wfirst"));  // first step to save
    p->addArgument( new Argument<int>("step_wlast"));   // first step to save

    // Disk dimensions
    Argument<int>* idim = new Argument<int>("idim");    // number of layers
    idim->setRequired(true);
    p->addArgument(idim);
    Argument<int>* jdim = new Argument<int>("jdim");    // number of cells in each layer
    jdim->setRequired(true);
    p->addArgument(jdim);

    // Input / output defs
    Argument<std::string>* outdir = new Argument<std::string>("outdir");        // data output directory
    outdir->setRequired(true);
    p->addArgument(outdir);                                     
    Argument<std::string>* mass_file = new Argument<std::string>("mass_file");  // input mass_file
    mass_file->setRequired(false);
    p->addArgument(mass_file);                                     
    p->addArgument( new Argument<std::string>("mass_dkey"));                    // initial data key of input mass file
    p->addArgument( new Argument<std::string>("blob_file"));                    // input blob json file

    // Simulated system parameters
    p->addArgument( new Argument<double>("m_primary", 0.8));
    p->addArgument( new Argument<double>("r_in", 5e8));
    p->addArgument( new Argument<double>("r_out", 50.0 * 5e8));

    // inner / outer mass influx
    p->addArgument(new Argument<double>("Q", 1e17));        // global disc mass influx
    p->addArgument(new Argument<double>("q", 0.5));         // local model mass influx

    // Radiation wavelength specification (range, step)
    p->addArgument( new Argument<double>("lam_low", 1e-5));
    p->addArgument( new Argument<double>("lam_high", 9e-5));
    p->addArgument( new Argument<double>("lam_step", 1e-7));

    // Command line args. parser
    if (!p->parse(argc, argv)) {
        if (rank == MASTER) {
            std::cout << *p; // prints out help msg.
        }
        return 0;
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
