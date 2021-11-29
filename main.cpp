#include <iostream>
#include <mpi.h>

#include "argparse-cpp/Argument.h"
#include "argparse-cpp/ArgumentParser.h"

#include "functions.cpp"
#include "Simulation.cpp"
#include "Radiation.cpp"

// simulation task
void sim(int rank, int n_workers, ArgumentParser* p) {
    // Comms dimensions 
    size_t n_jobs                   = p->i("idim") * p->i("jdim"); // number of jobs (cells)
    std::vector<size_t> comm_dim    = {n_jobs / n_workers, 11};
    while (comm_dim[0] * n_workers < n_jobs) { // int. rounding correction
        comm_dim[0]++;
    }
    
    // Process specific task
    if (rank == MASTER) {
        Simulation::master(comm_dim, n_workers, p);
    } else {
        Simulation::slave(comm_dim);
    }
}

// radiation task
void rad(int rank, int n_workers, ArgumentParser* p) {
    // Comms dimensions
    std::vector<size_t> comm_dim = {10, 10};

    // Process specific task
    if (rank == MASTER) {
        Radiation::master(comm_dim, n_workers, p);
    } else {
        Radiation::slave(comm_dim);
    }
}

int main(int argc, char **argv) {
    // create argParser
    ArgumentParser* p = new ArgumentParser();

    // add integer args.
    p->addArgument( new Argument<int>("step_n", 5e5));          // number of simulation steps
    p->addArgument( new Argument<int>("step_first"));           // first step in range
    p->addArgument( new Argument<int>("step_last"));            // last step in range
    p->addArgument( new Argument<int>("idim"));                 // number of layers(rings)
    p->addArgument( new Argument<int>("jdim"));                 // number of cell is each layer

    // add string args.
    p->addArgument( new Argument<std::string>("task", "sim"));  // select specific task (sim, rad, ...)
    p->addArgument( new Argument<std::string>("outdir"));       // data output directory
    p->addArgument( new Argument<std::string>("mass_file"));    // input mass data file
    p->addArgument( new Argument<std::string>("mass_dkey"));    // initial data key of input mass file
    p->addArgument( new Argument<std::string>("blob_file"));    // input blob json file
    p->addArgument( new Argument<std::string>("drain_file"));   // input drain data file

    // add double args.
    p->addArgument( new Argument<double>("dx", 0.01));
    p->addArgument( new Argument<double>("x", 0.0));
    p->addArgument( new Argument<double>("m_primary", 1.0));
    p->addArgument( new Argument<double>("r_in", 6.96e8));
    p->addArgument( new Argument<double>("r_out", 50.0 * 6.96e8));


    // Command line args. parser
    if (!p->parse(argc, argv)) {
        std::cout << p; // prints out help msg.
        return 0;
    }

    // MPI init
    int rank, size, n_workers; // mpi process rank, num of mpi processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    n_workers = size-1;

    // task specific call
    if (p->s("task") == "sim") {
        sim(rank, n_workers, p);
    } else if (p->s("task") == "rad") {
        rad(rank, n_workers, p);
    }
    
    // terminate mpi execution env
    MPI_Finalize();

    // Destroy Arg. parser object
    delete p;

    return 0;
}
