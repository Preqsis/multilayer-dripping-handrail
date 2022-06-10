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
    // ---------------
    // Run parameters
    // ---------------

    /**
     * Verbosity
     * default: false
     */
    Argument<bool>* verbose = new Argument<bool>("verbose", false);
    verbose->setShorthand("v");
    verbose->setHelp("Enables verbose mode.");
    p->addArgument(verbose);

    /**
     * Run simulation
     * default: false
     */
    Argument<bool>* action_sim = new Argument<bool>("sim", false);
    action_sim->setHelp("Run simulation task.");
    p->addArgument(action_sim);

    /**
     * Run radiation
     * default: false
     */
    Argument<bool>* action_rad = new Argument<bool>("rad", false);
    action_rad->setHelp("Run radiation task.");
    p->addArgument(action_rad);

    /**
     * Run synthetic observation
     * default: false
     */
    Argument<bool>* action_obs = new Argument<bool>("obs", false);
    action_obs->setHelp("Run mass observation task.");
    p->addArgument(action_obs);


    // ----------------------
    // Simulation parameters
    // ----------------------

    p->addArgument(new Argument<int>("n", 5e5, "Number of simulation steps.")); // number of simulation steps
    p->addArgument(new Argument<int>("step_first", 0)); // first step in range
    p->addArgument(new Argument<int>("step_last", 0)); // last step in range
    
    /**
     * Number of layers(rings)
     */
    Argument<int>* idim = new Argument<int>("idim");
    idim->setRequired(true);
    idim->setShorthand("I");
    idim->setHelp("Number of layers(aka. rings).");
    p->addArgument(idim);
    
    /**
     * Number of cell in one layer(ring)
     */
    Argument<int>* jdim = new Argument<int>("jdim");
    jdim->setRequired(true);
    jdim->setShorthand("J");
    jdim->setHelp("Number of cells in each ring.");
    p->addArgument(jdim);
    
    /**
     * ODE stepper
     * - used internaly in MSMM model
     */
    Argument<std::string>* stepper = new Argument<std::string>("stepper", "fehlberg78");
    stepper->setHelp("ODE stepper method.");
    p->addArgument(stepper);

    // --------------
    // Input / Output
    // --------------

    /**
     * Data output directory
     */
    Argument<std::string>* outdir = new Argument<std::string>("outdir");
    outdir->setRequired(true);
    outdir->setShorthand("o");
    outdir->setHelp("Output data directory.");
    p->addArgument(outdir);
    
    /**
     * External initialization sim. file
     */
    Argument<std::string>* init_file = new Argument<std::string>("init_file");
    init_file->setRequired(false);
    init_file->setHelp("Input HDF5 mass init/simulation file.");
    p->addArgument(init_file);                                     
    
    /**
     * Data key in external initialization sim. file
     */
    Argument<std::string>* init_dkey = new Argument<std::string>("init_dkey");
    init_dkey->setHelp("Initial data key in HDF5 input mass init/simulation file.");
    p->addArgument(init_dkey);

    /**
     * Blob json file
     */
    Argument<std::string>* blob_file = new Argument<std::string>("blob_file");
    blob_file->setHelp("Blobs json file.");
    p->addArgument(blob_file);
    
    /**
     * Specific rad. file (--obs modul)
     */
    Argument<std::string>* rad_file = new Argument<std::string>("rad_file");
    rad_file->setRequired(false);
    rad_file->setHelp("Input HDF5 spectrum data file.");
    p->addArgument(rad_file);
    
    /**
     * First simulation step to save
     * default: all
     */
    Argument<int>* save_start = new Argument<int>("save_start");
    save_start->setHelp("First step to save.");
    p->addArgument(save_start);
    
    /**
     * Last simulation step to save
     * default: all
     */
    Argument<int>* save_end = new Argument<int>("save_end");
    save_end->setHelp("Last step to save.");
    p->addArgument(save_end);

    // ------------------
    // System parameters
    // ------------------
    
    /**
     * Primary object mass
     */
    Argument<double>* m_primary = new Argument<double>("m_primary", 0.6);
    m_primary->setShorthand("M");
    m_primary->setHelp("Primary object mass.");
    p->addArgument(m_primary);
    
    /**
     * Influx gas gemteprature
     */
    Argument<double>* T_flow = new Argument<double>("T_flow", 4500);
    T_flow->setHelp("Influx temperature.");
    p->addArgument(T_flow);
    
    p->addArgument(new Argument<double>("r_in", 0.01));
    p->addArgument(new Argument<double>("r_out", 2.0));

    p->addArgument(new Argument<double>("Q", 1e14));        // global disc mass influx
    p->addArgument(new Argument<double>("q", 0.9));         // local model mass influx
    
    p->addArgument(new Argument<double>("wl_low", 1e-5));
    p->addArgument(new Argument<double>("wl_high", 9e-5));
    p->addArgument(new Argument<double>("wl_step", 1e-7));
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

    // Simulation
    // - mass flow simulation
    if (p->b("sim")) {
        sim(rank, n_workers, p);
    }

    // Radiation
    // - generates spectrum for all cells
    if (p->b("rad")) {
        rad(rank, n_workers, p);
    }

    // Observation
    // - synthetic filtred light curve 
    if (p->b("obs")) {
        obs(rank, n_workers, p);
    }
    
    // terminate mpi execution env
    MPI_Finalize();

    // Destroy Arg. parser object
    delete p;

    return 0;
}
