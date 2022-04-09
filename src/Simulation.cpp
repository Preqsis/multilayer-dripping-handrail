#include <iostream>
#include <math.h>
#include <mpi.h>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include "argparse-cpp/ArgumentParser.h"

#include "MSMM.h"

#include "Distributor.hpp"

#include "Functions.hpp"
namespace fn = Functions;

#include "Constants.hpp"
namespace cs = Constants;

#include "Simulation.hpp"

// dataset write check function
bool Simulation::writable(bool wf, bool wl, uint s, uint sf, uint sl) {
    return (!wf && !wl) || (wf && !wl && sf <= s) || (!wf && wl && s <= sl) || (wf && wl && sf <= s && s <= sl);
}

// 'Empty' simulation grid initalization
void Simulation::grid_init(double*** data, std::vector<size_t> dim, ArgumentParser* p) {
    double r_in     = p->d("r_in");
    double r_out    = p->d("r_out");
    //double T_in     = std::pow((3.0 * cs::G * p->d("m_primary") * cs::m_sun * 1e17 / (8.0 * M_PI * cs::k * std::pow(r_in, 3.0))), 0.25);
    double T_in     = p->d("temp_in");

    for (uint i=0; i< dim[0]; i++) { // rings
        for (uint j=0; j<dim[1]; j++){ // cells
            data[i][j][0]   = i; // ring coordinate
            data[i][j][1]   = j; // cell coordinate
            data[i][j][2]   = 0.0; // t
            data[i][j][3]   = 0.0; // z
            data[i][j][4]   = 0.0; // v
            data[i][j][5]   = 0.0; // m
            data[i][j][6]   = 0.0; // dm
            data[i][j][7]   = r_out - i * (r_out - r_in) / ((double) dim[0]); // cell specific radius
            data[i][j][8]   = 2.0 * M_PI * ((double) j) / ((double) dim[1]); // azimuth (cell specific rotation angle)
            data[i][j][9]   = 0.0; // drain
            data[i][j][10]  = T_in * std::pow((data[i][j][7] / r_in ), -0.75); // radius dependent disk body temperature
            data[i][j][11]  = std::pow(r_out, 2.0) / std::pow(((double)dim[0] - i - 1) * (r_out - r_in) / ((double)dim[0] - 1) , 2.0); // local 'g'
            //data[i][j][11]  = 1.0; // local 'g'
            data[i][j][12]  = 1.0; // compute flag (0.0 --> do not compute)
        }
    }
}

void Simulation::terminate(std::vector<size_t> dim, int n_workers) {
    double** data = fn::alloc_2D_double(dim);
    for (int slave = 1; slave <= n_workers; slave++) {
        MPI_Send(&data[0][0], dim[0]*dim[1], MPI_DOUBLE, slave, cs::mpi::STOP, MPI_COMM_WORLD); 
    }
}

// Function for "sim" MPI slave processes
void Simulation::slave(std::vector<size_t> dim, ArgumentParser* p) {
    /** MPI status flag **/
    MPI_Status status;

    // dx --> DORESIT!!!
    double dx = 0.1; 

    /** Comms matrix allocation */
    double** data   = fn::alloc_2D_double(dim);
    size_t len_data = dim[0] * dim[1];
    
    // mass-spring model
    MSMM2* model = new MSMM2();
        
    // solver / stepper
    boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;
    if (p->s("stepper") == "rk4") {
        boost::numeric::odeint::runge_kutta4<state_type> stepper;
    } else if (p->s("stepper") == "dopri5") {
        boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    }

    /** Recieve until STOP */
    while (true) {
        /** Recieve data */
        MPI_Recv(&data[0][0], len_data, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        /** If STOP --> end slave computation */
        if (status.MPI_TAG == cs::mpi::STOP) {
            break;
        }

        // pres vsechny zadane 'bunky'
        for (uint i = 0; i < dim[0]; i++) {
            // bunku nepocitat!
            if (data[i][dim[1]-1] == 0.0) {continue;}

            // bunka s nulovou hmotnosti
            // rovnice nedavaji smysl
            if (data[i][5] == 0.0) {continue;}

            // pocatectni podminky
            double      x = data[i][2];
            state_type  y = {data[i][3], data[i][4]};

            // parametry do modelu
            model->set_m(data[i][5]);
            model->set_dm(data[i][6] / dx); // pritok urceny vstupem, podeleno krokem aby davalo smysl pro rce.

            model->set_g(data[i][11]); // local 'g'

            // krok
            stepper.do_step(*model, y, x, dx);

            // vysledky zpet do data
            data[i][2] += dx; // inkrementace casu
            data[i][3] = y[0]; // z
            data[i][4] = y[1]; // v
        }


        /** Send results to master */
        MPI_Send(&data[0][0], len_data, MPI_DOUBLE, 0, cs::mpi::COMPUTE, MPI_COMM_WORLD);
    }

    // cleanup
    delete model;
}

// Function for "sim" MPI master process
void Simulation::master(std::vector<size_t> comm_dim, int n_workers, ArgumentParser* p) {
    double m_primary, r_in, r_out; // system params
    // from CLAs to easilly accesible vars
    m_primary       = p->d("m_primary") * cs::m_sun;
    r_in            = p->d("r_in");
    r_out           = p->d("r_out");
    
    // MPI status flag holder
    MPI_Status status;

    // verbosity?
    bool v = p->b("verbose");
    
    // number of simulation steps
    uint n = p->i("n");
    int one_percent = n / 100;
    

    // write checking vars
    bool wf = p->isSet("save_start");
    uint sf = wf ? p->i("save_start") : 0;
    bool wl = p->isSet("save_end");
    uint sl = wl ? p->i("save_end") : 0;

    // "real" sim dimensions
    // {rings, cells, cell_data}
    std::vector<size_t> dim = {(size_t)p->i("idim"), (size_t)p->i("jdim"), comm_dim[1]};

    // grid allocation
    double*** grid = fn::alloc_3D_double(dim);
    if (p->isSet("sim_file") && p->isSet("sim_dkey")) { // sim start of external data
        H5::File* sim_infile = new H5::File(p->s("sim_file"), H5::File::ReadOnly);
        sim_infile->getDataSet(p->s("sim_dkey")).read((double***) grid[0][0]);
        delete sim_infile;
    } else { // sim starts by setting parameters
        grid_init(grid, dim, p);
    }

    // if not exist create outdir
    fn::mkdir(p->s("outdir"));

    // mass output
    H5::File* sim_file;
    sim_file = new H5::File(p->s("outdir") + "/sim.h5", H5::File::Overwrite);
    sim_file->createAttribute<int>("idim", H5::DataSpace::From(dim[0])).write(dim[0]);
    sim_file->createAttribute<int>("jdim", H5::DataSpace::From(dim[1])).write(dim[1]);

    // initial state save and write check
    if (writable(wf, wl, 0, sf, sl)) {
        fn::writeDataSet(sim_file, grid, dim, "data_0"); 
    }

    // Rotation profile
    // - specifies angular velocity for each ring relative to i == 0
    std::vector<double> profile;
    double T_out = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r_out, 3.0)) / (cs::G * m_primary) , 1.0/2.0);
    double r, T;
    for (uint i=0; i<dim[0]; i++) {
        r = r_in + (dim[0] - i - 1) * (r_out - r_in) / (dim[0] - 1); // ring specific radius
        T = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r, 3.0)) / (cs::G * m_primary) , 1.0/2.0); // ring spe
        profile.push_back(std::round((std::pow(10.0, 4.0) * T_out * 2.0 * M_PI) / (T * ((double)dim[1]))) / std::pow(10.0, 4.0));
    }

    // Distributor object
    // handles mass redistribution
    Distributor* dst;
    if (p->isSet("blob_file")) 
        dst = new Distributor(grid, dim, p->s("blob_file"));
    else 
        dst = new Distributor(grid, dim);
    dst->setParams(p->d("m_primary"), p->d("r_in"), p->d("r_out"), p->d("Q"), p->d("q"));
    dst->setRotationProfile(profile);
    
    double dt = dst->get_dt();
    sim_file->createAttribute<double>("dt", H5::DataSpace::From(dt)).write(dt);

    // Communication data 'matrix' allocation
    double** data   = fn::alloc_2D_double(comm_dim);
    size_t len_data = comm_dim[0] * comm_dim[1];

    // Compute all simulation steps
    uint i, j, k, c;
    int slave;
    for (uint s = 1; s < n; s++) {
        // Sort, mark (compute flag) and send jobs to specific slave processes
        i = 0;
        j = 0;
        for (slave = 1; slave <= n_workers; slave++) {
            for (c = 0; c < comm_dim[0]; c++) {
                if (i < dim[0]) {
                    for (k = 0; k < comm_dim[1]; k++) {
                        data[c][k] = grid[i][j][k];
                    }
                }

                // compute flag
                data[c][comm_dim[1]-1] = (i < dim[0]) ? 1.0 : 0.0;

                // cell indexes
                j = (j < dim[1]-1) ? j + 1 : 0;
                i = (j == 0) ? i + 1 : i;
            }
            MPI_Send(&data[0][0], len_data, MPI_DOUBLE, slave, cs::mpi::COMPUTE, MPI_COMM_WORLD); // send data to slave
        }

        // Recieve and sort data back to grid
        for (slave = 1; slave <= n_workers; slave++) {
            MPI_Recv(&data[0][0], len_data, MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // recv. data

            for (c = 0; c < comm_dim[0]; c++) { // back to grid
                if (data[c][comm_dim[1]-1] == 0.0) {continue;} // compute == 1.0 only

                i = (uint) data[c][0];
                j = (uint) data[c][1];

                for (k = 0; k < comm_dim[1]; k++) {
                    grid[i][j][k] = data[c][k];
                }
            }
        }
        
        // Run distribution handler on grid
        dst->run(s);

        // write mass
        if (writable(wf, wl, s, sf, sl)) {
            fn::writeDataSet(sim_file, grid, dim, "data_" + std::to_string(s));
        }

        // info msg.
        if (v) {
            if (one_percent == 0 || s % one_percent == 0) {
                std::cout << "Mass distribution ... " << (int) (100.0 * s / (n)) << "%\t\r" << std::flush;
            }
        }
    }

    // Stop all slave processes
    terminate(comm_dim, n_workers);

    // cleanup
    delete sim_file;
    delete dst;

    // info msg.
    if (v) {
        std::cout << std::endl << "Done ..." << std::endl;
    }
}

