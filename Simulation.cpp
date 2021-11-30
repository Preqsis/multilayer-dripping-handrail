#ifndef SIMULATION_CPP
#define SIMULATION_CPP

#include <iostream>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include "MSMM.h"
#include "BlobScheduler.h"
#include "Distributor.h"

#include "functions.cpp"
namespace fn = Functions;

namespace Simulation {

// dataset write check function
bool writable(bool wf, bool wl, uint s, uint sf, uint sl) {
    return (!wf && !wl) || (wf && !wl && sf <= s) || (!wf && wl && s <= sl) || (wf && wl && sf <= s && s <= sl);
}

// 'Empty' simulation grid initalization
void grid_init(double*** data, std::vector<size_t> dim, double r_in, double r_out) {
    for (uint i=0; i< dim[0]; i++) { // rings
        for (uint j=0; j<dim[1]; j++){ // cells
            //k = i * jdim + j; // serialized communication data coordinate
            data[i][j][0] = i; // ring coordinate
            data[i][j][1] = j; // cell coordinate
            data[i][j][2] = 0.0; // t
            data[i][j][3] = 0.0; // z
            data[i][j][4] = 0.0; // v
            data[i][j][5] = 0.0; // m
            data[i][j][6] = 0.0; // dm
            data[i][j][7] = std::pow(r_out, 2.0) / std::pow(r_in + (dim[0] - i - 1) * (r_out - r_in) / (dim[0] - 1), 2.0); // r (ring specific radius)
            data[i][j][8] = 2.0 * M_PI * ((double) j) / ((double) dim[1]); // azimuth (cell specific rotation angle)
            data[i][j][9] = 0.0; // drain
            data[i][j][10] = 1000.0; // T ---> DORESIT!!! 
            data[i][j][11] = 1.0; // compute flag (0.0 --> do not compute)
        }
    }
}

// Function for "sim" MPI slave processes
void slave(std::vector<size_t> dim) {
    /** MPI status flag **/
    MPI_Status status;

    // dx --> DORESIT!!!
    double dx = 0.01; 

    /** Comms matrix allocation */
    double** data   = fn::alloc_2D_double(dim);
    size_t len_data = dim[0] * dim[1];
    
    // mass-spring model
    MSMM2* model = new MSMM2();
        
    // solver / stepper
    //boost::numeric::odeint::runge_kutta4<state_type> stepper;
    //boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;

    /** Recieve until STOP */
    while (true) {
        /** Recieve data */
        MPI_Recv(&data[0][0], len_data, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        /** If STOP --> end slave computation */
        if (status.MPI_TAG == STOP) {
            free(model);
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
            model->set_g(data[i][7]);

            // krok
            stepper.do_step(*model, y, x, dx);

            // vysledky zpet do datau
            data[i][2] += dx; // inkrementace casu
            data[i][3] = y[0]; // z
            data[i][4] = y[1]; // v
        }

        /** Send results to master */
        MPI_Send(&data[0][0], len_data, MPI_DOUBLE, 0, COMPUTE, MPI_COMM_WORLD);
    }
}

// Function for "sim" MPI master process
void master(std::vector<size_t> comm_dim, int n_workers, ArgumentParser* p) {
    double m_primary, r_in, r_out; // system params
    // from CLAs to easilly accesible vars
    m_primary       = p->d("m_primary") * M_SUN;
    r_in            = p->d("r_in");
    r_out           = p->d("r_out");

    // MPI status flag holder
    MPI_Status status;

    // number of simulation steps
    uint steps = p->i("step_n");

    // write checking vars
    bool wf = p->isSet("step_wfirst");
    uint sf = wf ? p->i("step_wfirst") : 0;
    bool wl = p->isSet("step_wlast");
    uint sl = wl ? p->i("step_wlast") : 0;

    // "real" sim dimensions
    // {rings, cells, cell_data}
    std::vector<size_t> dim = {(luint)p->i("idim"), (luint)p->i("jdim"), comm_dim[1]};

    // grid allocation
    double*** grid = fn::alloc_3D_double(dim);
    if (p->isSet("mass_file") && p->isSet("mass_dkey")) { // sim start of external data
        H5::File* mass_infile = new H5::File(p->s("mass_file"), H5::File::ReadOnly);
        mass_infile->getDataSet(p->s("mass_dkey")).read((double***) grid[0][0]);
    } else { // sim starts by setting parameters
        grid_init(grid, dim, r_in, r_out);
    }

    // if not exist create outdir
    fn::mkdir(p->s("outdir"));

    // mass output
    H5::File* mass_file;
    mass_file = new H5::File(p->s("outdir") + "/mass.h5", H5::File::Overwrite);
    mass_file->createAttribute<int>("idim", H5::DataSpace::From(dim[0])).write(dim[0]);
    mass_file->createAttribute<int>("jdim", H5::DataSpace::From(dim[1])).write(dim[1]);

    // initial state save and write check
    if (writable(wf, wl, 0, sf, sl)) {
        fn::writeDataSet(mass_file, grid, dim, "data_0"); 
    }

    // Rotation profile
    // - specifies angular velocity for each ring relative to i == 0
    std::vector<double> profile;
    double T_out = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r_out, 3.0)) / (G * m_primary) , 1.0/2.0);
    double r, T;
    for (uint i=0; i<dim[0]; i++) {
        r = r_in + (dim[0] - i - 1) * (r_out - r_in) / (dim[0] - 1); // ring specific radius
        T = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r, 3.0)) / (G * m_primary) , 1.0/2.0); // ring spe
        profile.push_back(std::round((std::pow(10.0, 4.0) * T_out * 2.0 * M_PI) / (T * ((double)dim[1]))) / std::pow(10.0, 4.0));
    }

    // Distributor object
    // handles mass redistribution
    Distributor* dist = new Distributor(grid, dim);
    if (p->isSet("blob_file"))
        dist->setBlobScheduler(new BlobScheduler(grid, dim, p->s("blob_file")));
    dist->setRotationProfile(profile);

    // Communication data 'matrix' allocation
    double** data = fn::alloc_2D_double(comm_dim);

    // Percentage info msg. variables
    uint percent        = 0;
    uint one_percent    = steps / 100;

    // Compute all simulation steps
    uint i, j, k, c;
    int slave;
    for (uint s = 1; s < steps; s++) {
        // Sort, mark (compute flag) and send jobs to specific slave processes
        i = 0;
        j = 0;
        for (slave = 1; slave <= n_workers; slave++) {
            for (c = 0; c < comm_dim[0]; c++) {
                for (k = 0; k < comm_dim[1]; k++) {
                    data[c][k] = grid[i][j][k];
                }

                // compute flag
                data[c][comm_dim[1]-1] = (i < dim[0]) ? 1.0 : 0.0;

                // cell indexes
                // rollover je z nejakeho duvodu rychlejsi ???
                i = (i < dim[0]-1) ? i + 1 : 0;
                j = (j < dim[1]-1) ? j + 1 : 0;
            }
            MPI_Send(&data[0][0], comm_dim[0]*comm_dim[1], MPI_DOUBLE, slave, COMPUTE, MPI_COMM_WORLD); // send data to slave
        }

        // Recieve and sort data back to grid
        for (slave = 1; slave <= n_workers; slave++) {
            MPI_Recv(&data[0][0], comm_dim[0]*comm_dim[1], MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // recv. data

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
        dist->run(s);

        // write mass
        if (writable(wf, wl, s, sf, sl)) {
            fn::writeDataSet(mass_file, grid, dim, "data_" + std::to_string(s));
        }

        // Print percentage msg.
        if (steps >= 100) {
            if (percent != s / one_percent) {
                percent = s / one_percent;
                std::cout << "Running ... " << percent << "%\t\r" << std::flush;
            }
        }
    }

    // Stop all workers (slave) by sending STOP flag
    for (slave = 1; slave <= n_workers; slave++) {
        MPI_Send(&data[0][0], comm_dim[0]*comm_dim[1], MPI_DOUBLE, slave, STOP, MPI_COMM_WORLD); 
    }
}

}

#endif

