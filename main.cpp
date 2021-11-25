#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cstdlib>
#include <filesystem>

#include <H5Attribute.hpp>
#include <H5DataSet.hpp>
#include <H5DataSpace.hpp>
#include <H5File.hpp>
namespace H5 = HighFive;

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include "Argument.h"
#include "ArgumentParser.h"

#include "MSMM.h"
#include "BlobScheduler.h"
#include "Distributor.h"

typedef boost::array<double, 2> state_type;
typedef long unsigned int lint;
typedef unsigned int uint;

const double M_SUN  = 1.9891e30;
const double G      = 6.6743e-11;

// Compute flags
const int MASTER    = 0;
const int COMPUTE   = 1;
const int STOP      = 2;

// alloc 2d double pointer array
// dims specified by individual ints
double** alloc_2D_double(int idim, int jdim) {
    double** arr;
    arr     = (double**)calloc(idim, sizeof(double*));
    arr[0]  = (double*)calloc(idim*jdim, sizeof(double));
    for (int i = 0; i < idim; i++) {
        arr[i] = arr[0] + i * jdim;
    }
    return arr;
}

double*** alloc_3D_double(size_t idim, size_t jdim, size_t kdim) {
    double ***aaa, **aa, *a;
    uint i;

    aaa     = (double***) calloc(idim, sizeof(double**));
    aa      = (double**) calloc(idim * jdim, sizeof(double*));
    aaa[0]  = aa;

    for (i = 1; i < idim; i++) {
        aaa[i] = aaa[i - 1] + jdim;
    }
    
    a       = (double *) calloc(idim * jdim * kdim, sizeof(double));
    aa[0]   = a;
    for (i = 1; i < idim * jdim; i++) {
        aa[i] = aa[i - 1] + kdim;
    }
    return aaa;
}

double*** alloc_3D_double(std::vector<size_t> dim) {
    return alloc_3D_double(dim[0], dim[1], dim[2]);
}

// alloc 2d double pointer array
// dims specified by vector
double** alloc_2D_double(std::vector<size_t> dim) {
    return alloc_2D_double(dim[0], dim[1]);
}

// 'Empty' grid initalization
void grid_init(uint idim, uint jdim, double r_in, double r_out, double dx, double** data) {
    uint k; 
    for (uint i=0; i<idim; i++) { // rings
        for (uint j=0; j<jdim; j++){ // cells
            k = i * jdim + j; // serialized communication data coordinate
            data[k][0] = i; // ring coordinate
            data[k][1] = j; // cell coordinate
            data[k][2] = 0.0; // t
            data[k][3] = 2.0; // z
            data[k][4] = 0.0; // v
            data[k][5] = 0.5; // m
            data[k][6] = 0.0; // dm
            data[k][7] = std::pow(r_out, 2.0) / std::pow(r_in + (idim - i - 1) * (r_out - r_in) / (idim - 1), 2.0); // r (ring specific radius)
            data[k][8] = dx; // dx (inner MSMM parameter)
            data[k][9] = 2.0 * M_PI * ((double) j) / ((double) jdim); // azimuth (cell specific rotation angle)
            data[k][10] = 1.0; // compute flag (0.0 --> do not compute)
        }
    }
}

void writeDataSet(H5::File* file, double** data, std::vector<size_t> dim, std::string key) {
    if (!file->exist(key)) {
        H5::DataSet* ds = new H5::DataSet(file->createDataSet<double>(key, H5::DataSpace(dim)));
        ds->write((double**) data[0]);
        delete ds;
    }
}

// Worker function
void worker(double** data, MSMM2* model, std::vector<size_t> cdim) {
    // solver / stepper
    //boost::numeric::odeint::runge_kutta4<state_type> stepper;
    //boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;

    // pres vsechny zadane 'bunky'
    for (uint i = 0; i < cdim[0]; i++) {
        // bunku nepocitat!
        if (data[i][10] == 0.0) {continue;}

        // bunka s nulovou hmotnosti
        // rovnice nedavaji smysl
        if (data[i][5] == 0.0) {continue;}

        // pocatectni podminky
        double x        = data[i][2];
        state_type y    = {data[i][3], data[i][4]};

        // velikost kroku
        double dx       = data[i][8];

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
}

// Function for MPI slave processes
void MPI_slave(std::vector<size_t> cdim) {
    /** MPI status flag **/
    MPI_Status status;

    /** Comms matrix allocation */
    double** data_recv_slave = alloc_2D_double(cdim);
    double** data_send_slave = alloc_2D_double(cdim);
    
    // mass-spring model
    MSMM2* model = new MSMM2();

    /** Recieve until STOP */
    while (true) {
        /** Recieve data */
        MPI_Recv(&data_recv_slave[0][0], cdim[0]*cdim[1], MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        /** If STOP --> end slave computation */
        if (status.MPI_TAG == STOP) {
            free(model);
            break;
        }

        /** Run worker on recieved data */
        worker(data_recv_slave, model, cdim);

        /** Send results to master */
        MPI_Send(&data_recv_slave[0][0], cdim[0]*cdim[1], MPI_DOUBLE, 0, COMPUTE, MPI_COMM_WORLD);
    }
}

// Function for MPI master process
void MPI_master(std::vector<size_t> cdim, int n_workers, ArgumentParser* p) {
    double m_primary, dx, x, r_in, r_out; // system params
    // from CLAs to easilly accesible vars
    m_primary       = p->d("m_primary") * M_SUN;
    dx              = p->d("dx");
    x               = p->d("x");
    r_in            = p->d("r_in");
    r_out           = p->d("r_out");




    MPI_Status status;
    lint s                      = (lint)p->i("s");
    lint idim                   = (lint)p->i("idim");
    lint jdim                   = (lint)p->i("jdim");

    // "real" sim dimensions
    // {rings, cells}
    std::vector<size_t> dim = {idim, jdim};

    // grid allocation
    // empty vs. external data
    std::vector<size_t> grid_dim    = {idim*jdim, cdim[1]};
    double** grid                   = alloc_2D_double(grid_dim);
    if (p->isSet("mass_file") && p->isSet("mass_dkey")) { // sim start of external data
        H5::File* mass_infile = new H5::File(p->s("mass_file"), H5::File::ReadOnly);
        mass_infile->getDataSet(p->s("mass_dkey")).read((double**) grid[0]);
    } else { // sim starts by setting parameters
        grid_init(idim, jdim, r_in, r_out, dx, grid);
    }

    // if not exist create outdir
    if (!std::filesystem::is_directory(p->s("outdir"))) {
        std::filesystem::create_directory(p->s("outdir"));
    }

    // mass output
    H5::File* mass_file;
    mass_file = new H5::File(p->s("outdir") + "/mass.h5", H5::File::Overwrite);
    mass_file->createAttribute<int>("idim", H5::DataSpace::From(idim)).write(idim);
    mass_file->createAttribute<int>("jdim", H5::DataSpace::From(jdim)).write(jdim);
    writeDataSet(mass_file, grid, grid_dim, "data_0"); // initial state save

    // drain allocation
    std::vector<size_t> drain_dim   = {idim, jdim};
    double** drain                  = alloc_2D_double(drain_dim);
    for (uint i=0; i<drain_dim[0]; i++)
        for (uint j=0; j<drain_dim[1]; j++)
            drain[i][j] = 0.0;

    // drain output
    H5::File* drain_file;
    drain_file = new H5::File(p->s("outdir") + "/drain.h5", H5::File::Overwrite);
    drain_file->createAttribute<int>("idim", H5::DataSpace::From(idim)).write(idim);
    drain_file->createAttribute<int>("jdim", H5::DataSpace::From(jdim)).write(jdim);
    writeDataSet(drain_file, drain, drain_dim, "data_0"); // initial state save

    /**
        * 'Rotation profile' 
        * - specifies angular velocity for each ring relative to ring i == 0
        */
    std::vector<double> profile;
    double T_out = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r_out, 3.0)) / (G * m_primary) , 1.0/2.0);
    double r, T;
    for (uint i=0; i<idim; i++) {
        r = r_in + (idim - i - 1) * (r_out - r_in) / (idim - 1); // ring specific radius
        T = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r, 3.0)) / (G * m_primary) , 1.0/2.0); // ring spe
        profile.push_back(std::round((std::pow(10.0, 4.0) * T_out * 2.0 * M_PI) / (T * ((double)jdim))) / std::pow(10.0, 4.0));
    }

    /**
        * 'Distributor' object
        * - handles mass redistribution (free flow and 'dripping')
        */
    Distributor* dist = new Distributor(dim, cdim, drain);
    if (p->isSet("blobfile")) {
        dist->setBlobScheduler(new BlobScheduler(dim, grid, p->s("blobfile")));
    }
    dist->setRotationProfile(profile);

    /* Create jobs set for each worker */
    uint k;
    std::vector<std::vector<uint>> jobs = {};
    for (uint w = 0; w < n_workers; w++) {
        jobs.push_back({});
        for (uint i = 0; i < cdim[0]; i++) {
            k = w * cdim[0] + i;
            jobs[w].push_back(k);
        }
    }

    /** Communcdim[0]ation matrix allocation */
    double** data_send_master = alloc_2D_double(cdim);
    double** data_recv_master = alloc_2D_double(cdim);

    /** Percentage info msg. variables */
    uint percent = 0;
    uint one_percent = s / 100;


    /** Compute all simulation steps */
    for (uint step=1; step <= s; step++) {
        /** Sort, mark (compute flag) and send jobs to specifcdim[0] workers */
        uint w;
        for (uint slave=1; slave<=n_workers; slave++) {
            w = slave-1;
            for (uint i=0; i<jobs[w].size(); i++) {
                k = jobs[w][i];

                if (k >= idim*jdim) { // outside of grid range --> compute = 0.0
                    data_send_master[i][cdim[1]-1] = 0.0;
                    continue;
                }

                /** actual grid / cell data to communcdim[0]ation matrix */
                for (uint l=0; l<cdim[1]; l++) {
                    data_send_master[i][l] = grid[k][l];
                }
                data_send_master[i][cdim[1]-1] = 1.0; // set compute flag --> compute = 1.0 (to be computed)
            }
            MPI_Send(&data_send_master[0][0], cdim[0]*cdim[1], MPI_DOUBLE, slave, COMPUTE, MPI_COMM_WORLD); // send data to slave
        }

        /** Recieve and sort data back to grid */
        for (uint slave=1; slave <= n_workers; slave++) {
            MPI_Recv(&data_recv_master[0][0], cdim[0]*cdim[1], MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // recv. data

            for (uint i=0; i<cdim[0]; i++) { // back to grid
                if (data_recv_master[i][10] == 0.0) {continue;} // compute == 1.0 only
                w = slave-1;
                k = w * cdim[0] + i;
                for (uint l=0; l<cdim[1]; l++) {
                    grid[k][l] = data_recv_master[i][l];
                }
            }
        }

        // write mass
        writeDataSet(mass_file, grid, grid_dim, "data_" + std::to_string(step));

        /** Run distribution handler object on grid */
        dist->step(grid, step);

        // write drain
        writeDataSet(drain_file, drain, drain_dim, "data_" + std::to_string(step));

        /** Print percentage msg. */
        if (s >= 100) {
            if (percent != step / one_percent) {
                percent = step / one_percent;
                std::cout << "Running ... " << percent << "%\t\r" << std::flush;
            }
        }
    }

    /** Stop all workers (slave) by sending STOP flag (and dummy data) */
    //double** dummy = alloc_2D_double(cdim);
    for (uint i=1; i <= n_workers; i++) {
        MPI_Send(&data_send_master[0][0], cdim[0]*cdim[1], MPI_DOUBLE, i, STOP, MPI_COMM_WORLD); 
    }
}

// Creates argument parser instance 
ArgumentParser* createArgumentParser() {
    ArgumentParser* p = new ArgumentParser();

    // add integer args.
    p->addArgument( new Argument<int>("s", 5e5));
    p->addArgument( new Argument<int>("idim"));
    p->addArgument( new Argument<int>("jdim"));

    // add double args.
    p->addArgument( new Argument<double>("dx", 0.01));
    p->addArgument( new Argument<double>("x", 0.0));
    p->addArgument( new Argument<double>("m_primary", 1.0));
    p->addArgument( new Argument<double>("r_in", 6.96e8));
    p->addArgument( new Argument<double>("r_out", 50.0 * 6.96e8));

    // add string args.
    p->addArgument( new Argument<std::string>("outdir"));
    p->addArgument( new Argument<std::string>("mass_file"));
    p->addArgument( new Argument<std::string>("mass_dkey"));
    p->addArgument( new Argument<std::string>("blob_file"));

    return p;
}

int main(int argc, char **argv) {
    // Command line args. parser
    ArgumentParser* p = createArgumentParser();
    if (!p->parse(argc, argv)) {
        std::cout << p; // prints out help msg.
        return 0;
    }

    // MPI init
    int p_rank; // mpi process rank
    int size; // num of mpi processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Comms dimensions 
    size_t n_jobs               = p->i("idim") * p->i("jdim"); // number of jobs (cells)
    size_t n_workers            = size-1; // number of MPI worker processes 
    std::vector<size_t> cdim    = {n_jobs / n_workers, 11};
    while (cdim[0] * n_workers < n_jobs) { // int. rounding correction
        cdim[0]++;
    }
    
    // Process specific task
    if (p_rank == MASTER) {
        MPI_master(cdim, n_workers, p);
    } else {
        MPI_slave(cdim);
    }

    // terminate mpi execution enviroment
    MPI_Finalize();

    // Destroy Arg. parser object
    delete p;

    return 0;
}
