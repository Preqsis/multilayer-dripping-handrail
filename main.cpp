#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cstdlib>

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

    MPI_Status status;
    int first_step, last_step; 
    H5::File* infile; 
    H5::File* outfile;
    double** grid;

    bool append                 = p->b("append");
    lint s                      = (lint)p->i("s");
    lint idim                   = (lint)p->i("idim");
    lint jdim                   = (lint)p->i("jdim");
    std::vector<size_t> dims    = {idim*jdim, cdim[1]}; // combined grid dimensions

    // decide on input type
    if (p->isSet("input")) { // sim start of external data
        // input file
        infile = new H5::File(p->s("input"), H5::File::ReadWrite);

        // decide first simulation step
        // and initial dataset
        std::string dkey;
        if (p->isSet("first_step")) {
            first_step = p->i("first_step");
            dkey = "data_" + std::to_string(first_step);
        } else {
            infile->getAttribute("last_step").read(first_step);
            dkey = "data_" + std::to_string(first_step);
            first_step += 1;
        }

        //decide last simulation step
        if (p->isSet("last_step")) {
            last_step = p->i("last_step");
        } else {
            last_step = first_step + s;
        }

        // load attributes to vars from file
        // ... pozdeji muzu udelat nejakou validaci ale zatim predpokladam ze jsou korektni 
        infile->getAttribute("dx").read(dx);
        infile->getAttribute("m_primary").read(m_primary);
        infile->getAttribute("r_in").read(r_in);
        infile->getAttribute("r_out").read(r_out);

        // allocate grid
        grid            = alloc_2D_double(dims);

        // read data to grid
        infile->getDataSet(dkey).read((double**) grid[0]);
    } else { // sim starts by setting parameters
        // from CLAs to easilly accesible vars
        m_primary       = p->d("m_primary") * M_SUN;
        dx              = p->d("dx");
        x               = p->d("x");
        r_in            = p->d("r_in");
        r_out           = p->d("r_out");

        // first and last step
        // CLAs ignored if sim starts without external data
        first_step      = 1;
        last_step       = first_step + s;
        
        // allocate grid
        grid            = alloc_2D_double(dims);
        
        // init 'empty' grid
        grid_init(idim, jdim, r_in, r_out, dx, grid);
    }

    // File attributes (sim metadata)
    if (append) {
        outfile = infile;
    } else {
        outfile = new H5::File(p->s("output"), H5::File::Overwrite);
        outfile->createAttribute<int>("idim", H5::DataSpace::From(idim)).write(idim);
        outfile->createAttribute<int>("jdim", H5::DataSpace::From(jdim)).write(jdim);
        outfile->createAttribute<double>("r_in", H5::DataSpace::From(r_in)).write(r_in);
        outfile->createAttribute<double>("r_out", H5::DataSpace::From(r_out)).write(r_out);
        outfile->createAttribute<double>("m_primary", H5::DataSpace::From(m_primary)).write(m_primary);
        outfile->createAttribute<double>("dx", H5::DataSpace::From(dx)).write(dx);
    }

    // allocate drain 
    std::vector<size_t> drain_dims  = {s, idim};
    double** drain                  = alloc_2D_double(drain_dims);
    for (uint stp=0; stp<drain_dims[0]; stp++)
        for (uint i=0; i<drain_dims[1]; i++)
            drain[stp][i] = 0.0;

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
    Distributor* dist = new Distributor(idim, jdim, cdim[0], cdim[1], drain, dx);
    dist->setRotationProfile(profile);

    /** Save initial state to main dataset */
    if (!append) {
        H5::DataSet* ds = new H5::DataSet(outfile->createDataSet<double>("data_0", H5::DataSpace(dims)));
        ds->write((double**) grid[0]);
        delete ds;
    }

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
    double** data_send_master = alloc_2D_double(cdim[0], cdim[1]);
    double** data_recv_master = alloc_2D_double(cdim[0], cdim[1]);

    /** Percentage info msg. variables */
    uint percent = 0;
    uint one_percent = s / 100;

    /** Compute all simulation steps */
    for (uint step=first_step; step <= last_step; step++) {
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

        /** Save recieved data to HDF file in step specifcdim[0]h dataset */
        std::ostringstream ss;
        ss << "data_" << step;
        std::string dname = ss.str();
        if (!outfile->exist(ss.str())) {
            outfile->createDataSet<double>(dname, H5::DataSpace(dims));
        }
        H5::DataSet* ds = new H5::DataSet(outfile->getDataSet(dname));
        ds->write((double**) grid[0]);
        delete ds;

        /** Run distribution handler object on grid */
        dist->step(grid);

        /** Print percentage msg. */
        if (s >= 100) {
            if (percent != (step - first_step) / one_percent) {
                percent = (step - first_step) / one_percent;
                std::cout << "Running ... " << percent << "%\t\r" << std::flush;
            }
        }
    }

    // save last step attr
    if (!outfile->hasAttribute("last_step"))
        outfile->createAttribute<uint>("last_step", H5::DataSpace::From(last_step)).write(last_step);
    outfile->getAttribute("last_step").write(last_step);

    // drain dataset
    // write to first non-existent dataset name
    int d = 0;
    while (outfile->exist("drain_" + std::to_string(d)))
        d += 1;
    outfile->createDataSet<double>("drain_" + std::to_string(d), H5::DataSpace(drain_dims));
    H5::DataSet* drain_dataset = new H5::DataSet(outfile->getDataSet("drain_" + std::to_string(d)));
    drain_dataset->write((double**) drain[0]);

    /** Stop all workers (slave) by sending STOP flag (and dummy data) */
    double** dummy = alloc_2D_double(cdim[0], cdim[1]);
    for (unsigned int i=1; i <= n_workers; i++) {
        MPI_Send(&dummy[0][0], cdim[0]*cdim[1], MPI_DOUBLE, i, STOP, MPI_COMM_WORLD); 
    }
}

ArgumentParser* createArgumentParser() {
    ArgumentParser* p = new ArgumentParser();

    // add boolean args.
    p->addArgument( new Argument<bool>("append", false));

    // add integer args.
    p->addArgument( new Argument<int>("s", 5e5));
    p->addArgument( new Argument<int>("idim"));
    p->addArgument( new Argument<int>("jdim"));
    p->addArgument( new Argument<int>("first_step", 0));
    p->addArgument( new Argument<int>("last_step"));

    // add double args.
    p->addArgument( new Argument<double>("dx", 0.01));
    p->addArgument( new Argument<double>("x", 0.0));
    p->addArgument( new Argument<double>("m_primary", 1.0));
    p->addArgument( new Argument<double>("r_in", 6.96e8));
    p->addArgument( new Argument<double>("r_out", 50.0 * 6.96e8));

    // add string args.
    p->addArgument( new Argument<std::string>("output"));
    p->addArgument( new Argument<std::string>("input"));

    return p;
}

int main(int argc, char **argv) {
    // Command line args. parser
    ArgumentParser* p = createArgumentParser();
    if (!p->parse(argc, argv))
        std::cout << p; // prints out help msg.
        return 0;

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
    while (cdim[0] * n_workers < n_jobs) // int. rounding correction
        cdim[0]++;
    
    /* Process specific task */
    if (p_rank == MASTER) {
        MPI_master(cdim, n_workers, p);
    } else {
        MPI_slave(cdim);
    }
    
    MPI_Finalize();

    delete p;

    return 0;
}
