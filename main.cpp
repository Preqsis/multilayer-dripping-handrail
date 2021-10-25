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

#include "InputParser.h"
#include "Argument.h"
#include "MSMM.h"
#include "Distributor.h"

typedef boost::array<double, 2> state_type;
typedef long unsigned int lint;
typedef unsigned int uint;

// Physical constants
const double M_SUN = 1.9891e30;
const double G = 6.6743e-11;

// Compute flags
const int MASTER = 0, COMPUTE=1, STOP=2;

// alloc 2d double pointer array
double** alloc_2D_double(int idim, int jdim) {
    double** arr;
    arr = (double**)calloc(idim, sizeof(double*));
    arr[0] = (double*)calloc(idim*jdim, sizeof(double));
    for (int i = 0; i < idim; i++) {
        arr[i] = arr[0] + i * jdim;
    }
    return arr;
}

// Worker function
void worker(double** data, MSMM2* model, uint ic) {
    // solver / stepper
    //boost::numeric::odeint::runge_kutta4<state_type> stepper;
    //boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;

    // pres vsechny zadane 'bunky'
    for (uint i=0; i<ic; i++) {
        // bunku nepocitat!
        if (data[i][10] == 0.0) {continue;}

        // bunka s nulovou hmotnosti
        // rovnice nedavaji smysl
        if (data[i][5] == 0.0) {continue;}

        // pocatectni podminky
        double x = data[i][2];
        state_type y = {data[i][3], data[i][4]};

        // velikost kroku
        double dx = data[i][8];

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
        data[i][5] += data[i][6]; // zapocteni pritoku k hmotnosti 
    }
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
            data[k][5] = 0.0; // m
            data[k][6] = 0.0; // dm
            data[k][7] = std::pow(r_out, 2.0) / std::pow(r_in + (idim - i - 1) * (r_out - r_in) / (idim - 1), 2.0); // r (ring specific radius)
            data[k][8] = dx; // dx (inner MSMM parameter)
            data[k][9] = 2.0 * M_PI * ((double) j) / ((double) jdim); // azimuth (cell specific rotation angle)
            data[k][10] = 1.0; // compute flag (0.0 --> do not compute)
        }
    }
}

// Commad line arguments defs.
Argument<lint>* _n = new Argument<lint>("n", "Number of simulation steps (default 5E5)", 5e5);
Argument<lint>* _idim = new Argument<lint>("idim", "Number of rings");
Argument<lint>* _jdim = new Argument<lint>("jdim", "Number of cells in one ring");
Argument<double>* _dx = new Argument<double>("dx", "Independent var step (default 0.01)", 0.01);
Argument<double>* _x = new Argument<double>("x", "Independent var initial value (default 0.0)", 0.0);
Argument<double>* _m_primary = new Argument<double>("m_primary", "Primary mass in Sun masses (default 1.0)", 1.0);
Argument<double>* _r_in = new Argument<double>("r_in", "Inner disc boundary radius (default 6.96e8 m)", 6.96e8);
Argument<double>* _r_out = new Argument<double>("r_out", "Outer disc boundary radius (default 50 * 6.96e8 m)", 50.0 * 6.96e8);
Argument<std::string>* _owner = new Argument<std::string>("owner", "Specified owner/creator of the outputed data file (default 'without owner')", "no owner");
Argument<std::string>* _output = new Argument<std::string>("output", "Output data file");
Argument<std::string>* _input = new Argument<std::string>("input", "Input data file");
Argument<std::string>* _dkey = new Argument<std::string>("dkey", "Selected input dataset key/name");

int main(int argc, char **argv) {
    // used vars
    uint    idim; // disc layer dimension
    uint    jdim; // disc cell dimension
    uint    ic; // mpi comms matrix i dimension
    uint    jc; // mpi comms matrix j dimension
    uint    n; // num of sim steps
    int     p_rank; // mpi process rank
    int     size; // num of mpi processes
    int     n_workers; // num of worker (aka. SLAVE) procesess
    double m_primary, dx, x, r_in, r_out; // system params
    std::vector<size_t> dims; //HDF5 dataset dims
    std::string output, input, dkey, owner;
    double** grid;

    // CLA parser
    InputParser parser(
        {_n, _idim, _jdim}, // integers 
        {_dx, _x, _m_primary, _r_in, _r_out},  // doubles
        {_output, _owner, _input, _dkey} // strings
    );
    if (!parser.run(argc, argv)) {
        std::cout << parser; // prints out help msg.
        return 0;
    }

    // common vars.
    owner   = _owner->getValue();
    n       = _n->getValue();

    // decide on input type
    if (_input->hasValue() && _dkey->hasValue()) { // sim start of external data
        input   = _input->getValue();
        dkey    = _dkey->getValue();

        // input file
        H5::File infile(input, H5::File::ReadOnly);

        // dataset and grid (HDF5) dims
        H5::DataSet dset = infile.getDataSet(dkey);
        dims = dset.getDimensions();

        // load attributes to vars from file
        // ... pozdeji muzu udelat nejakou validaci ale zatim predpokladam ze jsou korektni 
        infile.getAttribute("dx").read(dx);
        infile.getAttribute("idim").read(idim);
        infile.getAttribute("jdim").read(jdim);
        infile.getAttribute("m_primary").read(m_primary);
        infile.getAttribute("owner").read(owner);
        infile.getAttribute("r_in").read(r_in);
        infile.getAttribute("r_out").read(r_out);

        // allocate grid
        grid = alloc_2D_double(dims[0], dims[1]);

        // read data to grid
        dset.read((double**) grid[0]);
    } else { // sim starts by setting parameters
        // requird vars!!
        if (!_idim->hasValue() || !_jdim->hasValue() || !_output->hasValue() || !_n->hasValue()) {
            // pozdeji muzu udelat elegantnejsi reseni
            std::cout << "Missing some command line arguments!" << std::endl;
            return 0;
        }

        // from CLAs to easilly accesible vars
        idim            = _idim->getValue();
        jdim            = _jdim->getValue();
        m_primary       = _m_primary->getValue() * M_SUN;
        dx              = _dx->getValue();
        x               = _x->getValue();
        r_in            = _r_in->getValue();
        r_out           = _r_out->getValue();

        // grid dimms
        dims = {idim*jdim, 11};

        // allocate grid
        grid = alloc_2D_double(dims[0], dims[1]);

        // init 'empty' grid
        grid_init(idim, jdim, r_in, r_out, dx, grid);
    }

    // output check and value
    if (!_output->hasValue()) {
        std::cout << "Output HDF5 file not specified!" << std::endl;
        return 0;
    }
    output  = _output->getValue();

    // MPI init
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    n_workers = size-1; // number of MPI worker processes 
    ic      = idim * jdim / n_workers; // comm dim i
    jc      = dims[1]; // comm dim j

    /* Process specific task */
    if (p_rank == MASTER) { // main control process
        /* create HDF output file */
        H5::File outfile(output, H5::File::ReadWrite | H5::File::Create | H5::File::Overwrite);

        // File attributes (sim metadata)
        outfile.createAttribute<uint>("n", H5::DataSpace::From(n)).write(n);
        outfile.createAttribute<uint>("idim", H5::DataSpace::From(idim)).write(idim);
        outfile.createAttribute<uint>("jdim", H5::DataSpace::From(jdim)).write(jdim);
        outfile.createAttribute<double>("r_in", H5::DataSpace::From(r_in)).write(r_in);
        outfile.createAttribute<double>("r_out", H5::DataSpace::From(r_out)).write(r_out);
        outfile.createAttribute<double>("m_primary", H5::DataSpace::From(m_primary)).write(m_primary);
        outfile.createAttribute<double>("dx", H5::DataSpace::From(dx)).write(dx);
        outfile.createAttribute<std::string>("owner", H5::DataSpace::From(owner)).write(owner);

        /**
         * 'Drain' dataset
         * - records step by step mass 'drops'
         * - can be used to compute ring specific energy output
         **/
        double** drain = alloc_2D_double(n, idim);
        for (uint stp=0; stp < n; stp++) {
            for (uint i=0; i < idim; i++) {
                drain[stp][i] = 0.0;
            }
        }

        /**
         * 'Rotation profile' 
         * - specifies angular velocity for each ring relative to ring i == 0
         */
        std::vector<double> profile;
        double T_out = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r_out, 3.0)) / (G * m_primary) , 1.0/2.0);
        double r, T;
        for (uint i=0; i < idim; i++) {
            r = r_in + (idim - i - 1) * (r_out - r_in) / (idim - 1); // ring specific radius
            T = std::pow((4.0 * std::pow(M_PI, 2.0) * std::pow(r, 3.0)) / (G * m_primary) , 1.0/2.0); // ring spe
            profile.push_back(std::round((std::pow(10.0, 4.0) * T_out * 2.0 * M_PI) / (T * ((double)jdim))) / std::pow(10.0, 4.0));
        }

        /**
         * 'Distributor' object
         * - handles mass redistribution (free flow and 'dripping')
         */
        Distributor* dist = new Distributor(idim, jdim, ic, jc, drain, dx);
        dist->setRotationProfile(profile);

        /** Save initial state to main dataset */
        H5::DataSet* ds = new H5::DataSet(outfile.createDataSet<double>("data_0", H5::DataSpace(dims)));
        ds->write((double**) grid[0]);
        delete ds;

        /* Create jobs set for each worker */
        uint k;
        std::vector<std::vector<uint>> jobs = {};
        for (uint w = 0; w < n_workers; w++) {
            jobs.push_back({});
            for (uint i = 0; i < ic; i++) {
                k = w * ic + i;
                jobs[w].push_back(k);
            }
        }

        /** Communication matrix allocation */
        double** data_send_master = alloc_2D_double(ic, jc);
        double** data_recv_master = alloc_2D_double(ic, jc);
        
        /** Percentage info msg. variables */
        uint percent = 0;
        uint one_percent = n / 100;

        /** Compute all simulation steps */
        for (uint step=1; step <= n; step++) {
            /** Sort, mark (compute flag) and send jobs to specific workers */
            uint w;
            for (uint slave=1; slave<=n_workers; slave++) {
                w = slave-1;
                for (uint i=0; i<jobs[w].size(); i++) {
                    k = jobs[w][i];

                    if (k >= idim*jdim) { // outside of grid range --> compute = 0.0
                        data_send_master[i][jc-1] = 0.0;
                        continue;
                    }

                    /** actual grid / cell data to communication matrix */
                    for (uint l=0; l<jc; l++) {
                        data_send_master[i][l] = grid[k][l];
                    }
                    data_send_master[i][jc-1] = 1.0; // set compute flag --> compute = 1.0 (to be computed)
                }
                MPI_Send(&data_send_master[0][0], ic*jc, MPI_DOUBLE, slave, COMPUTE, MPI_COMM_WORLD); // send data to slave
            }

            /** Recieve and sort data back to grid */
            for (uint slave=1; slave <= n_workers; slave++) {
                MPI_Recv(&data_recv_master[0][0], ic*jc, MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // recv. data

                for (uint i=0; i<ic; i++) { // back to grid
                    if (data_recv_master[i][10] == 0.0) {continue;} // compute == 1.0 only
                    w = slave-1;
                    k = w * ic + i;
                    for (uint l=0; l<jc; l++) {
                        grid[k][l] = data_recv_master[i][l];
                    }
                }
            }

            /** Save recieved data to HDF file in step specifich dataset */
            std::ostringstream ss;
            ss << "data_" << step;
            H5::DataSet* ds = new H5::DataSet(outfile.createDataSet<double>(ss.str(), H5::DataSpace(dims)));
            ds->write((double**) grid[0]);
            delete ds;

            /** Run distribution handler object on grid */
            dist->step(grid);

            /** Print percentage msg. */
            if (n >= 100) {
                if (percent != step / one_percent) {
                    percent = step / one_percent;
                    std::cout << "Running ... " << percent << "%\t\r" << std::flush;
                }
            }
        }

        /** Save 'drain' dataset to HDF file */
        std::vector<size_t> drain_dims = {n, idim};
        H5::DataSet* DrainDataSet = new H5::DataSet(outfile.createDataSet<double>("drain", H5::DataSpace(drain_dims)));
        DrainDataSet->write((double**) drain[0]);
        delete DrainDataSet;

        /** Stop all workers (slave) by sending STOP flag (and dummy data) */
        double** dummy = alloc_2D_double(ic, jc);
        for (unsigned int i=1; i <= n_workers; i++) {
            MPI_Send(&dummy[0][0], ic*jc, MPI_DOUBLE, i, STOP, MPI_COMM_WORLD); 
        }
    } else {
        /** Communication matrix allocation */
        double** data_recv_slave = alloc_2D_double(ic, jc);
        double** data_send_slave = alloc_2D_double(ic, jc);
        
        // mass-spring model
        MSMM2* model = new MSMM2();

        /** Recieve until STOP */
        while (true) {
            /** Recieve data */
            MPI_Recv(&data_recv_slave[0][0], ic*jc, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            /** If STOP --> end slave computation */
            if (status.MPI_TAG == STOP) {
                free(model);
                break;
            }

            /** Run worker on recieved data */
            worker(data_recv_slave, model, ic);

            /** Send results to master */
            MPI_Send(&data_recv_slave[0][0], ic*jc, MPI_DOUBLE, 0, COMPUTE, MPI_COMM_WORLD);
        }
    }
    
    MPI_Finalize();

    return 0;
}
