#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cstdlib>

#include <H5DataSet.hpp>
#include <H5DataSpace.hpp>
#include <H5File.hpp>
namespace H5 = HighFive;

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include "InputParser.h"
#include "Argument.h"
#include "MSMM.h"

typedef boost::array<double, 2> state_type;
typedef long unsigned int lint;
typedef unsigned int uint;


/**
 * 'Communication' double** matrix allocation
 *
 * @param idim Number of rows in communication matrix.
 * @param jdim Number of columns in communication matrix.
 * @return allocated double** array 'matrix'.
 */
double** alloc_2D_double(int idim, int jdim) {
    double** arr;
    arr = (double**)calloc(idim, sizeof(double*));
    arr[0] = (double*)calloc(idim*jdim, sizeof(double));
    for (int i = 0; i < idim; i++) {
        arr[i] = arr[0] + i * jdim;
    }
    return arr;
}

/**
 * Worker function
 *
 * @param container 
 *
 */
double** worker(double** container, uint ic) {
    // solver / stepper
    //boost::numeric::odeint::runge_kutta4<state_type> stepper;
    //boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;

    // model
    MSMM2* model = new MSMM2();

    // pres vsechny zadane 'bunky'
    for (uint i=0; i<ic; i++) {
        // bunku nepocitat!
        if (container[i][10] == 0.0) {continue;}

        // pocatectni podminky
        double x = container[i][2];
        state_type y = {container[i][3], container[i][4]};

        // velikost kroku
        double dx = container[i][8];

        // parametry do modelu
        model->set_m(container[i][5]);
        model->set_dm(container[i][6] / dx); // pritok urceny vstupem, podeleno krokem aby davalo smysl pro rce.
        model->set_g(container[i][7]);

        // krok
        stepper.do_step(*model, y, x, dx);

        // vysledky zpet do containeru
        container[i][2] += dx; // inkrementace casu
        container[i][3] = y[0]; // z
        container[i][4] = y[1]; // v
        container[i][5] += container[i][6]; // zapocteni pritoku k hmotnosti 
    }
    
    free(model);

    return container;
}

class Distributor {
private:
    uint _idim;
    uint _jdim;
    uint _ic;
    uint _jc;
    uint _n;
    double _dx;
    double _q;
    double _zc;
    double** _drain;

    std::vector<uint> _probs;
    
    std::vector<double> RingTemp;

    std::vector<double> _RotationProfile;
public:
    Distributor(uint idim, uint jdim, uint ic, uint jc, double** drain, double dx) {
        _idim = idim;
        _jdim = jdim;
        _ic = ic;
        _jc = jc;
        _dx = dx;
        _q = 0.5;
        _zc = 5.5;
        _drain = drain;
        _n = 0;

        // sestaveni pravdepodobnostniho pole
        uint n;
        uint i = 0;
        for (uint j=0; j < _jdim; j++) {
            n = (uint) ((1.0 - (double) j / _jdim) * 100.0);
            for (uint k=0; k < n; k++) {
                _probs.push_back(j);
                i++;
            }
        }

        // 'random' seed
        srand(time(NULL));
    }

    void setRotationProfile(std::vector<double> profile) {
        _RotationProfile.clear();
        _RotationProfile = profile;
    }

    void step(double** grid) {
        // rotace vnejsiho prstence o "jeden pohyb"
        // zabranuje driftovani vliven nepresnosti datovych typu
        uint k;
        for (uint j = 0; j < _jdim; j++) {
            k = j;
            RingTemp.push_back(grid[k][9]);
        }
        std::rotate(RingTemp.begin(), RingTemp.begin() + 1, RingTemp.end());
        for (uint j = 0; j < _jdim; j++) {
            k = j;
            grid[k][9] = RingTemp[j];
        }
        RingTemp.clear();

        // rotace ostatnich podle profilu
        for (uint i = 1; i < _idim; i++) {
            for (uint j = 0; j < _jdim; j++) {
                k = i * _jdim + j;
                grid[k][9] = std::fmod(grid[k][9] + _RotationProfile[i], 2.0 * M_PI);
            }
        }

        // nalezeni pritokove bunky + pritok
        uint idx = 0;
        for (uint j = 0; j < _jdim; j++) {
            k = j;
            if (grid[k][9] == 0.0) {
                idx = k;
                break;
            }
        }
        grid[idx][5] += _q;

        double mb;
        double dp;
        uint k_out, k_in;
        uint k_lower_right, k_lower_left;
        double j_in, j_out;
        double part_lower_left, part_lower_right;
        for (uint i=_idim-1; i < _idim; i--) { // od stredu ven, unsigned preskakuje na maximalni hodnotu!!!
            if (i < _idim - 1) {
                k_out = i * _jdim + 0;
                k_in = (i + 1) * _jdim + 0;
                dp = (double)_jdim * (grid[k_out][9] - grid[k_in][9]) / (2.0 * M_PI);

                for (uint j = 0; j < _jdim; j++) {
                    k = i * _jdim + j;

                    if (grid[k][3] < _zc) { // k odtreni nedochazi
                        continue;
                    }

                    // urceni vnitrnich prilehajicih bunek
                    j_in = (double)j + dp;
                    j_in = (j_in > 0.0) ? j_in : (double)_jdim + j_in;

                    k_lower_left = (i + 1) * _jdim + (uint)std::floor(j_in) % _jdim;
                    k_lower_right = (i + 1) * _jdim + (uint)std::floor(j_in + 1.0) % _jdim;
                    part_lower_left = 1.0 - std::fmod(j_in, 1.0);
                    part_lower_right = std::fmod(j_in, 1.0);
                    
                    //std::cout << j << "\t" << k_lower_left << "\t" << k_lower_right << "\t" << part_lower_left << "\t" << part_lower_right << std::endl;
                    
                    // odtrzena hmotnost
                    mb = grid[k][5] - (0.2 * grid[k][5] + 0.3);

                    // odebrat od zdrojove bunky
                    grid[k][5] -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    grid[k][4] = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    grid[k][3] = 2.0;

                    grid[k_lower_left][5] += part_lower_left * mb;
                    grid[k_lower_right][5] += part_lower_right* mb;

                    // zapsat do drain datasetu
                    _drain[_n][i] += mb;
                }
            } else {
                for (uint j = 0; j < _jdim; j++) {
                    k = i * _jdim + j;

                    if (grid[k][3] < _zc) { // k odtreni nedochazi
                        continue;
                    }

                    // odtrzena hmotnost
                    mb = grid[k][5] - (0.2 * grid[k][5] + 0.3);

                    // odebrat od zdrojove bunky
                    grid[k][5] -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    grid[k][4] = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    grid[k][3] = 2.0;
                    
                    // hmotu na drain
                    _drain[_n][i] += mb;
                }
            }
        }

        _n += 1;
        
    }
};

/**
 * Used physical constants
 */
const double AU = 1.495978707e11;
const double M_SUN = 1.9891e30;
const double G = 6.6743e-11;

/**
 * Command line arguments defs
 **/
// long integers
Argument<lint>* _n = new Argument<lint>("n", 5e5, "Number of simulation steps (default 5E5)");
Argument<lint>* _idim = new Argument<lint>("idim", 10, "Number of rings (default 10)");
Argument<lint>* _jdim = new Argument<lint>("jdim", 10, "Number of cells in one ring (default 10)");
// doubles
Argument<double>* _dx = new Argument<double>("dx", 0.01, "Independent var step (default 0.01)");
Argument<double>* _x = new Argument<double>("x", 0.0, "Independent var initial value (default 0.0)");
Argument<double>* _m_primary = new Argument<double>("m_primary", 27.0, "Primary mass in Sun masses (default 27.0)");
Argument<double>* _m_secondary = new Argument<double>("m_secondary", 16.0, "Secondary mas in Sun masses (default 16.0)");
Argument<double>* _d = new Argument<double>("d", 0.2, "Primary to secondary distance in AU (default 0.2)"); 
Argument<double>* _r_in = new Argument<double>("r_in", 6.96e8, "Inner disc boundary radius (default 6.96e8 m)");
// strings
Argument<std::string>* _output = new Argument<std::string>("output", "./data.h", "Output data file (default ./data.h)");

int main(int argc, char **argv) {
    /* Command line prarser definition */
    InputParser parser(
        {_n, _idim, _jdim}, // integers 
        {_dx, _x, _m_primary, _m_secondary, _d, _r_in},  // doubles
        {_output} // strings
    );
    if (!parser.run(argc, argv)) {
        std::cout << parser; // prints out help msg.
        return 0;
    }

    /* MPI initialization */
    int p_rank;
    int size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Number of paralel workers */
    int n_workers = size-1; 

    /* Flags defs. */
    const int MASTER = 0; // for simplicity
    const int COMPUTE=1, STOP=2;

    /* 'Real' simulation grid dimensions */ 
    uint idim = _idim->getValue();
    uint jdim = _jdim->getValue();

    /**
     * 'Comunication' grid dimension
     * - serialized and divided by workers
     */
    uint ic = idim*jdim  / n_workers + 1;
    lint jc = 11; // number of parameters for each cell

    /**
     * HDF dataset dimensions 
     * - whole grid
     * - constant size throughout the simulation! 
     */
    std::vector<size_t> dims = {idim*jdim, jc};
    
    /* Process specific task */
    if (p_rank == MASTER) { // main control process
        /* create HDF output file */
        H5::File file(_output->getValue(), H5::File::ReadWrite | H5::File::Create | H5::File::Overwrite);
        
        /**
         *  Independent parameters
         *  - cla input
            */
        double m_primary = _m_primary->getValue() * M_SUN;
        double m_secondary = _m_secondary->getValue() * M_SUN;
        double r_in = _r_in->getValue();
        double d = _d->getValue() * AU;

        /**
         * Dependent paramaters
         * - computed
         */
        double r_out = d * std::cbrt(m_secondary / (3.0 * (m_primary + m_secondary))); // L1 point

        /**
         * Simulation start info msg.
         */
        std::cout << "--- Multi-Layer Dripping Handrail Parameters ---" << std::endl << std::endl;
        std::cout << "m_primary \t\t= " << _m_primary->getValue() << " * M_SUN (" << m_primary << " kg)" << std::endl;
        std::cout << "m_secondary \t\t= " << _m_secondary->getValue() << " * M_SUN (" << m_secondary << " kg)" << std::endl;
        std::cout << "d \t\t\t= " << _d->getValue() << " AU (" << d << " m)" << std::endl;
        std::cout << "r_in \t\t\t= " << r_in << " m" << std::endl;
        std::cout << "r_out \t\t\t= " << r_out << " m" << std::endl;
        std::cout << "n \t\t\t= " << _n->getValue() << " steps" << std::endl;
        std::cout << "output \t\t\t= " << _output->getValue() << std::endl << std::endl;
        std::cout << "------------------------------------------------" << std::endl;

        /**
         * 'Empty' (initial state) communication grid
         */
        double** grid = alloc_2D_double(idim*jdim, jc);
        uint k; 
        for (uint i=0; i<idim; i++) { // rings
            for (uint j=0; j<jdim; j++){ // cells
                k = i * jdim + j; // serialized communication grid coordinate
                grid[k][0] = i; // ring coordinate
                grid[k][1] = j; // cell coordinate
                grid[k][2] = 0.0; // t
                grid[k][3] = 2.0; // z
                grid[k][4] = 0.0; // v
                grid[k][5] = 0.2; // m
                grid[k][6] = 0.0; // dm
                grid[k][7] = std::pow(r_out, 2.0) / std::pow(r_in + (idim - i - 1) * (r_out - r_in) / (idim - 1), 2.0); // r (ring specific radius)
                grid[k][8] = _dx->getValue(); // dx (inner MSMM parameter)
                grid[k][9] = 2.0 * M_PI * ((double) j) / ((double) jdim); // azimuth (cell specifich rotation angle)
                grid[k][10] = 1.0; // compute flag (0.0 --> do not compute)
            }
        }

        /**
         * 'Drain' dataset
         * - records step by step mass 'drops'
         * - can be used to compute ring specific energy output
         **/
        double** drain = alloc_2D_double(_n->getValue(), idim);
        for (uint stp=0; stp < _n->getValue(); stp++) {
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
        Distributor* dist = new Distributor(idim, jdim, ic, jc, drain, _dx->getValue());
        dist->setRotationProfile(profile);

        /** Save initial state to main dataset */
        H5::DataSet* ds = new H5::DataSet(file.createDataSet<double>("data_0", H5::DataSpace(dims)));
        ds->write((double**) grid[0]);
        delete ds;

        /* Create jobs set for each worker */
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
        uint n = _n->getValue();
        uint one_percent = n / 100;

        /** Compute all simulation steps */
        for (uint step=1; step <= _n->getValue(); step++) {
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
            H5::DataSet* ds = new H5::DataSet(file.createDataSet<double>(ss.str(), H5::DataSpace(dims)));
            ds->write((double**) grid[0]);
            delete ds;

            /** Run distribution handler object on grid */
            dist->step(grid);

            /** Print percentage msg. */
            if (_n->getValue() >= 100) {
                if (percent != step / one_percent) {
                    percent = step / one_percent;
                    std::cout << "Running ... " << percent << "%\t\r" << std::flush;
                }
            }
        }

        /** Save 'drain' dataset to HDF file */
        std::vector<size_t> drain_dims = {_n->getValue(), idim};
        H5::DataSet* DrainDataSet = new H5::DataSet(file.createDataSet<double>("drain", H5::DataSpace(drain_dims)));
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

        /** Recieve until STOP */
        while (true) {
            /** Recieve data */
            MPI_Recv(&data_recv_slave[0][0], ic*jc, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            /** If STOP --> end slave computation */
            if (status.MPI_TAG == STOP) {
                break;
            }

            /** Run worker on recieved data */
            data_send_slave = worker(data_recv_slave, ic);

            /** Send results to master */
            MPI_Send(&data_recv_slave[0][0], ic*jc, MPI_DOUBLE, 0, COMPUTE, MPI_COMM_WORLD);
        }
    }
    
    MPI_Finalize();

    return 0;
}
