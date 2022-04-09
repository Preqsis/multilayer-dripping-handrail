#include <iostream>
#include <mpi.h>
#include <math.h>

#include "argparse-cpp/ArgumentParser.h"

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include "Functions.hpp"
namespace fn = Functions;

#include "Constants.hpp"
namespace cs = Constants;

#include "Radiation.hpp"

typedef std::vector<std::string> stype;
typedef std::vector<stype> jtype;

void Radiation::terminate(std::vector<size_t> dim, int n_workers) {
    double*** data = fn::alloc_3D_double(dim);
    for (int slave = 1; slave <= n_workers; slave++) {
        MPI_Send(&data[0][0][0], dim[0]*dim[1]*dim[2], MPI_DOUBLE, slave, cs::mpi::STOP, MPI_COMM_WORLD); 
    }
}

double Radiation::planck(double wl, double T) {
    double a    = 2.0 * cs::h * std::pow(cs::c, 2.0);
    double b    = cs::h * cs::c / (wl * cs::k * T);
    double val  =  a / (std::pow(wl, 5.0) * (std::exp(b) - 1.0 ));
    return val;
}

void Radiation::slave(std::vector<size_t> dim_sim, std::vector<size_t> dim_rad, ArgumentParser* p) {
    /** MPI status flag **/
    MPI_Status status;

    // Comm length
    size_t len_mass = dim_sim[0] * dim_sim[1] * dim_sim[2];
    size_t len_spec = dim_rad[0] * dim_rad[1] * dim_rad[2] * dim_rad[3];

    // Comms data allocation
    double*** data_sim     = fn::alloc_3D_double(dim_sim);
    double**** data_rad    = fn::alloc_4D_double(dim_rad);

    // Spectrum range init
    double wl_low  = p->d("wl_low");
    double wl_step = p->d("wl_step");
    for (size_t i = 0; i < dim_rad[0]; i++){
        for (size_t j = 0; j < dim_rad[1]; j++){
            for (size_t k = 0; k < dim_rad[2]; k++){
                data_rad[i][j][k][0] = wl_low + ((double) k) * wl_step;
            }
        }
    }

    double r_in = p->d("r_in");
    double r_out = p->d("r_out");

    double dR = (r_out - r_in) / dim_sim[0];
    
    double wl, T, frq, S, A, E, T_a;
    double q = p->d("q");
    double Q = p->d("Q");

    double T_in = p->d("temp_in"); // central object temperature
    double T_atm = p->d("temp_atm"); // acreetion disk atmosphere temperature

    double r, B, T_r, E_r, L_a;

    while (true) {
        MPI_Recv(&data_sim[0][0][0], len_mass, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == cs::mpi::STOP) {break;}

        if (status.MPI_TAG == cs::mpi::COMPUTE) {
            for (size_t i = 0; i < dim_rad[0]; i++) {          // pres vsechny prstence
                r = data_sim[i][0][7]; // ring radius


                //S   = M_PI * (pow(r[i], 2.0) - pow(r[i+1], 2.0)) / dim_sim[1];
                // Plocha prstence
                // ... jedna strana
                S = 2.0 * M_PI * r * dR / dim_sim[1];

                T_r = T_in * pow(r / r_in ,-0.75); // radius dependent disk body temperature

                //std::cout << i << " -- " << T_r << std::endl;

                //A   = 0.5 * cs::G * cs::m_sun * p->d("m_primary") * dR * Q / q;
                A   = 0.5 * cs::G * cs::m_sun * p->d("m_primary") * dR;

                for (size_t j = 0; j < dim_rad[1]; j++) {      // pres vsechny bunky v prstenci
                    for (size_t k = 0; k < dim_rad[2]; k++) {  // pres rozsah vl. delek

                        // vlnova delka -> frekvence
                        wl  = data_rad[i][j][k][0];
                        frq = cs::c / wl;
                        
                        /**
                         * Direct Radiation 
                         * - free-free emission
                         * - 
                         */

                        E_r = A * data_sim[i][j][9] / std::pow(r, 2.0); // energine ztracena drainem

                        L_a = 4.9e-11 * (E_r / T_atm) * std::exp(-1.0 * cs::h * frq/ (cs::k * T_atm)); // zareni z drainu -> flickering

                        B = planck(wl, T_r);
                        
                        //data_rad[i][j][k][1] = /*2 * S * B +*/ L_a;
                        data_rad[i][j][k][1] = B;

                        
                        // Re-radiation
                        // on data index ...2


                        // "tepelne" razeni
                        //data_rad[i][j][k][1]   = 2.0 * S * planck(wl, T); // vln. delka, teplota

                        // flickering vlivem ztraty energie
                        //E                       = A * data_sim[i][j][9] / (r[i] * r[i+1]); // energine ztracena drainem
                        //if (data_sim[i][j][9] > 0.0) {
                        //    std::cout << A << ", " << E << std::endl; 
                        //}
                        //data_rad[i][j][k][2]   = 4.9e-11 * (E / T) * std::exp(-1.0 * cs::h * frq/ (cs::k * T)); // zareni z drainu -> flickering
                    }
                }
            }
        }
        
        MPI_Send(&(data_rad[0][0][0][0]), len_spec, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
    }
}

void Radiation::master(std::vector<size_t> dim_sim, std::vector<size_t> dim_rad, int n_workers, ArgumentParser* p) {
    // MPI status flag holder
    MPI_Status status;

    std::string dkey;
    bool skip;
    int slave;
    int shift;
    int step_first  = p->i("step_first");
    int s           = step_first;
    int step_last   = p->i("step_last");
    int one_percent = (step_last - step_first) / 100;

    // verbosity?
    bool v = p->b("verbose");

    // input drain file / output spectrum File
    std::string sim_fpath = (p->isSet("sim_file")) ? p->s("sim_file") : p->s("outdir") + "/sim.h5";
    H5::File* sim_file = new H5::File(sim_fpath, H5::File::ReadOnly);
    H5::File* rad_file = new H5::File(p->s("outdir") + "/rad.h5", H5::File::Overwrite);

    // Comm length
    size_t len_mass = dim_sim[0] * dim_sim[1] * dim_sim[2];
    size_t len_spec = dim_rad[0] * dim_rad[1] * dim_rad[2] * dim_rad[3];

    // Comms data allocation
    double*** data_sim     = fn::alloc_3D_double(dim_sim);
    double**** data_rad    = fn::alloc_4D_double(dim_rad);

    while (true) {
        // out
        for (slave = 1; slave <= n_workers; slave++) {
            shift   = slave - 1;
            dkey    = "data_" + std::to_string(s + shift);
            skip    = !sim_file->exist(dkey) || s + shift > step_last;

            if (!skip) {
                sim_file->getDataSet(dkey).read((double***) data_sim[0][0]);
            }

            MPI_Send(&data_sim[0][0][0], len_mass, MPI_DOUBLE, slave, (skip) ? cs::mpi::SKIP : cs::mpi::COMPUTE, MPI_COMM_WORLD); 
        }

        // in
        for (slave = 1; slave <= n_workers; slave++) {
            shift   = slave - 1;
            dkey    = "data_" + std::to_string(s + shift);

            MPI_Recv(&data_rad[0][0][0][0], len_spec, MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == cs::mpi::COMPUTE) {
                fn::writeDataSet(rad_file, data_rad, dim_rad, dkey); // initial state save
            }
        }
        
        // info msg.
        if (v) {
            if (one_percent == 0 || (s - step_first) % one_percent == 0) {
                std::cout << "Radiation ... " << (int) (100.0 * (double) (s - step_first) / (step_last - step_first)) << "%\t\r" << std::flush;
            }
        }
        
        s += n_workers;
        if (s > step_last) {break;}
    }
    
    // end slaves
    terminate(dim_sim, n_workers);

    delete sim_file;
    delete rad_file;

    // info msg.
    if (v) {
        std::cout << std::endl << "Done ..." << std::endl;
    }
}
