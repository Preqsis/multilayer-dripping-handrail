#include <iostream>
#include <mpi.h>
#include <math.h>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include "argparse-cpp/ArgumentParser.h"

#include "Functions.hpp"
namespace fn = Functions;

#include "Constants.hpp"
namespace cs = Constants;

#include "Observation.hpp"

void Observation::terminate(std::vector<size_t> dim, int n_workers) {
    double**** data = fn::alloc_4D_double(dim);
    for (int slave = 1; slave <= n_workers; slave++) {
        MPI_Send(&data[0][0][0][0], dim[0]*dim[1]*dim[2]*dim[3], MPI_DOUBLE, slave, cs::mpi::STOP, MPI_COMM_WORLD); 
    }
}

double Observation::filter_gauss(double**** data, std::vector<size_t> dim, double mu, double fwhm) {
    double val  = 0.0;
    for (int i=0; i<dim[0]; i++) {
        for (int j=0; j<dim[1]; j++) {
            for (int k=0; k<dim[2]; k++) {
                //val += data[i][j][k][1] * std::exp(-0.5 * pow((data[i][j][k][0] - mu) * 2.355 / fwhm, 2.0));
                val += data[i][j][k][1] * std::exp(-0.5 * pow((data[i][j][k][0] - mu) / fwhm, 2.0));
            }
        }
    }
    return val;
}

// Function for "obs" MPI slave processes
void Observation::slave(std::vector<size_t> dim_rad, std::vector<size_t> dim_obs, ArgumentParser* p) {
    // MPI status flag holder
    MPI_Status status;

    // Comm length
    size_t len_rad = dim_rad[0] * dim_rad[1] * dim_rad[2] * dim_rad[3];
    size_t len_obs = dim_obs[0];

    // Comms data allocation
    double**** data_rad = fn::alloc_4D_double(dim_rad);
    double* data_obs     = fn::alloc_1D_double(dim_obs);

    while (true) {
        MPI_Recv(&data_rad[0][0][0][0], len_rad, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == cs::mpi::STOP) {break;}

        data_obs[0] = filter_gauss(data_rad, dim_rad, 365.6e-7, 34.0e-7); // U
        data_obs[1] = filter_gauss(data_rad, dim_rad, 435.3e-7, 78.1e-7); // B
        data_obs[2] = filter_gauss(data_rad, dim_rad, 547.7e-7, 99.1e-7); // V
        data_obs[3] = filter_gauss(data_rad, dim_rad, 634.9e-7, 106.56e-7); // R
        //data_obs[3] = filter_gauss(data_rad, dim_rad, 634.9e-9, 106.56e-9); // I

        MPI_Send(&(data_obs[0]), len_obs, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
    }
}

// Function for "obs" MPI master process
void Observation::master(std::vector<size_t> dim_rad, std::vector<size_t> dim_obs, int n_workers, ArgumentParser* p) {
    // MPI status flag holder
    MPI_Status status;

    std::string dkey;
    bool skip;
    int slave;
    int shift;
    int step_first  = p->i("step_first");
    int s           = step_first;
    int step_last   = p->i("step_last") > 0 ? p->i("step_last") : p->i("n")-1;
    int one_percent = (step_last - step_first) / 100;
    
    std::vector<size_t> dim_out = {step_last - step_first + 1, dim_obs[0]+1};
    int i = 0;

    // verbosity?
    bool v = p->b("verbose");

    // input rad file / output obs file
    std::string rad_fpath = (p->isSet("rad_file")) ? p->s("rad_file") : p->s("outdir") + "/rad.h5";
    H5::File* rad_file = new H5::File(rad_fpath, H5::File::ReadOnly);
    H5::File* obs_file = new H5::File(p->s("outdir") + "/obs.h5", H5::File::Overwrite);

    // Comm length
    size_t len_rad = dim_rad[0] * dim_rad[1] * dim_rad[2] * dim_rad[3];
    size_t len_obs = dim_obs[0];

    // Comms data allocation
    double**** data_rad    = fn::alloc_4D_double(dim_rad);
    double* data_obs        = fn::alloc_1D_double(dim_obs);
    double** data           = fn::alloc_2D_double(dim_out);

    while (true) {
        // out
        for (slave = 1; slave <= n_workers; slave++) {
            shift   = slave - 1;
            dkey    = "d" + std::to_string(s + shift);
            skip    = !rad_file->exist(dkey) || s + shift > step_last;
            
            if (!skip) {
                rad_file->getDataSet(dkey).read((double****) data_rad[0][0][0]);
            }

            MPI_Send(&data_rad[0][0][0][0], len_rad, MPI_DOUBLE, slave, (skip) ? cs::mpi::SKIP : cs::mpi::COMPUTE, MPI_COMM_WORLD); 
        }

        // in
        for (slave = 1; slave <= n_workers; slave++) {
            shift   = slave - 1;
            dkey    = "d" + std::to_string(s + shift);

            MPI_Recv(&data_obs[0], len_obs, MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == cs::mpi::COMPUTE) {
                data[i+shift][0] = s+shift;
                data[i+shift][1] = data_obs[0]; // U
                data[i+shift][2] = data_obs[1]; // B
                data[i+shift][3] = data_obs[2]; // V
                data[i+shift][4] = data_obs[3]; // R
                data[i+shift][5] = data_obs[4]; // I
            }
        }

        // info msg.
        if (v) {
            if (one_percent == 0 || (s - step_first) % one_percent == 0) {
                std::cout << "Observation ... " << (int) (100.0 * (double) (s - step_first) / (step_last - step_first)) << "%\t\r" << std::flush;
            }
        }
        
        s += n_workers;
        i += n_workers;
        if (s > step_last) {break;}
    }

    fn::writeDataSet(obs_file, data, dim_out, "data"); // initial state save

    terminate(dim_rad, n_workers);

    delete rad_file;
    delete obs_file;

    // info msg.
    if (v) {
        std::cout << std::endl << "Done ..." << std::endl;
    }
}

