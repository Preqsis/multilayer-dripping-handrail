#ifndef RADIATION_CPP
#define RADIATION_CPP

#include <iostream>
#include <mpi.h>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include "functions.cpp"

namespace Radiation {

typedef std::vector<std::string> stype;
typedef std::vector<stype> jtype;

void terminate(std::vector<size_t> dim, int n_workers) {
    double*** data = fn::alloc_3D_double(dim);
    for (int slave = 1; slave <= n_workers; slave++) {
        MPI_Send(&data[0][0][0], dim[0]*dim[1]*dim[2], MPI_DOUBLE, slave, STOP, MPI_COMM_WORLD); 
    }
}

bool is_valid(ArgumentParser* p) {
    return p->isSet("step_first") && p->isSet("step_last") && p->isSet("mass_file");
}

void spectrum(double*** data_mass, double**** data_spec, std::vector<size_t> dim_mass, std::vector<size_t> dim_spec) {
    for (size_t i = 0; i < dim_spec[0]; i++) {
        for (size_t j = 0; j < dim_spec[1]; j++) {
            for (size_t k = 0; k < dim_spec[2]; k++) {
                

                data_spec[i][j][k][1] = 1.0;


            }
        }
    }
}

void slave(std::vector<size_t> dim_mass, std::vector<size_t> dim_spec) {
    /** MPI status flag **/
    MPI_Status status;

    // Comm length
    size_t len_mass = dim_mass[0] * dim_mass[1] * dim_mass[2];
    size_t len_spec = dim_spec[0] * dim_spec[1] * dim_spec[2] * dim_spec[3];

    // Comms data allocation
    double*** data_mass     = fn::alloc_3D_double(dim_mass);
    double**** data_spec    = fn::alloc_4D_double(dim_spec);

    while (true) {
        MPI_Recv(&data_mass[0][0][0], len_mass, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == STOP) {break;}

        if (status.MPI_TAG == COMPUTE) {
            spectrum(data_mass, data_spec, dim_mass, dim_spec);
        }
        
        MPI_Send(&(data_spec[0][0][0][0]), len_spec, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
    }
}

void master(std::vector<size_t> dim_mass, std::vector<size_t> dim_spec, int n_workers, ArgumentParser* p) {
    // MPI status flag holder
    MPI_Status status;

    if (is_valid(p)) {
        std::string dkey;
        bool skip;
        int slave;
        int shift;
        int s           = p->i("step_first");
        int step_last   = p->i("step_last");

        // input drain file / output spectrum file
        H5::File* mass_file = new H5::File(p->s("mass_file"), H5::File::ReadOnly);
        H5::File* spec_file = new H5::File(p->s("outdir") + "/spectrum.h5", H5::File::Overwrite);

        // Comm length
        size_t len_mass = dim_mass[0] * dim_mass[1] * dim_mass[2];
        size_t len_spec = dim_spec[0] * dim_spec[1] * dim_spec[2] * dim_spec[3];

        // Comms data allocation
        double*** data_mass     = fn::alloc_3D_double(dim_mass);
        double**** data_spec    = fn::alloc_4D_double(dim_spec);

        // Spectrum range init
        double lam_low  = p->d("lam_low");
        double lam_step = p->d("lam_step");
        for (size_t i = 0; i < dim_spec[0]; i++){
            for (size_t j = 0; j < dim_spec[1]; j++){
                for (size_t k = 0; k < dim_spec[2]; k++){
                    data_spec[i][j][k][0] = lam_low + ((double) k) * lam_step;
                }
            }
        }

        while (true) {
            // out
            for (slave = 1; slave <= n_workers; slave++) {
                shift   = slave - 1;
                dkey    = "data_" + std::to_string(s + shift);
                skip    = !mass_file->exist(dkey) || s + shift > step_last;

                if (!skip) {
                    mass_file->getDataSet(dkey).read((double***) data_mass[0][0]);
                }

                MPI_Send(&data_mass[0][0][0], len_mass, MPI_DOUBLE, slave, (skip) ? SKIP : COMPUTE, MPI_COMM_WORLD); 
            }

            // in
            for (slave = 1; slave <= n_workers; slave++) {
                shift   = slave - 1;
                dkey    = "data_" + std::to_string(s + shift);

                MPI_Recv(&data_spec[0][0][0][0], len_spec, MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                if (status.MPI_TAG == COMPUTE) {
                    fn::writeDataSet(spec_file, data_spec, dim_spec, dkey); // initial state save
                }
            }

            s += n_workers;

            if (s > step_last) {break;}
        }
    }

    terminate(dim_mass, n_workers);
}

}

#endif

