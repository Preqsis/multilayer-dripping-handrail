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

void terminate(std::vector<size_t> comm_dim, int n_workers) {
    double** data = fn::alloc_2D_double(comm_dim);
    for (int slave = 1; slave <= n_workers; slave++) {
        MPI_Send(&data[0][0], comm_dim[0]*comm_dim[1], MPI_DOUBLE, slave, STOP, MPI_COMM_WORLD); 
    }
}

bool is_valid(ArgumentParser* p) {
    return p->isSet("step_first") && p->isSet("step_last") && p->isSet("drain_file");
}

void slave(std::vector<size_t> comm_dim) {
    /** MPI status flag **/
    MPI_Status status;

    std::vector<size_t> dim_in  = {25, 157};
    std::vector<size_t> dim_out = {25, 157, 1000};

    /** Comms matrix allocation */
    double** data_in    = fn::alloc_2D_double(dim_in);
    double*** data_out  = fn::alloc_3D_double(dim_out);


    while (true) {
        MPI_Recv(&data_in[0][0], dim_in[0]*dim_in[1], MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == STOP) {break;}

        if (status.MPI_TAG == COMPUTE) {

            for (uint i=0; i<dim_out[0]; i++) {
                for (uint j=0; j<dim_out[1]; j++) {
                    for (uint k=0; k<dim_out[2]; k++) {
                        data_out[i][j][k] = 0.0;
                    }
                }
            }


        }
        
        MPI_Send(&(data_out[0][0][0]), dim_out[0]*dim_out[1]*dim_out[2], MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
    }
}

void master(std::vector<size_t> comm_dim, int n_workers, ArgumentParser* p) {
    // MPI status flag holder
    MPI_Status status;

    if (is_valid(p)) {
        std::string dkey;
        bool skip;
        int slave;
        int shift;
        int s                   = p->i("step_first");
        int step_last           = p->i("step_last");

        H5::File* drain_file    = new H5::File(p->s("drain_file"), H5::File::ReadOnly);
        H5::File* spec_file     = new H5::File(p->s("outdir") + "/spectrum.h5", H5::File::Overwrite);

        std::vector<size_t> dim_out = {25, 157};
        std::vector<size_t> dim_in  = {25, 157, 1000};

        double** data_out   = fn::alloc_2D_double(dim_out);
        double*** data_in   = fn::alloc_3D_double(dim_in);
        



        while (true) {
            // out
            for (slave = 1; slave <= n_workers; slave++) {
                shift   = slave - 1;
                dkey    = "data_" + std::to_string(s + shift);
                skip    = !drain_file->exist(dkey) || s + shift > step_last;

                if (!skip) {
                    drain_file->getDataSet(dkey).read((double**) data_out[0]);
                }

                MPI_Send(&data_out[0][0], dim_out[0]*dim_out[1], MPI_DOUBLE, slave, (skip) ? SKIP : COMPUTE, MPI_COMM_WORLD); 
            }
            

            // in
            for (slave = 1; slave <= n_workers; slave++) {
                shift   = slave - 1;
                dkey    = "data_" + std::to_string(s + shift);

                std::cout << s+shift << std::endl;

                //std::cout << "ahoj" << std::endl;

                MPI_Recv(&data_in[0][0][0], dim_in[0]*dim_in[1]*dim_in[2], MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                if (status.MPI_TAG == COMPUTE) {
                    fn::writeDataSet(spec_file, data_in, dim_in, dkey); // initial state save
                }
            }

            s += n_workers;

            if (s > step_last) {break;}
        }
    }

    terminate(comm_dim, n_workers);
}

}

#endif

