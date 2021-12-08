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
    double a = 2.0 * cs::h * std::pow(cs::c, 2.0);
    double b = cs::h * cs::c / (wl * cs::k * T);
    return a / (std::pow(wl, 5.0) * (std::exp(b) - 1.0 ));
}

void Radiation::slave(std::vector<size_t> dim_mass, std::vector<size_t> dim_spec, ArgumentParser* p) {
    /** MPI status flag **/
    MPI_Status status;

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

    double dR = (p->d("r_out") - p->d("r_in")) / dim_mass[0];
    double* r = new double[dim_mass[0]+1];
    for (size_t i=0; i <= dim_mass[0]; i++) {
        r[i] = p->d("r_out") - i * dR;
    }

    double E, T, dL, B, S, L, wl, frq, drn;
    double A = cs::G * p->d("m_primary") * cs::m_sun; // Bulhar --> G * Mp = 1.0 
    double Q = 1e17;

    while (true) {
        MPI_Recv(&data_mass[0][0][0], len_mass, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == cs::mpi::STOP) {break;}

        if (status.MPI_TAG == cs::mpi::COMPUTE) {
            for (size_t i = 0; i < dim_spec[0]; i++) {          // pres vsechny prstence
                for (size_t j = 0; j < dim_spec[1]; j++) {      // pres vsechny bunky v prstenci
                    for (size_t k = 0; k < dim_spec[2]; k++) {  // pres rozsah vl. delek
                        
                        drn = Q * data_mass[i][j][9]; // drain

                        //if (drn > 0.0) {
                        //std::cout << drn << ", " << dR << std::endl;
                        //}

                        E   = 0.5 * A * dR * drn / (r[i] * r[i+1]);

                        T   = data_mass[i][j][10];

                        wl  = data_spec[i][j][k][0];
                        frq = cs::c / wl;

                        dL  = 4.9e-11 * (E / 1e5) * std::exp(-1.0 * cs::h * frq/ (cs::k * 1e5));

                        B   = planck(wl, T);
                        
                        S   = M_PI * (pow(r[i], 2.0) - pow(r[i+1], 2.0)) / dim_mass[1];

                        //L   = 2.0 * S * B + dL;
                        L   = 2.0 * S * B + dL;

                        data_spec[i][j][k][1] = L;
                    }
                }
            }

            //spectrum(data_mass, data_spec, dim_mass, dim_spec);
        }
        
        MPI_Send(&(data_spec[0][0][0][0]), len_spec, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
    }
}

void Radiation::master(std::vector<size_t> dim_mass, std::vector<size_t> dim_spec, int n_workers, ArgumentParser* p) {
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
    bool v = p->b("v");


    // input drain file / output spectrum File
    std::string mass_fpath = (p->isSet("mass_file")) ? p->s("mass_file") : p->s("outdir") + "/mass.h5";
    H5::File* mass_file = new H5::File(mass_fpath, H5::File::ReadOnly);
    H5::File* spec_file = new H5::File(p->s("outdir") + "/spectrum.h5", H5::File::Overwrite);

    // Comm length
    size_t len_mass = dim_mass[0] * dim_mass[1] * dim_mass[2];
    size_t len_spec = dim_spec[0] * dim_spec[1] * dim_spec[2] * dim_spec[3];

    // Comms data allocation
    double*** data_mass     = fn::alloc_3D_double(dim_mass);
    double**** data_spec    = fn::alloc_4D_double(dim_spec);

    while (true) {
        // out
        for (slave = 1; slave <= n_workers; slave++) {
            shift   = slave - 1;
            dkey    = "data_" + std::to_string(s + shift);
            skip    = !mass_file->exist(dkey) || s + shift > step_last;

            if (!skip) {
                mass_file->getDataSet(dkey).read((double***) data_mass[0][0]);
            }

            MPI_Send(&data_mass[0][0][0], len_mass, MPI_DOUBLE, slave, (skip) ? cs::mpi::SKIP : cs::mpi::COMPUTE, MPI_COMM_WORLD); 
        }

        // in
        for (slave = 1; slave <= n_workers; slave++) {
            shift   = slave - 1;
            dkey    = "data_" + std::to_string(s + shift);

            MPI_Recv(&data_spec[0][0][0][0], len_spec, MPI_DOUBLE, slave, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == cs::mpi::COMPUTE) {
                fn::writeDataSet(spec_file, data_spec, dim_spec, dkey); // initial state save
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
    terminate(dim_mass, n_workers);

    delete mass_file;
    delete spec_file;

    // info msg.
    if (v) {
        std::cout << std::endl << "Done ..." << std::endl;
    }
}
