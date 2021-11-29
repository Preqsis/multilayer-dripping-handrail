#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP

#include <iostream>
#include <filesystem>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

typedef boost::array<double, 2> state_type;
typedef long unsigned int luint;
typedef unsigned int uint;

const double M_SUN  = 1.9891e30;
const double G      = 6.6743e-11;

// Compute flags
const int MASTER    = 0;
const int COMPUTE   = 1;
const int STOP      = 2;
const int SKIP      = 3;

namespace Functions {

// alloc 2D double pointer array
double** alloc_2D_double(int idim, int jdim) {
    double** arr;
    arr     = (double**)calloc(idim, sizeof(double*));
    arr[0]  = (double*)calloc(idim*jdim, sizeof(double));
    for (int i = 0; i < idim; i++) {
        arr[i] = arr[0] + i * jdim;
    }
    return arr;
}

// alloc 2D double pointer array
double** alloc_2D_double(std::vector<size_t> dim) {
    return alloc_2D_double(dim[0], dim[1]);
}

double ***alloc_3D_double(int l, int m, int n) {
    double *data = new double [l*m*n];
    double ***array = new double **[l];
    for (int i=0; i<l; i++) {
        array[i] = new double *[m];
        for (int j=0; j<m; j++) {
            array[i][j] = &(data[(i*m+j)*n]);
        }
    }
    return array;
}

// alloc 3D double pointer array
double*** alloc_3D_double(std::vector<size_t> dim) {
    return alloc_3D_double(dim[0], dim[1], dim[2]);
}

void writeDataSet(H5::File* file, double** data, std::vector<size_t> dim, std::string key) {
    if (!file->exist(key)) {
        H5::DataSet* ds = new H5::DataSet(file->createDataSet<double>(key, H5::DataSpace(dim)));
        ds->write((double**) data[0]);
        delete ds;
    }
}

void writeDataSet(H5::File* file, double*** data, std::vector<size_t> dim, std::string key) {
    if (!file->exist(key)) {
        H5::DataSet* ds = new H5::DataSet(file->createDataSet<double>(key, H5::DataSpace(dim)));
        ds->write((double***) data[0][0]);
        delete ds;
    }
}

// check if dir path exists
bool isdir(std::string path) {
    return std::filesystem::is_directory(path);
}

// create dir if not exists
void mkdir(std::string path) {
    // if not exist create outdir
    if (!isdir(path)) {
        std::filesystem::create_directory(path);
    }
}

}

#endif

