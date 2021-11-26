#ifndef FUNCTIONS_CPP
#define FUNCTINOS_CPP

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

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

// alloc 3D double pointer array
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

#endif

