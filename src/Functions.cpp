#include <filesystem>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include "Functions.hpp"

typedef boost::array<double, 2> state_type;
typedef unsigned int uint;

double* Functions::alloc_1D_double(int idim) {
    double* arr;
    arr = (double*)calloc(idim, sizeof(double));
    return arr;
}

double* Functions::alloc_1D_double(std::vector<size_t> dim) {
    return alloc_1D_double(dim[0]);
}

// alloc 2D double pointer array
double** Functions::alloc_2D_double(int idim, int jdim) {
    double** arr;
    arr     = (double**)calloc(idim, sizeof(double*));
    arr[0]  = (double*)calloc(idim*jdim, sizeof(double));
    for (int i = 0; i < idim; i++) {
        arr[i] = arr[0] + i * jdim;
    }
    return arr;
}

// alloc 2D double pointer array
double** Functions::alloc_2D_double(std::vector<size_t> dim) {
    return alloc_2D_double(dim[0], dim[1]);
}

// alloc 3D double pointer array
double*** Functions::alloc_3D_double(size_t idim, size_t jdim, size_t kdim) {
    double *data = new double [idim * jdim * kdim];
    double ***array = new double **[idim];
    for (size_t i=0; i<idim; i++) {
        array[i] = new double *[jdim];
        for (size_t j = 0; j < jdim; j++) {
            array[i][j] = &(data[(i*jdim+j)*kdim]);
        }
    }
    return array;
}

// alloc 3D double pointer array
double*** Functions::alloc_3D_double(std::vector<size_t> dim) {
    return alloc_3D_double(dim[0], dim[1], dim[2]);
}

// alloc 4D double pointer array
double**** Functions::alloc_4D_double(size_t idim, size_t jdim, size_t kdim, size_t ldim) {
    double *data = new double [idim * jdim * kdim * ldim];
    double ****array = new double ***[idim];
    for (size_t i=0; i<idim; i++) {
        array[i] = new double **[jdim];
        for (size_t j = 0; j < jdim; j++) {
            array[i][j] = new double *[kdim];
            for (size_t k = 0; k < kdim; k++) {
                array[i][j][k] = &(data[((i*jdim+j)*kdim+k)*ldim]);
            }
        }
    }
    return array;
}

// alloc 4D double pointer array
double**** Functions::alloc_4D_double(std::vector<size_t> dim) {
    return alloc_4D_double(dim[0], dim[1], dim[2], dim[3]);
}

void Functions::writeDataSet(H5::File* file, double** data, std::vector<size_t> dim, std::string key) {
    if (!file->exist(key)) {
        H5::DataSet* ds = new H5::DataSet(file->createDataSet<double>(key, H5::DataSpace(dim)));
        ds->write((double**) data[0]);
        delete ds;
    }
}

void Functions::writeDataSet(H5::File* file, double*** data, std::vector<size_t> dim, std::string key) {
    if (!file->exist(key)) {
        H5::DataSet* ds = new H5::DataSet(file->createDataSet<double>(key, H5::DataSpace(dim)));
        ds->write((double***) data[0][0]);
        delete ds;
    }
}

void Functions::writeDataSet(H5::File* file, double**** data, std::vector<size_t> dim, std::string key) {
    if (!file->exist(key)) {
        H5::DataSet* ds = new H5::DataSet(file->createDataSet<double>(key, H5::DataSpace(dim)));
        ds->write((double****) data[0][0][0]);
        delete ds;
    }
}

// check if dir path exists
bool Functions::isdir(std::string path) {
    return std::filesystem::is_directory(path);
}

// create dir if not exists
void Functions::mkdir(std::string path) {
    // if not exist create outdir
    if (!isdir(path)) {
        std::filesystem::create_directory(path);
    }
}
