#ifndef FUNCTIONS_CPP
#define FUNCTINOS_CPP

#include <filesystem>

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

// 'Empty' simulation grid initalization
void sim_grid_init(double*** data, std::vector<size_t> dim, double r_in, double r_out, double dx) {
    for (uint i=0; i< dim[0]; i++) { // rings
        for (uint j=0; j<dim[1]; j++){ // cells
            //k = i * jdim + j; // serialized communication data coordinate
            data[i][j][0] = i; // ring coordinate
            data[i][j][1] = j; // cell coordinate
            
            data[i][j][2] = 0.0; // t
            data[i][j][3] = 0.0; // z
            data[i][j][4] = 0.0; // v
            data[i][j][5] = 0.0; // m
            data[i][j][6] = 0.0; // dm

            data[i][j][7] = std::pow(r_out, 2.0) / std::pow(r_in + (dim[0] - i - 1) * (r_out - r_in) / (dim[0] - 1), 2.0); // r (ring specific radius)
            
            data[i][j][8] = dx; // dx (inner MSMM parameter)

            data[i][j][9] = 2.0 * M_PI * ((double) j) / ((double) dim[1]); // azimuth (cell specific rotation angle)

            data[i][j][10] = 1.0; // compute flag (0.0 --> do not compute)
        }
    }
}

#endif

