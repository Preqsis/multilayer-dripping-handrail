#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

namespace Functions {

// alloc 2D double pointer array
double** alloc_2D_double(int idim, int jdim);

// alloc 2D double pointer array
double** alloc_2D_double(std::vector<size_t> dim);

// alloc 3D double pointer array
double ***alloc_3D_double(size_t idim, size_t jdim, size_t kdim);

// alloc 3D double pointer array
double*** alloc_3D_double(std::vector<size_t> dim);

// alloc 4D double pointer array
double**** alloc_4D_double(size_t idim, size_t jdim, size_t kdim, size_t ldim);

// alloc 4D double pointer array
double**** alloc_4D_double(std::vector<size_t> dim);

void writeDataSet(H5::File* file, double** data, std::vector<size_t> dim, std::string key);

void writeDataSet(H5::File* file, double*** data, std::vector<size_t> dim, std::string key);

void writeDataSet(H5::File* file, double**** data, std::vector<size_t> dim, std::string key);

// check if dir path exists
bool isdir(std::string path);

// create dir if not exists
void mkdir(std::string path);

}

#endif
