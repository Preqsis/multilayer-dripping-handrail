LIB_DIR=./lib

INC=$(LIB_DIR)/highfive /usr/include/hdf5/serial/ /usr/local/include/boost_1_76_0 $(LIB_DIR)
INC_DIRS=$(foreach d, $(INC), -I$d)

CFLGS=-c -Wall -pedantic -fopenmp

all:
	mpic++ main.cpp -std=c++17 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5 $(INC_DIRS) -o mldrun

