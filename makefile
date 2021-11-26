LIB_DIR=./lib

INC=$(LIB_DIR)/highfive /usr/include/hdf5/serial/
INC_DIRS=$(foreach d, $(INC), -I$d)

all:
	mpic++ mld_sim.cpp -std=c++17 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5 $(INC_DIRS) -Wall -o mld_sim

