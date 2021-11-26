all:
	mpic++ mld_sim.cpp -std=c++17 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5 -I/usr/include/hdf5/serial -Wall -o mld_sim
