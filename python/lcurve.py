#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import argparse

from scipy import constants as cs

M_sun = 1.9891E30 # hmostnost slunce [kg]
G = 6.67430E-11 # grav. konstanta
sigma = 5.670374419e-8 # Stefan-Boltzmann

M_primary = 1. * M_sun # kg ... white dwarf
dM = 10**14 # kg * s^-1 ... rychlost akrece

M_secondary = 0.5 * M_sun

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, help="HDF5 data file")
    args = p.parse_args()

    f = h5py.File(args.input, "r")
    idim, jdim = f.attrs["idim"], f.attrs["jdim"]
    data_drain = f["drain"][()]
    f.close()

    q = 1e14 * dt
    r_in = 9e6 # m
    r_out = 50. * r_in
    T = 1e6 # teplota atm. disku ??
    wl = 200e-9 # 2000 angstrom v m
    ny = cs.c / wl

    #r_out = d * (M_secondary / (3. * (M_primary+M_secondary)))**(1./3.) # vjensi polomer disku (v L1)






    i = np.linspace(0, idim-1, idim) # pole indexu prstencu
    r_sim = r_in + (idim - i - 1) * (r_out - r_in) / (idim - 1) # hodnoty vsech polomeru [m]

    dt = np.sqrt(4. * np.pi**2 * r_sim[0]**3 / (G * M_primary)) / jdim
    dR = r_sim[0] - r_sim[1]


    r = np.linspace(r_in, r_out, int(1e5))


    T_in        = (3. * G * M_primary * dM / (8. * np.pi * r_in**3. * sigma))**(1./4.)
    T_ef        = T_in * (r_in / r)**(3./4.) * (1. - np.sqrt(r_in / r))**(1./4.)
    T_ef_sim    = T_in * (r_in / r_sim)**(3./4.) * (1. - np.sqrt(r_in / r_sim))
    T_ef_sim    = T_ef_sim[:-1]


    e = np.zeros((data_drain.shape[0], data_drain.shape[1]-1), dtype="<f8")
    for idx in range(idim-1):
        e[:,idx] = q * 0.5 * G * M_primary * data_drain[:,idx] * dR / (r_sim[idx] * r_sim[idx+1])

    R = r_sim[:-1]

    E = np.sum(np.abs(e), axis=1)

    dL = 4.9e-11 * (E / T) * np.exp(-1. * ny * cs.h / (cs.k * T))

    t = np.linspace(0., (dL.shape[0]-1)*dt, dL.shape[0])



    B   = (2. * cs.h * ny**3. / cs.c**2.) / (np.exp(cs.h * ny / (cs.k * T_ef_sim)) - 1.)
    L1  = np.sum(2. * 2. * np.pi * R * dR * B) + dL


    data_lc         = np.empty((L1.shape[0], 2), dtype=np.float64)
    data_lc[:,0]    = np.linspace(1, data_lc.shape[0], data_lc.shape[0])
    data_lc[:,1]    = L1

    f = h5py.File(fname, "a")
    f.create_dataset("data_lc", data=data_lc)
    f.close()
    

if __name__ == "__main__":
    main()
