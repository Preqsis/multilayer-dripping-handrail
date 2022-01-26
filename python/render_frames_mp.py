#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("./render")

import argparse
import cairo
import h5py
import os

import numpy as np
import multiprocessing as mp

from Render import *

def worker(sim_file, obs_file, output, frames, w, h, i, lcdepth=200, frange=(0, 100)) -> None:
    with h5py.File(sim_file, "r") as f_sim, h5py.File(obs_file, "r") as f_obs:
        idim, jdim = f_sim.attrs["idim"], f_sim.attrs["jdim"]

        data_obs = f_obs["data"][()]
        y_max = data_obs[(data_obs[:,0] >= frange[0]) & (data_obs[:,0] <= frange[1])][:,2].max()
        lc_ylim = (- y_max * 0.05, y_max * 1.05)

        print(data_obs)

        return 

        for frame in frames:
            dkey = f"data_{frame}"
            print(i, dkey)

            m = (data_obs[:,0] <= frame) & (data_obs[:,0] > frame - lcdepth)

            img = render_frame(f_sim[dkey][()], data_obs[m], idim, jdim, w=w, h=h, lc_depth=200, lc_ylim=lc_ylim)
            img.save(f"{output}/frame_{frame:06}.png")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--sim_file", type=str, help="HDF5 simulation file")
    p.add_argument("--rad_file", type=str, help="HDF5 radiation file")
    p.add_argument("--obs_file", type=str, help="HDF5 observation file")

    p.add_argument("--output", default="/frames", type=str, help="Plots output dir.")
    p.add_argument("--clear", action="store_true", help="If set => clear plots output directory")
    p.add_argument("--width", type=int, default=1920, help="")
    p.add_argument("--height", type=int, default=1080, help="")
    p.add_argument("--first_frame", default=0, type=int, help="First frame to plot")
    p.add_argument("--last_frame", default=100, type=int, help="Last frame to plot")
    p.add_argument("--s", type=int, default=1, help="Plot every nth frame (default=1)")
    p.add_argument("--n", type=int, default=8, help="number of workers (default=8)")
    p.add_argument("--lcdepth", type=int, default=200, help="Light curve plot depth (default=200)")
    args = p.parse_args()

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    if args.clear:
        os.system(f"rm -rf {args.output}/*")

    # frame-sety pro workery
    fsets = np.array_split(np.arange(args.first_frame, args.last_frame+args.s, args.s), args.n)

    # definice procesu
    proc = [mp.Process(target=worker, args=(args.sim_file, args.obs_file, args.output, frames, args.width, args.height, i, args.lcdepth, (args.first_frame, args.last_frame))) for i, frames in enumerate(fsets)]

    # zpracovani
    for p in proc:
        p.start()
    for p in proc:
        p.join()
