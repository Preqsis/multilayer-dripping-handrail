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

from Render import render_disc

def worker(fpath, output, frames, w, h, i):
    with h5py.File(fpath, "r") as f:
        idim, jdim = f.attrs["idim"], f.attrs["jdim"]

        for frame in frames:
            dkey = f"data_{frame}" 
            print(i, dkey)
            with open(f"{output}/frame_{frame:06}.png", "wb") as out:
                out.write(render_disc(f[dkey][()], idim, jdim, w=w, h=h))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, help="HDF5 data file")
    p.add_argument("--output", default="/tmp", type=str, help="Plots output dir.")
    p.add_argument("--clear", action="store_true", help="If set => clear plots output directory")
    p.add_argument("--width", type=int, default=1500, help="")
    p.add_argument("--height", type=int, default=1500, help="")
    p.add_argument("--first_frame", default=0, type=int, help="First frame to plot")
    p.add_argument("--last_frame", default=100, type=int, help="Last frame to plot")
    p.add_argument("--s", type=int, default=1, help="Plot every nth frame (default=100)")
    p.add_argument("--n", type=int, default=8, help="number of workers (default=8)")
    args = p.parse_args()

    # frame-sety pro workery
    fsets = np.array_split(np.arange(args.first_frame, args.last_frame+args.s, args.s), args.n)
    
    # definice procesu
    proc = [mp.Process(target=worker, args=(args.input, args.output, frames, args.width, args.height, i)) for i, frames in enumerate(fsets)]

    # zpracovani
    for p in proc:
        p.start()
    for p in proc:
        p.join()
        


