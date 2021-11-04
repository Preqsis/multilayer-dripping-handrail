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

def worker(fpath, output, frames, w, h, i, lcdepth=200) -> None:
    data_lc         = np.empty((100, 2))
    data_lc[:,0]    = np.arange(0, data_lc.shape[0], 1)
    data_lc[:,1]    = np.random.rand(data_lc.shape[0], 1)[:,0]

    with h5py.File(fpath, "r") as f:
        idim, jdim = f.attrs["idim"], f.attrs["jdim"]

        data_lc = f["data_lc"][()]

        for frame in frames:
            dkey = f"data_{frame}"
            print(i, dkey)

            m = (data_lc[:,0] <= frame) & (data_lc[:,0] > frame - lcdepth)

            img = render_frame(f[dkey][()], data_lc[m], idim, jdim, w=w, h=h)
            img.save(f"{output}/frame_{frame:06}.png")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, help="HDF5 data file")
    p.add_argument("--output", default="/tmp", type=str, help="Plots output dir.")
    p.add_argument("--clear", action="store_true", help="If set => clear plots output directory")
    p.add_argument("--width", type=int, default=1920, help="")
    p.add_argument("--height", type=int, default=1080, help="")
    p.add_argument("--first_frame", default=0, type=int, help="First frame to plot")
    p.add_argument("--last_frame", default=100, type=int, help="Last frame to plot")
    p.add_argument("--s", type=int, default=1, help="Plot every nth frame (default=1)")
    p.add_argument("--n", type=int, default=8, help="number of workers (default=8)")
    p.add_argument("--lcdepth", type=int, default=200, help="Light curve plot depth (default=200)")
    args = p.parse_args()

    # frame-sety pro workery
    fsets = np.array_split(np.arange(args.first_frame, args.last_frame+args.s, args.s), args.n)

    # definice procesu
    proc = [mp.Process(target=worker, args=(args.input, args.output, frames, args.width, args.height, i, args.lcdepth)) for i, frames in enumerate(fsets)]

    # zpracovani
    for p in proc:
        p.start()
    for p in proc:
        p.join()
