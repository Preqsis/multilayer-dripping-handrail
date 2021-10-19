#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("./render")

import argparse
import math
import cairo
import h5py

import matplotlib.pyplot as plt

import multiprocessing as mp

import random

from Render import render_disc

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=str, help="HDF5 data file")
    p.add_argument("--output", default="/tmp", type=str, help="Plots output dir.")
    p.add_argument("--clear", action="store_true", help="If set => clear plots output directory")
    p.add_argument("--idim", type=int, default=50, help="")
    p.add_argument("--jdim", type=int, default=314, help="")
    p.add_argument("--width", type=int, default=1500, help="")
    p.add_argument("--height", type=int, default=1500, help="")
    p.add_argument("--first_frame", type=int, help="First frame to plot")
    p.add_argument("--last_frame", type=int, help="Last frame to plot")
    p.add_argument("--nth_frame", type=int, default=1, help="Plot every nth frame (default=100)")
    args = p.parse_args()

    w = args.width
    h = args.height
    idim = args.idim
    jdim = args.jdim

    r_in = 0.1
    r_out = 0.45

    # data
    f =  h5py.File(args.data, "r")
    #keys, data = list(f.keys()), {}
    #first_frame = 0 if args.first_frame is None else args.first_frame
    #last_frame = len(keys) if args.last_frame is None else args.last_frame+1

    first_frame = args.first_frame
    last_frame = args.last_frame

    mlimit = 16.

    for frame in range(first_frame, last_frame, args.nth_frame):

        #if "data_%d" % frame in keys:
        data = f["data_%d" % frame][()]
        #else:
        #    continue
        
        imdata = render_disc(data, idim, jdim, w, h, "magma", r_in, r_out, mlimit)
        
        with open("%s/frame_%d.png" % (args.output, frame), "wb") as f_out:
            f_out.write(imdata)
