#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import cairo
import h5py

import matplotlib.pyplot as plt

import random



if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=str, help="HDF5 data file")
    p.add_argument("--output", default="/tmp", type=str, help="Plots output dir.")
    p.add_argument("--clear", action="store_true", help="If set => clear plots output directory")
    p.add_argument("--idim", type=int, default=25, help="")
    p.add_argument("--jdim", type=int, default=157, help="")
    p.add_argument("--width", type=int, default=1500, help="")
    p.add_argument("--height", type=int, default=1500, help="")
    p.add_argument("--first_frame", type=int, help="First frame to plot")
    p.add_argument("--last_frame", type=int, help="Last frame to plot")
    p.add_argument("--nth_frame", type=int, default=1, help="Plot every nth frame (default=100)")
    args = p.parse_args()

    r_in = 0.1
    r_out = 0.45
    dr = (r_out - r_in) / args.idim

    # data
    f =  h5py.File(args.data, "r")
    #keys, data = list(f.keys()), {}
    #first_frame = 0 if args.first_frame is None else args.first_frame
    #last_frame = len(keys) if args.last_frame is None else args.last_frame+1

    first_frame = args.first_frame
    last_frame = args.last_frame

    mlimit = 16.
    cmap = plt.cm.get_cmap("magma")

    for frame in range(first_frame, last_frame, args.nth_frame):

        print(frame)

        #if "data_%d" % frame in keys:
        data = f["data_%d" % frame][()]
        #else:
        #    continue

        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, args.width, args.height)
        c = cairo.Context(surface)
        c.scale(args.width, args.height)

        # set bacground color (black)
        c.set_source_rgb(0, 0, 0)
        c.paint()
        c.save()
        c.restore()

        c.set_line_width(dr)
        
        dphi = 2. * math.pi / args.jdim

        r = r_out
        for i in range(args.idim):
            for j in range(args.jdim):
                k = i * args.jdim + j

                m, azm = data[k][5], data[k][9] % (2. * math.pi)

                rgba = cmap(m / mlimit)

                c.arc(0.5, 0.5, r, azm, azm+dphi)
                c.set_source_rgb(rgba[0], rgba[1], rgba[2])  # Solid color
                c.stroke()

            r -= dr

        surface.write_to_png("%s/frame_%06d.png" % (args.output, frame))
        surface.finish()
