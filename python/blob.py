#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("./render")

import numpy as np
import h5py
import argparse
import PIL

from Render import *

class parse_azimuth(argparse.Action):
    """Azimuth parser class"""
    def __call__(self, parser, args, value, option_strings=None):
        setattr(args, self.dest, value * 3.14159265359 / 180.)

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, help="HDF5 data file")
    p.add_argument("--output", type=str, help="HDF5 data file")
    p.add_argument("--dkey", type=str, help="Blob target frame key")
    p.add_argument("--ic", default=0, type=int)
    p.add_argument("--ac", default=0., action=parse_azimuth, type=float)
    p.add_argument("--r", default=200, type=float)
    p.add_argument("--m", default=10., type=float)
    args = p.parse_args()

    # nedefinovany datovy klic
    if args.dkey is None:
        return None

    # pridat blob podle zadanych parametru
    with h5py.File(args.input, "r+") as f:
        idim, jdim = f.attrs["idim"], f.attrs["jdim"]

        # data prislusneho kroku
        data    = f[args.dkey][()]

        # "stredova bunka" blobu
        m       = data[:,0] == args.ic
        jc      = np.argmin(np.abs(data[m][data[m,1].astype(np.int),9] - args.ac))

        # streadova bunka
        c       = data[jdim * args.ic + jc]
        rc      = idim - c[0]

        # ostatni bunky
        r       = idim - data[:,0]
        azm     = data[:,9]

        # vytvorit masku bunek v danem okruhu
        lx      = r * np.cos(azm) - rc * np.cos(c[9])
        ly      = r * np.sin(azm) - rc * np.sin(c[9])
        l       = np.sqrt(lx**2. + ly**2.)
        m_blob  = l <= args.r

        # pricist
        data[m_blob,5] += np.floor(args.m / (0.8 * l[m_blob]))
                
        img = render_disc(data, idim, jdim, return_surface=True)
        img.write_to_png("/home/preqsis/Plocha/test.png")

        d       = f[args.dkey]
        d[...]  = data

if __name__ == "__main__":
    main()

