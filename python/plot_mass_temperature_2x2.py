#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("./render")

import argparse
import h5py
import os

import numpy as np
from Render import *
from io import BytesIO
import PIL
import matplotlib.pyplot as plt
import matplotlib as mpl

class parse_csvarg(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, [v for v in values.split(",")])

def plot(data, m_range, T_range, idim, jdim, outfile, width=1500, height=1500, dpi=300, qs=1.0):
    r_in, r_out = 0.1, 0.45
    surfaces = {}

    # mass distribution
    surfaces["mass"] = render_disc(
        data,
        idim, jdim,
        r_in=r_in, r_out=r_out,
        w=width, h=height,
        return_surface=True,
        bg_rgb=(1,1,1),
        val_index=5,
        val_range=m_range,
        skip_empty=True,
        scf=qs
    )
    
    # temprature distribution
    surfaces["temperature"] = render_disc(
        data,
        idim, jdim,
        r_in=r_in, r_out=r_out,
        w=width, h=height,
        return_surface=True,
        bg_rgb=(1,1,1),
        val_index=10,
        val_range=T_range,
        cmap="jet",
        skip_empty=True,
        log_norm=True
    )

    fig, ax = plt.subplots(1, 2) 
    fig.set_size_inches(21, 9)
    plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)

    buff = BytesIO()
    surfaces["mass"].write_to_png(buff)
    ms = np.array(PIL.Image.open(buff))
    
    buff = BytesIO()
    surfaces["temperature"].write_to_png(buff)
    ts = np.array(PIL.Image.open(buff))
    
    ms_img = ax[0].imshow(ms, 
                            cmap="inferno", 
                            norm=mpl.colors.Normalize(vmin=m_range[0]*qs, vmax=m_range[1]*qs)
                            )
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    #ax[0].title.set_text("a)")
    
    ts_img = ax[1].imshow(ts, cmap="jet", norm=mpl.colors.LogNorm(vmin=T_range[0], vmax=T_range[1]))
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    #ax[1].title.set_text("b)")

    cb1 = fig.colorbar(ms_img, ax=[ax[0]], location="right", label="Cell mass $(g)$", shrink=0.9)
    cb2 = fig.colorbar(ts_img, ax=[ax[1]], location="right", label="Cell temperature $(K)$", shrink=0.9)

    fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.close(fig)

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sim_file", type=str, help="HDF5 simulation file")
    p.add_argument("--dkeys", action=parse_csvarg, default=[], help="")
    p.add_argument("--output", default="./plot", type=str, help="")
    p.add_argument("--width", type=int, default=1500, help="")
    p.add_argument("--height", type=int, default=1500, help="")
    p.add_argument("--dpi", type=int, default=300, help="")

    args = p.parse_args()

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    data = {}
    with h5py.File(args.sim_file, "r") as f:
        for dkey in args.dkeys:
            idim, jdim, qs = f.attrs["idim"], f.attrs["jdim"], f.attrs["qs"]
            data[dkey] = f[dkey][()]

    for dkey, d in data.items():
        m_range = (1e30, 0)
        T_range = (1e30, 0)

        m_min, m_max = d[:,:,5].min(), d[:,:,5].max()
        m_range = (m_range[0] if m_range[0] < m_min else m_min, m_range[1] if m_range[1] > m_max else m_max)
        
        T_min, T_max = d[:,:,10].min(), d[:,:,10].max()
        T_range = (T_range[0] if T_range[0] < T_min else T_min, T_range[1] if T_range[1] > T_max else T_max)
        T_range = (T_range[0] if T_range[0] > 0. else 10, T_range[1])

    for dkey, d in data.items():
        outfile = f"{args.output}/plot_{dkey}.png"
        plot(d, m_range, T_range, idim, jdim, outfile, width=args.width, height=args.height, dpi=args.dpi, qs=qs)

if __name__ == "__main__":
    main()

