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

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern"],
    "font.size": 10})

R_sun = 69634000000

class parse_csvarg(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, [v for v in values.split(",")])

class parse_range(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, [float(val) for val in values.split(",")[:2]])

def plot(data, density_range, T_range, idim, jdim, outfile, width=1500, height=1500, par=None, format="png") -> None:
    r_in, r_out = 0.1, 0.45
    surfaces = {}

    # mass distribution
    surfaces["density"] = render_disc(
        data,
        idim, jdim,
        r_in=r_in, r_out=r_out,
        w=width, h=height,
        return_surface=True,
        bg_rgb=(1,1,1),
        val_index=5,
        val_range=density_range,
        skip_empty=True
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
        cmap="plasma",
        skip_empty=True,
        log_norm=True
    )

    fig, ax = plt.subplots(1, 2) 
    fig.set_size_inches(7, 3)
    #plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)

    buff = BytesIO()
    surfaces["density"].write_to_png(buff)
    ms = np.array(PIL.Image.open(buff))
    
    buff = BytesIO()
    surfaces["temperature"].write_to_png(buff)
    ts = np.array(PIL.Image.open(buff))
    
    ms_img = ax[0].imshow(ms, 
                            cmap="inferno", 
                            norm=mpl.colors.Normalize(vmin=density_range[0], vmax=density_range[1])
                            )
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].axis("off")
    #ax[0].title.set_text("a)")
    
    ts_img = ax[1].imshow(ts, cmap="plasma", norm=mpl.colors.LogNorm(vmin=T_range[0], vmax=T_range[1]))
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].axis("off")
    #ax[1].title.set_text("b)")

    cb1 = fig.colorbar(ms_img, ax=[ax[0]], location="right", shrink=0.8)
    cb1.set_label(label="$\Sigma\ \mathrm{[g \cdot cm^{-2}]}$")
    cb2 = fig.colorbar(ts_img, ax=[ax[1]], location="right", shrink=0.8)
    cb2.set_label(label="$T\ \mathrm{[K]}$")

    fig.savefig(f"{outfile}.{format}", bbox_inches="tight", dpi=600, format=format)
    plt.close(fig)

def get_range(data, val_index, force_min=None, force_max=None) -> List:
    rmin, rmax = 1e30, 0.0

    for dkey, d in data.items():
        vmin, vmax = d[:,:,val_index].min(), d[:,:,val_index].max()
        rmin = rmin if rmin < vmin else vmin 
        rmax = rmax if rmax > vmax else vmax

    rmin = force_min if force_min is not None and force_min < rmin else rmin
    rmax = force_max if force_max is not None and force_max > rmax else rmax

    return [rmin, rmax]

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sim_file", type=str, help="HDF5 simulation file")
    p.add_argument("--dkeys", action=parse_csvarg, default=[], help="")
    p.add_argument("--output", default="./plot", type=str, help="")
    p.add_argument("--temperature_range", action=parse_range, help="")
    p.add_argument("--density_range", action=parse_range, help="")

    args = p.parse_args()

    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    
    data, par = {}, {}
    with h5py.File(args.sim_file, "r") as f:
        for key in f.attrs:
            par[key] = f.attrs[key]


        idim, jdim, qs = par["idim"], par["jdim"], par["qs"]
        r_in, r_out = par["r_in"] * R_sun, par["r_out"] * R_sun

        # pro prepocet hmotnosit na plosnou hustotu
        S   = np.zeros((idim, jdim))
        idx = np.arange(0, idim, 1)
        r   = r_in + (idim - idx - 1) * (r_out - r_in) / (idim - 1) 
        for i in idx:
            S[i,:] = r[i] * 2. * np.pi * (r[0] - r[1]) / jdim

        for dkey in args.dkeys:
            data[dkey] = f[dkey][()]
            data[dkey][:,:,5] = data[dkey][:,:,5] * qs / S


    density_range = args.density_range if args.density_range is not None else get_range(data, 5)
    temperature_range = args.temperature_range if args.temperature_range is not None else get_range(data, 10, force_min=1.0)
    print(f"Temperature range: {temperature_range}; Density range: {density_range}")

    for dkey, d in data.items():
        print(dkey)
        plot(d, density_range, temperature_range, idim, jdim, f"{args.output}/plot_{dkey}", width=1500, height=1500, par=par, format="pdf")

if __name__ == "__main__":
    main()

