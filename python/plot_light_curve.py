#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import h5py
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tikzplotlib as tz
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern"],
    "font.size": 10})

fmap = {"U": 0, "B": 1, "V": 2, "R":3, "I": 4}

def savepdf_tex(fig, name, **kwargs):
    import subprocess, os
    fig.savefig("temp.pdf", format="pdf", **kwargs)
    incmd = ["inkscape", "temp.pdf", "--export-pdf={}.pdf".format(name),
             "--export-latex"] #"--export-ignore-filters",
    subprocess.check_output(incmd)
    os.remove("temp.pdf")

class parse_csvarg(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, [v for v in values.split(",")])

def plot(data, outfile, filters, format="png") -> None:
    fig, axes = plt.subplots(len(filters), 1) 
    fig.set_size_inches(7, len(filters) * 2.5)

    for i, filter in enumerate(filters):
        print(i, filter, fmap[filter])

        ax = axes[i] if len(filters) > 1 else axes
        
        ax.xaxis.set_tick_params(width=1, length=7, direction="in")
        ax.yaxis.set_tick_params(width=1, length=7, direction="in", which="both")
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

        ax.set_ylabel("$L_{\mathrm{bol}}\ \mathrm{[erg \cdot s^{-1}]}$")

        if i < len(filters)-1:
            ax.axes.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel("$t\ \mathrm{[s]}$")

        ax.text(.5,.92, f'$\mathcal{filter}$ filter',
            horizontalalignment='center',
            transform=ax.transAxes)

        ax.plot(data[:,0] * 60, data[:,fmap[filter]], color="black", linewidth=1)

    fig.savefig(f"{outfile}.{format}", bbox_inches="tight", format=format)

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--obs_file", type=str, help="")
    p.add_argument("--filters", action=parse_csvarg, default=[], help="")
    p.add_argument("--output", default="./plot", type=str, help="")
    p.add_argument("--fname", default="light_curve", type=str, help="")
    args = p.parse_args()

    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    
    data, par = {}, {}
    with h5py.File(args.obs_file, "r") as f:
        data = f["data"][()]

    plot(data,  f"{args.output}/{args.fname}", filters=args.filters, format="pdf")

if __name__ == "__main__":
    main()

