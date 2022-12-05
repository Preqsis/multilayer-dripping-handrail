#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from h5py import File
from matplotlib.ticker import FormatStrFormatter
from pathlib import Path
from argparse import ArgumentParser, Action

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

class parse_csvarg(Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, [v for v in values.split(",")])

def plot(obs_file: str, plot_file: str, filters: list[str]) -> None:
    with File(obs_file, "r") as f:
        data = f["data"][()]

    fig, axes = plt.subplots(len(filters), 1) 
    fig.set_size_inches(7, len(filters) * 3.5)

    for i, filter in enumerate(filters):
        ax = axes[i] if len(filters) > 1 else axes 
        ax.xaxis.set_tick_params(width=1, length=7, direction="in")
        ax.yaxis.set_tick_params(width=1, length=7, direction="in", which="both")
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        ax.set_ylabel(r"$L_{\mathrm{bol}}\ \mathrm{[erg \cdot s^{-1}]}$")
        if i < len(filters)-1:
            ax.axes.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel(r"$t\ \mathrm{[s]}$")
        ax.set_title(r'$\mathcal{%s}$ filter' % (filter))
        ax.plot(data[:,0] * 60, data[:,fmap[filter]], color="black", linewidth=1)
    fig.savefig(plot_file, bbox_inches="tight", format=Path(plot_file).suffix[1:])

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--obs_file", type=str, help="")
    parser.add_argument("--plot_file", default="./plot", type=str, help="")
    parser.add_argument("--filters", action=parse_csvarg, default=[], help="")
    args = parser.parse_args()

    plot(args.obs_file, args.plot_file, args.filters)

if __name__ == "__main__":
    main()

