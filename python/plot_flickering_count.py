#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from pathlib import Path

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern"],
    "font.size": 10
})

def plot(input: str, output: str) -> None:
    data = np.load(input) 
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(7, 7)
    ax.xaxis.set_tick_params(width=1, length=7, direction="in")
    ax.yaxis.set_tick_params(width=1, length=7, direction="in", which="both")
    #  ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.set_yscale("log")
    ax.set_xlabel("layer index $i$")
    ax.set_ylabel("$N_\mathrm{f}\ [\mathrm{events} \cdot \mathrm{step}^{-1}]$")
    ax.plot(data, color="black", linewidth=1, marker="o", markersize=3)
    fig.savefig(output, format=Path(output).suffix[1:])

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--input", type=str)
    parser.add_argument("--output", type=str)
    args = parser.parse_args()

    plot(args.input, args.output)

if __name__ == "__main__":
    main()

