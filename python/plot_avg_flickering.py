import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from h5py import File
from pathlib import Path

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern"],
    "font.size": 10
})

def plot(power_file: str, event_file: str, plot_file: str) -> None:
    data_power, data_event = np.load(power_file), np.load(event_file) 

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches(7, 7)

    ax[0].set_title("Average event count")
    ax[0].xaxis.set_tick_params(width=1, length=7, direction="in")
    ax[0].yaxis.set_tick_params(width=1, length=7, direction="in", which="both")
    ax[0].set_yscale("log")
    ax[0].set_ylabel(r"$\bar{N}_i\ [\mathrm{events} \cdot \mathrm{step}^{-1}]$")
    ax[0].axes.xaxis.set_ticklabels([])
    ax[0].invert_xaxis()
    ax[0].plot(data_event, color="black", linewidth=1, marker="o", markersize=3)

    ax[1].set_title("Average radiation power")
    ax[1].xaxis.set_tick_params(width=1, length=7, direction="in")
    ax[1].yaxis.set_tick_params(width=1, length=7, direction="in", which="both")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("layer index $i$")
    ax[1].set_ylabel(r"$\bar{L}_{i}\ [\mathrm{erg} \cdot \mathrm{step}^{-1}]$")
    ax[1].invert_xaxis()
    ax[1].bar(np.arange(0, data_power.shape[0], 1), data_power, color="#5f9ea0", linewidth=1)

    fig.savefig(plot_file, format=Path(plot_file).suffix[1:])

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--power_file", type=str)
    parser.add_argument("--event_file", type=str)
    parser.add_argument("--plot_file", type=str)
    parser.add_argument("--n", type=int)
    args = parser.parse_args()

    plot(args.power_file, args.event_file, args.plot_file)

if __name__ == "__main__":
    main()
