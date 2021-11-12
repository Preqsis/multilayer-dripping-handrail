#!/usr/bin/env python

import h5py
import argparse

def listAttrs(fpath) -> None:
    with h5py.File(fpath, "r") as f:
        for key in f.keys():
            print(key)

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, required=True, help="HDF5 input data file")
    p.add_argument("--output", type=str, help="HDF5 output data file")
    p.add_argument("--dataset", type=str, help="HDF5 dataset name to extact")
    p.add_argument("--list", action="store_true", default=False, help="If set, lists all input file datasets")
    args = p.parse_args()

    if args.list:
        listAttrs(args.input)
        return

    if len([x for x in (args.output, args.dataset) if x is not None]) <= 1:
        p.error('both --output and --dataset must be specified')
    

    return

    with h5py.File(args.input, "r") as fin, h5py.File(args.output, "w") as fout:
        fout.create_dataset(args.dataset, data=fin[args.dataset])
        for key in fin.attrs.keys():
            fout.attrs[key] = 1 if key == "n" else fin.attrs[key]
        fout.attrs["note"] = "EXTRACTED DATASET"

if __name__ == "__main__":
    main()
