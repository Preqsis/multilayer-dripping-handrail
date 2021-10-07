#!/usr/bin/env python

import h5py
import argparse
import sys

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, help="HDF5 input data file")
    p.add_argument("--output", type=str, help="HDF5 output data file")
    p.add_argument("--dataset", type=str, help="HDF5 dataset name to extact")
    p.add_argument("--list", action="store_true", default=False, help="If set, lists all input file datasets")
    args = p.parse_args()

    if args.input is None:
        print("Input file not specified!")
        sys.exit(0)

    in_file = h5py.File(args.input, "r")

    if args.list:
        for key in in_file.keys():
            print(key)
        sys.exit()

    if args.output is None:
        print("Output file not specified!")
        sys.exit(0)

    if args.dataset is None:
        print("Dataset not specified!")
        sys.exit(0)

    out_file = h5py.File(args.output, "w")
    out_file.create_dataset(args.dataset, data=in_file[args.dataset])

    for key in in_file.attrs.keys():
        out_file.attrs[key] = 1 if key == "n" else in_file.attrs[key]
    out_file.attrs["note"] = "EXTRACTED DATASET"

    in_file.close()
    out_file.close()
