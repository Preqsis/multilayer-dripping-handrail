#!/usr/bin/env python

import h5py
import argparse
import re

class parse_ints(argparse.Action):
    RE_INTS         = re.compile(r'(?P<num>\d+)', re.DOTALL)
    RE_INTS_RANGE   = re.compile(r'\d+:\d+')

    def __call__(self, parser, args, value, option_strings=None) -> None:
        if (self.RE_INTS_RANGE.match(value)): # zadane jako range
            setattr(args, self.dest, range(*tuple(map(int, value.split(':')))))
        else: # zadane jako list
            setattr(args, self.dest, [int(i.group("num")) for i in self.RE_INTS.finditer(value)])

def listAttrs(fpath) -> None:
    with h5py.File(fpath, "r") as f:
        for key in f.keys():
            print(key)

def extractSteps(input, output, steps) -> None:
    with h5py.File(input, "r") as fi, h5py.File(output, "w") as fo:
        for s in steps:
            dkey = f"data_{s}"
            print(dkey)
            fo.create_dataset(dkey, data=fi[dkey])

def extractDrain(input, output, steps) -> None:
    with h5py.File(input, "r") as fi, h5py.File(output, "a") as fo:
        for dkey in fi.keys():
            if not re.match("^drain_([0.9]{1,})$", dkey):
                continue
            print(dkey)
            data = fi[dkey][()]
            data = data[tuple(steps),:].shape
            fo.create_dataset(dkey, data=data)
        
def copyAttrs(input, output, add_attrs = {}) -> None:
    with h5py.File(input, "r") as fi, h5py.File(output, "r+") as fo:
        for key in fi.attrs.keys():
            fo.attrs[key] = fi.attrs[key]

        for key in add_attrs:
            fo.attrs[key] = add_attrs[key]

        fo.attrs["note"] = "EXTRACTED"

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, required=True, help="HDF5 input data file")
    p.add_argument("--output", type=str, help="HDF5 output data file")
    p.add_argument("--steps", action=parse_ints, help="")
    p.add_argument("--list", action="store_true", default=False, help="If set, lists all input file datasets")
    args = p.parse_args()

    if args.list:
        listAttrs(args.input)
        return

    if len([x for x in (args.output, args.steps) if x is not None]) <= 1:
        p.error('both --output and --steps must be specified')

    extractSteps(args.input, args.output, args.steps)

    #extractDrain(args.input, args.output, args.steps)

    copyAttrs(args.input, args.output, add_attrs={"n": len(args.steps)})

if __name__ == "__main__":
    main()
