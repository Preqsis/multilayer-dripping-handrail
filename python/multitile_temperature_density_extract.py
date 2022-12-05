import numpy as np
from argparse import ArgumentParser, Action
from h5py import File

class parse_csvarg(Action):
    def __call__(self, parser, args, values, option_string=None):
        extract_list = []
        for value in values.split(","):
            fpath, frame = value.split(":")
            extract_list.append({"fpath": fpath, "frame": frame})
        setattr(args, self.dest, extract_list)

def extract(fpath: str, frame: str):
    with File(fpath, "r") as f:
        return f[frame][()]

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--extract", action=parse_csvarg, default=[])
    parser.add_argument("--output", type=str)
    args = parser.parse_args()

    with File(args.output, "w") as f:
        for e in args.extract:
            data = extract(e["fpath"], e["frame"])
            f.create_dataset(e["frame"], data=data)

if __name__ == "__main__":
    main()
