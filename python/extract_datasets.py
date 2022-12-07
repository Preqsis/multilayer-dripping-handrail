import json
import numpy as np
from argparse import ArgumentParser, Action
from h5py import File

def extract(fpath: str, frame: str):
    with File(fpath, "r") as f:

        for key, value in f.attrs.items():
            print(key, value)

        return f[frame][()]

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--conf", type=str)
    parser.add_argument("--output", type=str)
    args = parser.parse_args()

    with open(args.conf, "r") as f:
        conf = json.load(f)

    with File(args.output, "w") as f:
        for e in conf["extract"]:
            print(e)
            data = extract(e["fpath"], e["frame"])
            f.create_dataset(e["frame_out"], data=data)

if __name__ == "__main__":
    main()
