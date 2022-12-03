import numpy as np
import multiprocessing as mp
from argparse import ArgumentParser
from h5py import File

def power_sum_worker(proc_num: int, input_file: File, steps: list[int], return_dict: dict) -> None:
    I, J, K, L = input_file[f"d0"][()].shape
    data_power = np.zeros((I,))
    for step in steps:
        print(f"power_sum_worker: {proc_num}, step: {step}")
        data = input_file[f"d{step}"][()]
        data = data[:,:,:,1]
        data = np.sum(data, axis=2)
        data = np.sum(data, axis=1)
        data_power += data
    return_dict[proc_num] = data_power

def compute(rad_file: str, data_file: str, n: int = 4) -> None:
    with File(rad_file, "r") as f:
        n_steps = len(f.keys())
        wsets = np.array_split(np.arange(0, n_steps, 1), n)
        manager = mp.Manager()
        return_dict = manager.dict()
        proc = [mp.Process(target=power_sum_worker, args=(i, f, wset, return_dict)) for i, wset in enumerate(wsets)]
        for p in proc:
            p.start()
        for p in proc:
            p.join()

    data = np.zeros_like(return_dict[0])
    for i in range(args.n):
        data += return_dict[i]
    data /= n_steps
    np.save(data_file, data)

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--rad_file", type=str)
    parser.add_argument("--data_file", type=str)
    parser.add_argument("--n", type=int)
    args = parser.parse_args()

    compute(args.rad_file, args.data_file, n)

if __name__ == "__main__":
    main()
