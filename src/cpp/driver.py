"""driver.py

"""

from joblib import Parallel, delayed

import numpy as np
import itertools as it
import multiprocessing as mp
import subprocess


def parameters(n, d, lam):
    """Given (n, d, lam), try to return an a, b consistent.

    """
    b = d - lam * np.sqrt(d)
    a = 2.0 * lam * np.sqrt(d) + b

    if a < 0.0 or a > n:
        raise ValueError("a is not in range")
    if b < 0.0 or b > n:
        raise ValueError("b is not in range")

    if np.fabs((a + b)/2.0 - d) > 1e-8:
        raise ValueError("does not reconstruct d")
    if np.fabs((a - b)/np.sqrt(2.0 * (a + b)) - lam) > 1e-8:
        raise ValueError("does not reconstruct lambda")

    return a, b


def goOnce(idx, n, d, lam, trial, a, b, seedSample, seedOpt):
    cmd = "./instance --n {} --r {} --a {} --b {} --seed-gen {} --seed-opt {}".format(n, 100, a, b, seedSample, seedOpt)
    print cmd
    subprocess.call(cmd + " > {}.json".format(idx), shell=True)


def main():
    lambdas = [1.0, 1.005, 1.01, 1.015, 1.02, 1.025, 1.03]
    ds = [5.0]
    ns = [16000]
    trials = 1

    seenSampleSeeds, seenOptSeeds = set(), set()

    def grid():
        for idx, (n, d, lam, trial) in enumerate(it.product(ns, ds, lambdas, range(trials))):
            a, b = parameters(n, d, lam)

            seedSample = np.random.randint(low=0, high=0xFFFFFFFF)
            while seedSample in seenSampleSeeds:
                seedSample = np.random.randint(low=0, high=0xFFFFFFFF)
            seenSampleSeeds.add(seedSample)

            seedOpt = np.random.randint(low=0, high=0xFFFFFFFF)
            while seedOpt in seenOptSeeds:
                seedOpt = np.random.randint(low=0, high=0xFFFFFFFF)
            seenOptSeeds.add(seedOpt)

            yield idx, n, d, lam, trial, a, b, seedSample, seedOpt

    Parallel(n_jobs=mp.cpu_count() / 2, backend='threading')(
        delayed(goOnce)(*args) for args in grid())


if __name__ == '__main__':
    main()
