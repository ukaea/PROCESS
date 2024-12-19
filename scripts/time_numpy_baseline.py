"""Running this script will print the average runtime of a computationally intensive baseline.
This can then be used to normalise other runtime calculations to mitigate bias introduced
by hardware, system load, etc. The average runtime [s] will be printed to stdout."""

import time

import numpy as np


def numpy_baseline_runtime():
    """Runs a computational intense numeric baseline 10 times and
    returns the average runtime in seconds."""

    np.random.seed(42)

    N = 3500

    A = np.random.rand(N, N)
    B = np.random.rand(N, N)

    starts = []
    ends = []

    for _ in range(10):
        starts.append(time.time())

        A @ B

        ends.append(time.time())

    starts = np.array(starts)
    ends = np.array(ends)

    return (ends - starts).mean()


if __name__ == "__main__":
    print(f"Average runtime of the numpy baseline: {numpy_baseline_runtime()}")
