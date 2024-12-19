"""Running this script will run large tokamak and print its runtime [s] to stdout."""

import subprocess
import time
from pathlib import Path


def large_tokamak_runtime():
    """Runs large tokamak and returns the runtime in seconds"""
    script_dir = Path(__file__).resolve().parent
    large_tokamak_path = (
        script_dir / "../tests/regression/scenarios/large-tokamak/IN.DAT"
    )

    start = time.time()

    subprocess.run(
        f"process -i {large_tokamak_path.resolve()!s}", shell=True, check=False
    )

    return time.time() - start


if __name__ == "__main__":
    print(f"Runtime of large tokamak: {large_tokamak_runtime()}")
