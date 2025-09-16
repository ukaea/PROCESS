# Python Optimisation and Numba
PROCESS was originally written in Fortran, a very fast, compiled code. Python is dynamically 
interpreted and is comparatively slow. For most sections of PROCESS, this speed difference is 
negligible. However, certain sections of the code that perform numerically intensive calculations 
can cause noticeable increases in the runtime when these evaluations are happening thousands of 
time throughout the multiple iterations PROCESS does attempting to find a solution.

## Numba

Our policy for dealing with slow code revolves around [Numba](https://numba.readthedocs.io/en/stable/index.html) 
-- a project that just-in-time (jit) compiles Python code into machine code, just as Fortran was compiled.

This allows us to achieve near-Fortran speed while still writing Python code. The caveat is that our 
code must be compilable. Essentially, a function must be called with Numba-compliant types **only**. 
Further, Numba does not support all Python code, see [here](https://numba.readthedocs.io/en/stable/reference/pysupported.html) 
for more information.

!!! Info "First-run compilation"
    Because Numba is JIT compiled, it does not compile when PROCESS is installed, it is instead 
    compiled when that function is first called. This means that the run after a fresh installation 
    of PROCESS will be rather slow. Subsequent runs of PROCESS will use the cached compilation and 
    so will be much faster. This will also happen if you make a change to a Numba'd function while 
    in an editable install.

### Numba examples

```python
from numba import njit
import numpy as np

@njit
def my_function(a, b):
    c = a[0] + a[1]

    return c + b
```

Here, it is obvious that `a` is a `list` or Numpy `ndarray` while `b` is a `float` (could also be 
an `int`, both can be used interoperably). This means that we can call this function as follows:

```python
my_function(np.array([1, 2]), 3.0)
```

But, also **cannot** do the following:

```python
from numba import njit
from process.data_structure import superconducting_tf_coil_variables as sctfv

@njit
def my_other_function(n):
    return n + sctfv.n_tf_coils
```

because Numba does not know what `sctfv` is.

!!! Info "Numba benefits"
    The above examples are simple by design, however this also means they are poor candidates for 
    Numba'ing as there would no noticeable speed increase from doing so.

### Using Numba

Because of these limitations we advise Numba is used on pure mathematical functions that take 
numerical inputs, and return numerical outputs. This should not be done with object oriented code.

Numba has great support for mathematical operations and Numpy functions (including matrix 
operations). Generally, Numba will only be used by RSE's after a section of the code is profiled 
and found to be a source of slowness. Please consult us before using it as the error messages can 
be rather cryptic and difficult to debug.

### Debugging with Numba

Numba does not support debugging with `pdb`. To debug a Numba'd function, comment out the `@njit` 
decorator before debugging (making sure to uncomment it before committing any code).

## Profiling slow code

The following series of commands will generate a call trace graph of a PROCESS run the `large-tokamak` 
regression test example, along with a colour indicating its relative runtime. This can be used to 
quickly identify which parts of PROCESS are taking the most time, and are likely good candidates 
for optimisation. 

First ensure that `gprof2dot` is installed by issuing the command `pip install gprof2dot`.

!!! Warning "`dot`"
    `dot` is also required to run this command, but requires separate installation (cannot be done 
    via `pip`). For Linux users, the command to install `dot` is: `sudo apt-get install graphviz`.

Execute the following commands from the PROCESS root directory:

```bash
> python3 -m cProfile -o large-tok.pstats process/main.py -i tests/regression/scenarios/large-tokamak/IN.DAT

> gprof2dot -f pstats large-tok.pstats | dot -Tpng -o large_tokamak_profile.png
```

Dark blue nodes indicate this function had a low runtime. Light blue, green, or red nodes indicate 
a function that consumed an increasing amount of the total runtime. Each node lists two percentages:

1. The total calltime: the total percentage of the code runtime that happened during a call to 
   this function (including subsequent function calls it makes).
2. The self calltime: identified by this percentage being enclosed in brackets. It identifies how 
   long was spent in this function alone (so does not include the runtime functions which it calls).
   