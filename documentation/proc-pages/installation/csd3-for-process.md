# Using the CSD3 service to run PROCESS

[The Cambridge Service for Data Driven Discovery (CSD3)](https://www.hpc.cam.ac.uk/high-performance-computing) is one of the EPSRC Tier2 National HPC Services hosted by Research Computing Services at the University of Cambridge.

This facilty can be used to run batch jobs with PROCESS. This section shall explain how to register for use of the cluster, install PROCESS and submit batch jobs. Find the full documentation for CSD3 [here](https://docs.hpc.cam.ac.uk/hpc/index.html).

## Register to use CSD3

- Register for a [SAFE Account](https://epcced.github.io/safe-docs/safe-for-users/)
- Once your account is created request to the relevant project and wait for approval.
- When you have recieved an approval confirmation via email your account will be active.
- Now you can ssh into CSD3. Note that you must register for two-factor authentication upon your first login request.

    ```bash
    ssh username@login-cpu.hpc.cam.ac.uk
    ```

## Install PROCESS on CSD3

PROCESS can be installed using miniconda, a lightweight version of Anaconda which is a distribution of Python and R for scientific computing that simplifies package management.

### Install conda and set up the environment

First, install miniconda on CSD3 using the latest linux version from Anaconda.com

``` bash
wget <https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh>

chmod +x Miniconda3-latest-Linux-x86_64.sh

./Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc
```

You should now a `(base)` command-prompt prefix, if not, restart your shell.

Then create a conda environement and activate it.

``` bash
conda create --name my_process

conda activate my_process
```

It is worth checking that a compatible version of python has been installed (Python3.8 or 3.10). You can check this by entering `which python` into the terminal. You can install a compatible version of python with the command `conda install python=3.10`.

Next, install gfortran.

```bash
conda install -c conda-forge gfortran
```

### Build PROCESS

Now you should follow the installations for PROCESS. Set up ssh keys, clone the repository from gitlab and then install requirements from within the process folder:

```bash
pip install -r requirements.txt
```

Build then test:

```bash
cmake -S . -B build

cmake --build build

pytest tests/regression
```

## Submitting jobs

CSD3 uses a SLURM queuing system to manage job submission. Follow the instructions in the user guide for further details [here](https://docs.hpc.cam.ac.uk/hpc/user-guide/batch.html).

Modfify one of the SLURM job submission scripts, which are found in `/usr/local/Cluster-Docs/SLURM/` with your account information. You can add some commands for the job script to execute. For example, the job script can contain instructions to activate your environment:

``` bash
conda activate my_process
```

and then execute a command, such as a single run of PROCESS:

```bash
#!/bin/sh
process -i ~/process/tests/scenarios/regression/large-tokamak/IN.DAT
```

Note that you must use the correct project code for you work.

Use the `sbatch` command to submit your script. You will recieve an email notifying you that the job has begun, and one to notify you when it has ended successfully or failed.
