

# Installation using Singularity/Apptainer container

Singularity is a container environment similar to Docker. This means a user can run PROCESS with 
all required dependencies installed. Singularity, however, is designed to work with user-level 
permissions and, as such, is supported by many shared resource administrators (unlike Docker, 
which poses a security risk).

Singularity can convert OCI compliant containers into the Singularity Image Format (SIF) to run 
the Docker container above. Download, and convert the Docker container by running: 

```bash
apptainer pull process.sif docker://ghcr.io/ukaea/process-ci:latest
```

Singularity will then ask for a username and password, your CCFE GitLab short username and 
password. Singularity will write the container into your current directory, it can then be moved or 
copied like any file. Running the following command will load a Singularity shell with the 
dependencies for PROCESS installed:

```bash
apptainer shell process.sif
``` 

Singularity will automatically mount your home (`$HOME`) directory into the container. Therefore, 
if PROCESS lives in `~/process` on your system, it will also live inside of `~/process` in the 
shell environment.

It should also be noted that while the Singularity container has a Python 3.8 by default, it will 
be impossible to pip install any packages without getting a `Read-only file system` error. This is 
because you are treated as a non-admin user within the container, and, as such, you cannot update 
the system Python. For this reason, it is recommended that you still use a virtual environment 
within the Singularity container (as described above). `pip install <package> --user` will work; 
however, it will cause conflicts with existing Python packages you have installed outside of your container.

<h3>Steps to use Singularity to install PROCESS on Freia</h3>

Log onto Freia and make sure the Singularity module is loaded

```bash
module load singularity/3.7.1
```

Pull the Singularity image and enter your git username and password when prompted:

```bash
apptainer pull process.sif docker://ghcr.io/ukaea/process-ci:latest
```

Make sure you have created [ssh keys](https://docs.gitlab.com/ee/ssh/) and open the Singularity shell:

```bash
apptainer shell process.sif
```

Singularity is an environment which allows you to use dependencies not available on Freia, like 
the correct versions of Python and cmake. With the shell open, the installation of PROCESS can proceed as usual:

```bash
git clone git@git.ccfe.ac.uk:process/process.git
cd process
python3 -m venv env --without-pip --system-site-packages
source env/bin/activate
export PATH=$PATH:~USERNAME/.local/bin/
cmake -S . -B build
cmake --build build
```

Now you can run commands within the shell like `process -i tests/regression/scenarios/large-tokamak/IN.DAT` 
to verify installation and create [batch jobs](https://ukaea.github.io/PROCESS/io/utilities/).
