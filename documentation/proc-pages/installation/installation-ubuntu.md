# Installation on Ubuntu/Windows
PROCESS is developed using Ubuntu 22 and tested using Ubuntu 20. We cannot guarantee any other operating system will be able to compile PROCESS or reproduce results. We do unofficially support MacOS however PROCESS is **not** currently tested on this OS by the CI system. 

!!! Info "Windows User"
    Windows users should run PROCESS using WSL.

    To install Windows Subsystem for Linux (WSL) follow the 'Manual Installation Steps' 
    [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Choose WSL 2 and Ubuntu 20 (if 
    installing from the Microsoft store then Ubuntu 20 is installed by default). 

    The following command is used to install WSL (**you need admin privileges to perform this command**):

    ```bash
    wsl --install
    ```

    If the above procedure fails to work, there is a [Microsoft help page](https://learn.microsoft.com/en-us/windows/wsl/install).


GFortran version 9 or above is needed for successful installation and execution of PROCESS. Versions 
below GFortran-9 will be rejected by CMake by default since, while PROCESS might compile 
successfully with lower GFortran versions, other aspects of PROCESS (tests, coverage, etc.) will fail.

Open the terminal and install `cmake`, `gfortran`, `pip`, etc.:

```bash
sudo apt update
sudo apt install -y cmake gfortran python3-pip lcov poppler-utils python3-venv
```

Next, the code will need to be downloaded so you can work with it. The PROCESS code is stored in a 
GitHub repository and as such needs to be 'cloned' - i.e bought to your VSCode window from GitHub. 

Clone the PROCESS repository from GitHub, and navigate into the resulting directory:

```bash
git clone https://github.com/ukaea/PROCESS
cd PROCESS
```

Create and activate a virtual environment:

```bash
python3 -m venv env
source env/bin/activate
```

!!! Note "Virtual environment"
    While not mandatory, it is highly recommended users create a Python virtual environment in 
    order to use the PROCESS Python package, as this ensures that installations of required package 
    versions don't affect other packages that require different versions in your environment.

    In the installation instructions above, a virtual environment is mentioned for using to work 
    on PROCESS. It might not be so clear as to why a virtual environment is used however.

    Imagine a scenario where you are working on multiple projects and each project uses different 
    versions of different dependencies. This can be confusing and cause issues down the line. A 
    virtual environment is an aid that will allow maintenance of dependencies within each project. 
    Essentially, having a virtual environment for each project you will work on will store the 
    dependencies of those specific projects within itself, preventing crossover. This is the reason 
    we recommend installing one so your other projects don't influence the dependencies of PROCESS 
    and visa versa.

Now compile the Fortran and create the Python interface. This is done using `cmake` to configure the 
build and then `make` to build it. Finally start the build process:

```bash
cmake -S . -B build
cmake --build build
```

If you plan on developing code for PROCESS, please see the `pre-commit` documentation for 
installing this tool required by developers: [development/pre-commit](http://process.gitpages.ccfe.ac.uk/process/development/pre-commit/)

The build step may take some time when run for the first time (~3 mins) as the Fortran code is 
compiled and then wrapped using `f2py` to create the Python libraries. Once this is completed 
the PROCESS Python package is then automatically installed using `pip` and should be ready to use 
on Linux. If the installation was successful the command `process` should be available on the command line.

To rebuild, for example after making a change to the Fortran source, run `cmake --build build` again. 
Python-only changes are reflected immediately, as the `cmake` build script performs a development (editable) installation by default.

!!! warning "Users of multiple branches"
    For users of PROCESS that run on multiple branches, it is recommended that each time you switch branches you **fully rebuild** PROCESS:
    ```bash
    rm -r build
    cmake -S . -B build
    cmake --build build
    ```

The PROCESS test suite provides through tests that can be used to confirm a successful installation; 
the tests can then be used to verify changes you make have not affected the wider codebase.

Firstly, ensure you are in the PROCESS root directory.

```BASH
cd PROCESS
```

The test suite uses PyTest and can be fully run using:

```BASH
pytest
```

which runs unit, integration, and regression tests.

A more in-depth discussion of testing can be found [here](https://ukaea.github.io/PROCESS/development/testing/).

If everything passes, this indicates a successful installation. If anything fails, this indicates 
that your environment produces different results to what is expected.

!!! Question "Installation troubleshooting"

    === "PYTHONPATH"

        If you have previously modified your `$PYTHONPATH` environment variable to
        include `process/utilities`, perhaps in your `~/.bashrc` file, then please remove this modification.
        Re-start your terminal for the changes to take effect, and check this is not on your `$PYTHONPATH` with:

        ```bash
        echo $PYTHONPATH
        ```

        This modification is not required to run PROCESS now, and it may result in Ford failing during the build process otherwise.
        ```

    === "FORD Install"

        FORD installed to .local/bin

        An error may be encountered here because `ford` is installed in `.local/bin`, which is not on 
        the `PATH` in some environments, so you will need to add `.local/bin` to the `PATH` if this error 
        occurs. You can do this using `nano`:  

        First use:

        ```bash
        nano ~/.profile
        ```

        Then use the arrow keys to navigate to the bottom of the `nano` editor and type:

        ```bash
        export PATH=$PATH:/home/yourusername/.local/bin
        ``` 

        where you use your own username in place of `yourusername` above. Then use `Ctrl-X`, then type `Y`, then press enter. Then either close and reopen the terminal, or type:

        ```bash
        source ~/.profile
        ```

        Finally, remove the old build folder and try again:

        ```bash
        rm build -rf
        cmake -S . -B build
        cmake --build build
        ```

    === "CMake version"

        CMake version > 3.13.0

        CMake needs to be at least version `3.13.0`. This is so that the command `cmake -S . -B build` executes 
        correctly. Running this command on an earlier CMake version results in:  

        ```bash
        CMake Error: The source directory "/home/process/build" does not exist.
        Specify --help for usage, or press the help button on the CMake GUI.
        ```

        subsequently making the `build` directory and running the command again results in:

        ```bash
        CMake Error: The source directory "/home/process/build" does not appear to contain CMakeLists.txt.
        Specify --help for usage, or press the help button on the CMake GUI.
        ```

    === "Python dev version"

        Errors at the end of the build along the lines of `Python.h: No such file or directory #include "Python.h"` 
        could occur because the developer version of Python cannot be found. This can be fixed by 
        installing the developer version:

        ```bash
        sudo apt install libpython<python-version>-dev
        ```

        E.g.

        ```bash
        sudo apt install libpython3.8-dev
        ```
