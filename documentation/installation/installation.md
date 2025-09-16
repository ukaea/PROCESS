# Installing PROCESS
PROCESS is a command-line Python program. As such, it should be supported by a wide range of operating systems and Python versions above 3.10. However, we only test PROCESS on the latest Ubuntu with Python 3.10 and, therefore, we cannot guarantee its stability on other platforms. 

!!! Info "Windows User"
    While it is probably possible to run PROCESS on Windows, we recommend that Windows users use WSL as PROCESS is not tested for nor built for Windows use.

    To install Windows Subsystem for Linux (WSL) follow the 'Manual Installation Steps' 
    [here](https://learn.microsoft.com/en-us/windows/wsl/install). Choose WSL 2 and Ubuntu 20 (if 
    installing from the Microsoft store then Ubuntu 20 is installed by default). 

    The following command is used to install WSL (**you need admin privileges to perform this command**):

    ```bash
    wsl --install
    ```


Start by downloading the PROCESS source code using `git`:

```bash
git clone https://github.com/ukaea/PROCESS
cd PROCESS
```

!!! Info "Virtual Environment"
    It is software best-practice, especially for scientific software, to install software and its dependencies in a virtual environment. This
    will stop PROCESS and its dependencies interfering with other software installed on your system.

    In Python, a virtual environment can be easily created by running:

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```

    Please note that every time you start a new terminal you will need to run:

    ```
    source .venv/bin/activate
    ```

    before trying to run PROCESS.
    

PROCESS and its dependencies can then be installed by simply running:

```bash
pip install .
```

PROCESS also uses `poppler-utils` to create PDF summary files, which can be installed by running:
```bash
sudo apt update
sudo apt install -y poppler-utils
```

This dependency is not mandatory as most of the functionality of PROCESS is available without it.

## Examples

If you are new to PROCESS, you may want to run some of the examples in the `examples/` directory. These will introduce you to the basic functionality of reactor design with PROCESS. The examples require additional dependencies that can be installed using:
```
pip install '.[examples]'
```

## Developers
Developers of PROCESS will want to install PROCESS in a slightly different way and with a few additional dependencies.

First, developers may require an *editable* install of PROCESS: where any changes made to PROCESS source code are automatically reflected when using PROCESS (there will be no need to re-install PROCESS when changes are made to the Python source code). This can be done by running:

```bash
pip install -e .
```

For developers looking to contribute code back to our GitHub repository via a pull request, installing `pre-commit` will ensure your code meets our quality standards. Please see our dedicated pre-commit documentation for details on how to install it [development/pre-commit](https://ukaea.github.io/PROCESS/development/pre-commit/). 


Finally, developers will want to check that their work has not broken anything in PROCESS. Our test suite runs large portions of the code and compares the result to the "correct" value, reporting any differences. The following optional dependencies of PROCESS will need to be installed to run the test suite:

```
pip install -e '.[examples, test]'
```

Documentation for running the test suite can be found at [development/testing](https://ukaea.github.io/PROCESS/development/testing/).
