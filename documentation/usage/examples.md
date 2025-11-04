# PROCESS Examples

Examples of how to use PROCESS are provided in the form of Jupyter notebooks. These provide a 
convenient way of explaining PROCESS usage alongside runnable code cells, as well as displaying 
some types of output.

Renders of these notebooks can be found under the same sub-heading as this page. We would encourage users to run the notebooks themselves once they have read through the renders as this will help build experience with running and visualising PROCESS.


### Notebooks in VS Code

The recommended way to run notebooks is in VS Code; this has the additional advantage of being able to debug notebooks. Simply open the `.ipynb` file in VS Code to view and run it. You may be required to select a notebook kernel on first run; be sure to select the virtual environment where PROCESS is installed (e.g. `~/PROCESS/.venv`).

### Notebooks via the Jupyter server

Another way of running Jupyter notebooks is via a web browser and the Jupyter server. Start by navigating to your PROCESS directory and activate your virtual environment (e.g. `source .venv/bin/activate`). 

You should then ensure that the `jupyter` and `notebook` packages are installed by running:

```bash
pip install jupyter
```

Next, navigate to the `examples` directory within PROCESS

```bash
cd examples/
```

and then run:

```bash
python -m notebook
```

A web browser will open and the notebook can be run from there. If you're using WSL, you'll have to `ctrl + click` the link.

### Notebooks via Binder
Another way of running the PROCESS example notebooks is to use Binder. This uses a JupyterHub server to host the contents of
PROCESS, allowing the examples to be run via a web browser and without installation on your computer. You can click 
[here](https://mybinder.org/v2/gh/ukaea/PROCESS/HEAD) to try this out.
The Binder may take some time to load, but once loaded you will see options down the left hand side - double click on `examples` and then `examples.ipynb`, then you can run this notebook in your web browser.

## Maintaining notebooks

Notebooks are located in the `examples` directory and are tested using `testbook` to ensure they keep working.


## Issues running notebooks

If you encounter a `PDFInfoNotInstalledError` when running a notebook, ensure poppler utilities are properly installed:

```bash
apt-get update
apt-get install poppler-utils
```