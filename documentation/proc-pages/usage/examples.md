# PROCESS Examples

Examples of how to use PROCESS are provided in the form of Jupyter notebooks. These provide a 
convenient way of explaining PROCESS usage alongside runnable code cells, as well as displaying 
some types of output.

<h2>Notebooks in VS Code</h2>

The recommended way to run notebooks is in VS Code; this has the additional advantage of being able to debug notebooks. Simply open the `.ipynb` file in VS Code to view and run it. You may be required to select a notebook kernel on first run; be sure to select the virtual environment where PROCESS is installed (e.g. `~/process/env`).

<h2>Notebooks via the Jupyter server</h2>

Another way of running Jupyter notebooks is via a web browser and the Jupyter server. Start by navigating to your PROCESS directory and activate your virtual environment (e.g. `source env/bin/activate`). Next, navigate to the `examples` directory within PROCESS

```bash
cd examples/
```

and then run:

```bash
python -m notebook
```

A web browser will open and the notebook can be run from there. If you're using WSL, you'll have to `ctrl + click` the link.


<h3>Installation</h3>

`jupyter` (which includes the `notebook` package) should already be installed in the `process` 
virtual environment when the `cmake` installation is performed (it is included in the 
`requirements.txt`). 

If not, install with 

```bash
pip install jupyter
``` 

for all Jupyter packages or 

```bash
pip install notebook
``` 

for just the `notebook` package.


<h2>Maintaining notebooks and scripts</h2>

Notebooks are located in the `examples` directory and are tested using `testbook` to ensure they keep working.


<h2>Issues running notebooks</h2>

If you encounter a `PDFInfoNotInstalledError` when running a notebook, ensure poppler utilities are properly installed:

```bash
apt-get update
apt-get install poppler-utils
```