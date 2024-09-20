# Process Jupyter notebook examples

Examples of how to use Process are provided in the form of Jupyter notebooks. These provide a convenient way of explaining Process usage alongside runnable code cells, as well as displaying some types of output.

## Running notebooks

### Notebooks in VS Code

The recommended way to run notebooks is in VS Code; this has the additional advantage of being able to debug notebooks. Simply open the `.ipynb` file in VS Code to view and run it. You may be required to select a notebook kernel on first run; be sure to select the virtual environment where PROCESS is installed (e.g. `~/process/env`).

### Notebooks via the Jupyter server

Another way of running Jupyter notebooks is via a web browser and the Jupyter server. Start by navigating to your PROCESS directory and activate your virtual environment (e.g. `source env/bin/activate`). Next, navigate to the `examples` directory within PROCESS

```bash
cd examples/
```

and then run:

```bash
python -m notebook
```

A web browser will open and the notebook can be run from there. If you're using WSL, you'll have to `ctrl + click` the link.

### Installation

`jupyter` (which includes the `notebook` package) should already be installed in the `process` virtual environment when the `cmake` installation is performed (it is included in the `requirements.txt`). If not, install with `pip install jupyter` (for all Jupyter packages) or `pip install notebook` (for just the `notebook` package).

## Maintaining notebooks
Notebooks are located in the `examples` directory and are tested using `testbook` to ensure they keep working.
