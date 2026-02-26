# PROCESS Examples

Examples of how to use PROCESS are provided in the form of Jupyter notebooks. These provide a 
convenient way of explaining PROCESS usage alongside runnable code cells, as well as displaying 
some types of output.

Renders of these notebooks can be found under the same sub-heading as this page. We would encourage users to run the notebooks themselves once they have read through the renders as this will help build experience with running and visualising PROCESS.

Please note that when trying to run the examples, it is important to have PROCESS installed editably to ensure the data files can be located. Without an editable install, most examples will not run.


### Notebooks in VS Code

The recommended way to run notebooks is in VS Code; this has the additional advantage of being able to debug notebooks. Simply open the `.ex.py` file in VS Code and click run in interactive mode to view and run it. You may be required to select a notebook kernel on first run; be sure to select the virtual environment where PROCESS is installed (e.g. `~/PROCESS/.venv`).

### Notebooks via the Jupyter server

Another way of running Jupyter notebooks is via a web browser and the Jupyter server. Start by navigating to your PROCESS directory and activate your virtual environment (e.g. `source .venv/bin/activate`). 

You should then ensure that examples required dependencies are installed using the examples extra:

```bash
pip install -e .'[examples]'
```

Next, navigate to the `examples` directory within PROCESS

```bash
cd examples/
```

and then run:

```bash
jupyter lab
```

A web browser will open and the notebook can be run from there. If you're using WSL, you'll have to `ctrl + click` the link.

### Notebooks via Binder
Another way of running the PROCESS example notebooks is to use Binder. This uses a JupyterHub server to host the contents of
PROCESS, allowing the examples to be run via a web browser and without installation on your computer. You can click 
[here](https://mybinder.org/v2/gh/ukaea/PROCESS/HEAD?urlpath=%2Fdoc%2Ftree%2Fexamples%2F) to try this out.
The Binder may take some time to load, but once loaded you will be in the `examples` folder and can select example notebooks to run in your web browser. Simply click the notebook you wish to run (a `.ex.py` file) and select `Open with > Notebook`.
