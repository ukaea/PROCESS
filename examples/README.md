# Process Jupyter notebook examples

Examples of how to use Process are provided in the form of Jupyter notebooks. These provide a convenient way of explaining Process usage alongside runnable code cells, as well as displaying some types of output.

Please make sure the extra example dependencies are installed

```bash
pip install -e .'[examples]'
```

It is important that the install is editable, otherwise the examples will be unable to find the required data to run.

## Running notebooks

### Notebooks in Binder
Navigate to [our Binder](https://mybinder.org/v2/gh/ukaea/PROCESS/HEAD?urlpath=%2Fdoc%2Ftree%2Fexamples%2F), this should open a Binder in the `examples` directory. Right click on the `.ex.py` notebook you wish to open and select `Open with > Notebook`.

### Notebooks in VS Code

The recommended way to run notebooks is in VS Code; this has the additional advantage of being able to debug notebooks. Simply open the `.ex.py` file in VS Code and click run in interactive mode to view and run it. You may be required to select a notebook kernel on first run; be sure to select the virtual environment where PROCESS is installed (e.g. `~/PROCESS/.venv`).

If this is not possible or an error occurs, please consult [our docs](https://ukaea.github.io/PROCESS/usage/examples/).

## Notebooks in Jupyter Notebook >v7.x (nb7)

As described [here](https://jupytext.readthedocs.io/en/latest/text-notebooks.html) previous versions of Jupyter Notebook (a.k.a nb classic) will open python scripts and markdown documents as notebooks by default. From v7 onwards, to open and run a text notebook (.py or .md extension) as a notebook in Jupyter:
 - right-click on the file,
 - select "Notebook" in the drop down menu.

Alternatively, you can change the default viewer to "Jupytext Notebook" if you would like particular file types to open as a notebook by default via a double click. This can be set via the command line, e.g., for .py and .md files, using:

```bash
jupytext-config set-default-viewer python markdown
```

## Creating new examples

We do not store notebook (`ipynb`) files on the repository but the pypercent representation of them. This can be generated from an existing notebook using the jupytext drop down within a running jupyter notebook window as described [here](https://jupytext.readthedocs.io/en/latest/paired-notebooks.html). Our examples should all have the extension `.ex.py` to distinguish them from normal python files.

Please adjust the metadata at the top of the new example `.ex.py` file to be:

```
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,ipynb
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---
```
