import os
import platform
import site

from setuptools import find_packages, setup

MODULE_NAME = "process"
_install_loc = os.path.join(site.getsitepackages()[0], MODULE_NAME)
EXTRA_ARGS = []
if platform.system() == "Darwin":
    EXTRA_ARGS = ["-Wl,-rpath," + os.path.join(_install_loc, "lib")]

setup_kwargs = {
    "name": MODULE_NAME,
    "version": "3.2.1",
    "description": (
        "Power Reactor Optimisation Code for Environmental and Safety Studies"
    ),
    "url": "https://ccfe.ukaea.uk/resources/process/",
    "author": "UKAEA",
    "packages": find_packages(),
    "package_dir": {"process": "process"},
    "package_data": {
        "process": ["data/lz_non_corona_14_elements/*"],
        "process.data.impuritydata": ["*"],
        "process.io": ["*.png"],
    },
    "test_suite": "pytest",
    "python_requires": ">=3.10",
    "install_requires": [
        "numpy>=1.23",
        "scipy>=1.10",
        "cvxpy!=1.3.0,!=1.3.1",
        "osqp>=1.0",
        "pandas>=2.0",
        "numba>=0.57",
        "PyVMCON>=2.3.1,<3.0.0",
        "CoolProp>=6.4",
        "matplotlib>=2.1.1",
        "seaborn>=0.12.2",
        "tabulate",
    ],
    "extras_require": {
        "test": ["pytest>=5.4.1", "requests>=2.30", "testbook>=0.4"],
        "examples": ["pillow>=5.1.0", "jupyter==1.0.0", "pdf2image==1.16.0"],
        "plotly": ["plotly>=5.15.0,<6"],
    },
    "entry_points": {"console_scripts": ["process=process.main:main"]},
    "extra_link_args": EXTRA_ARGS,
}

if __name__ == "__main__":
    setup(**setup_kwargs)
