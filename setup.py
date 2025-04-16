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
    "version": "3.1.0",
    "description": (
        "Power Reactor Optimisation Code for Environmental and Safety Studies"
    ),
    "url": "https://ccfe.ukaea.uk/resources/process/",
    "author": "UKAEA",
    "packages": find_packages(),
    "package_dir": {"process": "process"},
    "package_data": {
        "process": [
            "lib/lib*",
            "fortran*.so",
            "data/lz_non_corona_14_elements/*",
            "utilities/*",
        ],
        "process.io": ["python_fortran_dicts.json"],
        "process.data.impuritydata": ["*"],
        "process.uncertainties": ["*.json"],
    },
    "test_suite": "pytest",
    "python_requires": ">=3.10",
    "install_requires": [
        "numpy>=1.23,<2",
        "scipy>=1.10",
        "numba-scipy @ git+https://github.com/numba/numba-scipy@23c3b33440ea1fe0f84d05d269fb4a3df4b92787",
        "cvxpy!=1.3.0,!=1.3.1",
        "osqp<1.0",
        "pandas>=2.0",
        "tables",
        "SALib",
        "numba>=0.57",
        "PyVMCON>=2.2.2,<3.0.0",
        "CoolProp>=6.4",
        "matplotlib>=2.1.1",
        "seaborn>=0.12.2",
        "tabulate",
    ],
    "extras_require": {
        "test": ["pytest>=5.4.1", "requests>=2.30", "testbook>=0.4"],
        "examples": ["pillow>=5.1.0", "jupyter==1.0.0", "pdf2image==1.16.0"],
    },
    "entry_points": {
        "console_scripts": ["process=process.main:main"],
        "numba_extensions": ["init = numba_scipy:_init_extension"],
    },
    "extra_link_args": EXTRA_ARGS,
}

if __name__ == "__main__":
    setup(**setup_kwargs)
