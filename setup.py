from setuptools import setup, find_packages
import site
import os
import platform

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
            "data/fluids/*",
            "data/h_data/*",
            "data/lz_non_corona/*",
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
        "numpy>=1.23",
        "scipy>=1.10",
        "cvxpy!=1.3.0,!=1.3.1",
        "pandas",
        "tables",
        "SALib",
        "numba>=0.57",
        "PyVMCON>=2.1.0,<3.0.0",
        "CoolProp>=6.4",
        "seaborn>=0.12.2",
    ],
    "extras_require": {"test": ["pytest"]},
    "entry_points": {
        "console_scripts": [
            "process_script=process.process_script_advanced:main",
            "process=process.main:main",
        ]
    },
    "extra_link_args": EXTRA_ARGS,
}

if __name__ == "__main__":
    setup(**setup_kwargs)
