#!/usr/bin/env python
"""
Code to run multi-dimensional parameter scans using PROCESS.

Author: S. Torrisi (storrisi@u.rochester.edu)

Date: August 2014

Input files:
ndscan.json (configuration file, per default in working directory)

Outputfiles:


Compatible with PROCESS version 382

"""

import process_io_lib.NCDFfromMFILE as NCD
import argparse

if __name__ == "__main__":
    ############################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to create a netcdf\
 file from a multi-dimensional parameter scan."
    )

    PARSER.add_argument(
        "-f",
        "--configfile",
        default="ndscan.json",
        help="configuration file, default = ndscan.json",
    )

    ARGS = PARSER.parse_args()

    ############################################################
    # main program

    NCONVERTER = NCD.NCDFconverter(ARGS.configfile)
    NCONVERTER.convert_mfilelibrary()
