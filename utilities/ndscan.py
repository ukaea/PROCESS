#!/usr/bin/env python
"""
Code to run multi-dimensional parameter scans using PROCESS.

Author: S. Torrisi (storrisi@u.rochester.edu)

Date: August 2014

Input files:
ndscan.json (configuration file, per default in working directory)

"""


from process_io_lib.process_config import NdScanConfig
import argparse

if __name__ == "__main__":
    ############################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to run a multi-\
    dimensional parameter scan using PROCESS."
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

    NDSCANNER = NdScanConfig(ARGS.configfile)
    NDSCANNER.start_scan()
