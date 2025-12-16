"""
Code to read from a PROCESS MFILE and write values into a csv

Author: R Chapman (rhian.chapman@ukaea.uk)
Date: 26/05/2021
Updated: 06/06/2022

Input files:
mfile (default MFILE.DAT) as output from PROCESS
variable list (default mfile_to_csv_vars.json) as defined by user

Instructions:
- command line call: python mfile_to_csv.py -f </path/to/mfile.dat> -v </path/to/varfile.json>

Output file:
.csv will be saved to the directory of the input file

"""

# == import modules ==
# standard python modules
import argparse
import csv
import json
from pathlib import Path, PurePath

# PROCESS-specific modules
from process.io.mfile import MFile

# == define functions ==


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Read from a PROCESS MFILE and write values into a csv."
    )
    parser.add_argument(
        "-f",
        "--mfile",
        type=str,
        default="MFILE.DAT",
        help="Specify input mfile name, default = MFILE.DAT",
    )
    parser.add_argument(
        "-v",
        "--varfile",
        type=str,
        default="mfile_to_csv_vars.json",
        help="Specify file holding variable names, default = mfile_to_csv_vars.json",
    )

    return parser.parse_args(args)


def get_vars(vfile="mfile_to_csv_vars.json"):
    """Returns variable names from identified file.

    :param args: input JSON filename
    :type args: string
    :return: variable names
    :rtype: list
    """
    print("Fetching list of variables from", vfile)

    return json.loads(Path(vfile).read_text())["vars"]


def read_mfile(mfilename="MFILE.DAT", variables=None):
    """Returns specified variable values from identified file.

    :param args: input filename, variable names
    :type args: string, list
    :return: variable descriptions, names, and values
    :rtype: list of tuples
    """
    if variables is None:
        variables = []
    print("Reading from MFILE:", mfilename)

    m_file = MFile(mfilename)

    output_vars = []

    # for each variable named in the input varfile, get the description and data value
    for var_name in variables:
        if var_name not in m_file.data:
            print(f"Variable '{var_name}' not in MFILE. Skipping and moving on...")
        else:
            # In case of a file containing multiple scans, (scan = -1) uses the last scan value
            var_val = m_file.data[var_name].get_scan(-1)
            description = m_file.data[var_name].var_description
            var_data = (description, var_name, var_val)
            output_vars.append(var_data)

    return output_vars


def get_savenamepath(mfilename="MFILE.DAT"):
    """Returns path/filename.csv for file saving.

    :param args: input filename
    :type args: string
    :return: output filename
    :rtype: pathlib.PurePosixPath
    """

    # Either save it locally or output the csv file to the directory of the input file
    dirname = Path.cwd() if mfilename == "MFILE.DAT" else PurePath(mfilename).parent

    csv_filename = PurePath(mfilename).stem
    return PurePath(dirname, csv_filename + ".csv")


def write_to_csv(csv_outfile, output_data=None):
    """Write to csv file.

    :param args: input filename, variable data
    :type args: string, list of tuples
    """
    if output_data is None:
        output_data = []
    with open(csv_outfile, "w") as csv_file:
        print("Writing to csv file:", csv_outfile)
        writer = csv.writer(csv_file, delimiter=",")
        writer.writerow(["Description", "Varname", "Value"])

        for vardesc in output_data:
            writer.writerow(vardesc)


def main(args=None):
    """Extract certain variables from an MFILE.DAT and output to CSV.

    :param args: optional command-line args for testing, defaults to None
    :type args: list, optional
    """
    # read from command line inputs
    args = parse_args(args)

    # read list of required variables from input json file
    jvars = get_vars(args.varfile)

    # read required data from input mfile
    output_data = read_mfile(args.mfile, jvars)

    # identify save location
    output_file = get_savenamepath(args.mfile)

    # write to csv
    write_to_csv(output_file, output_data)

    # write final line to screen
    print("Complete.")


# == program ==

if __name__ == "__main__":
    main()

# == end ==
