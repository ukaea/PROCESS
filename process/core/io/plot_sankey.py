"""
Code to display the power flow of a PROCESS run in a Sankey diagram

Input file:
MFILE.DAT
"""

import argparse
import pathlib

from pylab import savefig, show

from process.core.io.sankey_funcs import plot_sankey


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Program to plot the power flow in PROCESS using a Sankey diagram."
    )

    parser.add_argument("-e", "--end", default="pdf", help="file format, default = pdf")

    parser.add_argument(
        "-m", "--mfile", default="MFILE.DAT", help="mfile name, default = MFILE.DAT"
    )

    args = parser.parse_args(args)

    plot_sankey(args.mfile)

    # Get directory of mfile
    mfile_path = pathlib.Path(args.mfile).resolve()
    mfile_dir = mfile_path.parent
    output_path = mfile_dir / f"SankeyPowerFlow.{args.end}"
    savefig(str(output_path))

    show()


if __name__ == "__main__":
    main()
