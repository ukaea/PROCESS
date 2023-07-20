"""

    Script to convert IN.DAT into new format

    James Morris
    CCFE
    19/11/15

"""

# Third-party libraries
import argparse

# Local libraries
import process.io.in_dat as indat


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="PROCESS IN.DAT converter. "
        "Convert IN.DAT into new "
        "format.\n"
        "For info contact "
        "james.morris2@ccfe.ac.uk"
    )

    parser.add_argument(
        "-f",
        metavar="f",
        type=str,
        default="IN.DAT",
        help='File to read as IN.DAT (default="IN.DAT")',
    )

    parser.add_argument(
        "-o",
        metavar="o",
        type=str,
        default="new_IN.DAT",
        help='File to read as IN.DAT (default="new_IN.DAT")',
    )

    args = parser.parse_args()

    i = indat.InDat(filename=args.f)
    i.write_in_dat(output_filename=args.o)
