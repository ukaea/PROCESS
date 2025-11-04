"""Run the tracked files and move into tracking directory."""

import argparse
import shutil
import subprocess
from pathlib import Path

from tracking_data import ProcessTracker


def run_and_move_tracked_files(arguments):
    """
    Run PROCESS on the tracked files and move the MFiles to the tracking directory.
    """
    path_to_tracking_files = Path(arguments.input_locs)

    input_files = path_to_tracking_files.glob("*.IN.DAT")
    tracking_dir = Path(__file__).parent
    for file_path in input_files:
        try:
            subprocess.run(
                f"process -i {file_path} {arguments.options}", shell=True, check=True
            )
        except subprocess.CalledProcessError:
            continue
        created_mfile = file_path.with_name(
            file_path.name.replace(".IN.DAT", ".MFILE.DAT")
        )
        moved_mfile = created_mfile.with_name(
            created_mfile.name.replace(".MFILE.DAT", "_MFILE.DAT")
        )
        created_mfile.rename(tracking_dir / moved_mfile.name)


def tracking(arguments):
    """
    Call the ProcessTracker, move the MFiles to the database and rename.
    """
    path_to_tracking_mfiles = Path(__file__).parent
    mfiles = path_to_tracking_mfiles.glob("*_MFILE.DAT")
    for mfile_path in mfiles:
        ProcessTracker(
            mfile=mfile_path,
            database=arguments.db,
            message=arguments.commit,
            hashid=arguments.hash,
            tracking_variables_file=arguments.tracking_variables_file,
        )

        copied_mfile = shutil.copy(mfile_path, Path(arguments.db) / mfile_path.name)
        moved_mfile_name = copied_mfile.with_name(
            copied_mfile.name.replace("_MFILE.DAT", f"_MFILE_{arguments.hash}.DAT")
        )
        copied_mfile.rename(moved_mfile_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="Choose the usage mode", dest="command")

    # Run command
    subparser_run = subparsers.add_parser("run")
    subparser_run.add_argument("input_locs", type=Path)
    subparser_run.add_argument("--options", type=str, default="")

    # Tracking command
    subparser_trk = subparsers.add_parser("track")
    subparser_trk.add_argument("db", type=str)
    subparser_trk.add_argument("commit", type=str)
    subparser_trk.add_argument("hash", type=str)
    subparser_trk.add_argument(
        "--tracking-variables-file",
        type=Path,
        default=None,
        help="A JSON file containing a list of variables to track."
        "See the description of DEFAULT_TRACKING_VARIABLES for details on formatting the strings in the list.",
    )

    arguments = parser.parse_args()

    if arguments.command == "run":
        run_and_move_tracked_files(arguments)
    else:
        tracking(arguments)
