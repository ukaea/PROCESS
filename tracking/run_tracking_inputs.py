"""Run the tracked files and move into tracking directory."""

import argparse
from pathlib import Path
import subprocess

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
            subprocess.run(f"process -i {file_path}", shell=True, check=True)
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
    db = arguments.db
    commit = arguments.commit
    hash = arguments.hash
    path_to_tracking_mfiles = Path(__file__).parent
    mfiles = path_to_tracking_mfiles.glob("*_MFILE.DAT")
    for mfile_path in mfiles:
        ProcessTracker(mfile=mfile_path, database=db, message=commit, hashid=hash)
        moved_mfile_name = mfile_path.with_name(
            mfile_path.name.replace("_MFILE.DAT", f"_MFILE_{hash}.DAT")
        )
        copied_mfile_destination = Path(db) / moved_mfile_name.name
        copied_mfile_destination.write_bytes(mfile_path.read_bytes())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="Choose the usage mode", dest="command")
    subparser_run = subparsers.add_parser("run")
    subparser_trk = subparsers.add_parser("track")
    subparser_run.add_argument("input_locs", type=str)
    subparser_trk.add_argument("db", type=str)
    subparser_trk.add_argument("commit", type=str)
    subparser_trk.add_argument("hash", type=str)

    arguments = parser.parse_args()

    if arguments.command == "run":
        run_and_move_tracked_files(arguments)
    else:
        tracking(arguments)
