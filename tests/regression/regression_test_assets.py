import subprocess
import requests
import dataclasses
import re
import logging
from typing import Optional
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class TrackedMFile:
    hash: str
    experiment_name: str
    download_link: str


class RegressionTestAssetCollector:
    remote_repository_url = (
        "https://api.github.com/repos/timothy-nunn/process-tracking-data/contents/"
    )

    def __init__(self) -> None:
        self._hashes = self._git_commit_hashes()
        self._tracked_mfiles = self._get_tracked_mfiles()

    def get_reference_mfile(
        self, scenario_name: str, directory: Path, target_hash: Optional[str] = None
    ):
        """Finds the most recent reference MFile for `<scenario_name>.IN.DAT`
        and downloads it to the `directory` with the name `ref.<scenario_name>.MFILE.DAT`.

        Providing a `target_hash` will ONLY return a reference MFILE that exactly
        matches the requested commit hash.

        :returns: Path to the downloaded reference MFile, if no reference MFile can be found,
        `None` is returned.
        :rtype: Path
        """
        reference_mfile_location = directory / f"ref.{scenario_name}.MFILE.DAT"
        for mf in self._tracked_mfiles:
            if (mf.experiment_name == scenario_name and target_hash is None) or (
                mf.experiment_name == scenario_name and target_hash == mf.hash
            ):
                with open(reference_mfile_location, "w") as f:
                    f.write(requests.get(mf.download_link).content.decode())

                logger.info(f"Reference MFile found for commit {mf.hash}")
                return reference_mfile_location

        return None

    def _git_commit_hashes(self):
        return (
            subprocess.run(
                'git log --format="%H"',
                shell=True,
                capture_output=True,
                check=True,
            )
            .stdout.decode()
            .split("\n")
        )

    def _get_tracked_mfiles(self):
        repository_files = requests.get(self.remote_repository_url).json()

        return sorted(
            [
                mfile
                for f in repository_files
                if (mfile := self._get_tracked_mfile(f)) is not None
                and mfile.hash in self._hashes
            ],
            key=lambda m: self._hashes.index(m.hash),
        )

    @staticmethod
    def _get_tracked_mfile(json_data):
        """Converts JSON data of a file tracked on GitHub into a `TrackedMFile`, if appropriate"""
        rematch = re.match(r"([a-zA-Z0-9_.]+)_MFILE_([a-z0-9]+).DAT", json_data["name"])

        if rematch is None:
            return None
        return TrackedMFile(
            rematch.group(2), rematch.group(1), json_data["download_url"]
        )
