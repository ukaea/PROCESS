"""Module for calculating plasma scrape off layers physics"""

import logging

from process.core import constants
from process.core.model import Model

logger = logging.getLogger(__name__)


class ScrapeOffLayer(Model):
    """Model for calculating plasma scrape off layers physics."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        """Run the model. This model cannot yet be 'run'."""

    def output(self) -> None:
        """Output plasma scrape off layer physics information."""