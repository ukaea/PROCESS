import logging

from process.fortran import (
    constants,
)

logger = logging.getLogger(__name__)


class ResistiveTFCoil:
    def __init__(self):
        self.outfile = constants.nout
