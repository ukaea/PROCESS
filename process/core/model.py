"""Module containing routines to set up the data structure and models"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from process.core.data_structure.base import DataStructure


class Model(abc.ABC):
    """Set up Model base class"""

    data: DataStructure

    @abc.abstractmethod
    def run(self) -> None:
        """Run the model.

        The run method is responsible for 'running' the model, ensuring it updates
        the data structure with variables that subsequent models will require.
        """

    @abc.abstractmethod
    def output(self) -> None:
        """Output model data.

        This method will always be called after run method and should output the model
        data to the output files.
        """
