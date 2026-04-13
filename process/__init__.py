import logging
from importlib.metadata import version

# setup the package "root" logger, every other logger in PROCESS
# will inherit from this (implicitly, based on the namespace)
root_logger = logging.getLogger(__name__)  # noqa: RUF067
root_logger.setLevel(logging.DEBUG)  # noqa: RUF067

__version__ = version("process")
