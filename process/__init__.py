import logging
from importlib.metadata import version

# setup the package "root" logger, every other logger in PROCESS
# will inherit from this (implicity, based on the namespace)
root_logger = logging.getLogger(__name__)
root_logger.setLevel(logging.DEBUG)

__version__ = version("process")
