"""
PROCESS I/O library configuration API.

Designed to be used by utility developers when creating configuration
mechanisms for their tool.
"""

import json
import logging

logger = logging.getLogger(__name__)


class ConfigurationParser(object):
    """Abstract parser class. Must be subclassed to be used.

    The parser should always put read-in data in the data property.
    """

    def __init__(self):
        self._data = None

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        """Validate the configuration is provided in a specific format."""
        self.data_validate(value)
        self._data = value

    @data.deleter
    def data(self):
        del self._data

    def data_validate(self, value):
        """Check that value corresponds to a specific data format."""
        logger.info("type of value: {}".format(type(value)))
        if not isinstance(value, dict) and value is not None:
            raise ValueError("Configuration data must be specified as a dictionary")


class JsonConfigParser(ConfigurationParser):
    """JSON configuration parser."""

    def __init__(self, filename):
        self.data = None
        try:
            with open(filename) as fh:
                config_file_data = json.load(fh)
                self.data = config_file_data
        except FileNotFoundError:
            logger.error("Cannot find configuration file {}".format(filename))
            pass


class Config(object):
    """Generic configuration for PROCESS tools. Read-only."""

    def __init__(self, config_file, parser=JsonConfigParser):
        self.config_file = config_file
        parser = JsonConfigParser(config_file)
        self.config_data = self._lowercase(parser.data)

    def _lowercase(self, objekt):
        if isinstance(objekt, list):
            return [self._lowercase(item) for item in objekt]
        elif isinstance(objekt, dict):
            return dict(
                (key.lower(), self._lowercase(value)) for key, value in objekt.items()
            )
        else:
            return objekt

    def _search_config_for(self, config, *keys):
        """Recursively search config (a dict) for keys."""
        try:
            if isinstance(keys[0], str):
                search_key = keys[0].lower()
            else:
                search_key = keys[0]
            value = config[search_key]
        except IndexError:
            raise
        except KeyError:
            raise
        except TypeError:
            raise

        if isinstance(config, dict) and len(keys) > 1:
            return self._search_config_for(value, *keys[1:])
        elif not isinstance(value, dict) and len(keys) > 1:
            raise KeyError("{} cannot be found in {}".format(search_key, value))
        else:
            return self._lowercase(value)

    def get(self, *config_keys, default=None):
        """
        Return configured value corresponding to config_keys if possible.

        For nested configs, sequential items in config_keys can be used.
        For example, if the configuration is:
        {"a": "b", "c": {1: "hello", 2: {"x": "alpha", "z": "beta"}}}
        you can access the value of "z" by calling get("c", 2, "z").

        """

        try:
            return self._search_config_for(self.config_data, *config_keys)
        except KeyError:
            if default:
                logger.info("Using default for {}".format(config_keys))
                return default
            else:
                logger.exception(
                    "Cannot find value or default for {} in configuration".format(
                        config_keys
                    )
                )
        except (IndexError, TypeError):
            raise
