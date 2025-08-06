###############################################################################
#                                                                             #
#                      MFILE to Python Dictionary Conversion                  #
#                                                                             #
#   Converts data from a PROCESS MFILE to a Python Dictionary with the        #
#   option to then write the result to a file. The format of the output is    #
#   determined by the specified file type which can be JSON, TOML, YAML or    #
#   a Pickle file. If TOMLKit is available the output contains docstrings.    #
#                                                                             #
#   @author :   K. Zarebski <kristian.zarebski@ukaea.uk>                      #
#   @date   :   last modified 2021-02-22                                      #
#                                                                             #
###############################################################################
import logging
import os
import re
from collections import OrderedDict, abc
from typing import Any

MFILE_END = "# Copy of PROCESS Input Follows #"
VETO_LIST = [" # PROCESS found a feasible solution #"]
HEADER_MAPPING = {"Power Reactor Optimisation Code": "metadata"}


class MFILEParser(abc.MutableMapping):
    def __init__(self, input_mfile: str = ""):
        """Parse an MFILE and extract output values.

        Parameters
        ----------
        input_mfile : str, optional
            MFILE to search, by default ""
        """
        self._input_file = input_mfile
        self._mfile_data: OrderedDict = OrderedDict()
        self._logger = logging.getLogger(self.__class__.__name__)
        if self._input_file:
            self.parse(self._input_file)

    def __iter__(self):
        for group in self._mfile_data:
            yield from self._mfile_data[group]

    def __len__(self):
        return sum([len(self._mfile_data[g]) for g in self._mfile_data])

    def items(self):
        for group in self._mfile_data:
            for param, value in self._mfile_data[group].items():
                yield param, value["value"]

    def __getitem__(self, key):
        for group in self._mfile_data:
            if key in self._mfile_data[group]:
                return self._mfile_data[group][key]["value"]
        raise KeyError(f"No variable '{key}' found.")

    def get_info_dict(self):
        """Get complete information dictionary.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            retrieve the full information dictionary containing values and
            descriptions of extracted parameters
        """
        return self._mfile_data

    def __delitem__(self, key):
        for group in self._mfile_data:
            if key in self._mfile_data[group]:
                del self._mfile_data[group][key]
                return
        raise KeyError(f"No variable '{key}' found.")

    def __setitem__(self, key, value):
        for group in self._mfile_data:
            if key in self._mfile_data[group]:
                self._mfile_data[group][key]["value"] = value
                return
        raise KeyError(f"No variable '{key}' found.")

    def get_parameter_value(self, param_name: str) -> Any:
        for group in self._mfile_data:
            if param_name in self._mfile_data[group]:
                return self._mfile_data[group][param_name]["value"]
        raise KeyError(f"No variable '{param_name}' found.")

    def _line_string_search(self, lines: list[str], search_str: str) -> list[int]:
        """Search for substring in file lines.

        Parameters
        ----------
        lines : List[str]
            list of file lines to search
        search_str : str
            search term to look for

        Returns
        -------
        List[int]
            list of indexes for matching lines
        """
        return [i for i, line in enumerate(lines) if search_str in line]

    def _find_var_val_from_str(self, value_str: str) -> Any:
        """Convert a string variable to float, int etc.

        This function parsers values given within the MFILE removing other
        unneeded information such as the specific typing in PROCESS

        Parameters
        ----------
        value_str : str
            value as a string

        Returns
        -------
        Any
            the value in the appropriate form
        """
        for type_str in ["OP", "IP", "ITV"]:
            value_str = value_str.replace(type_str, "")
        try:
            return int(value_str)
        except ValueError:
            pass
        try:
            return float(value_str)
        except ValueError:
            return value_str

    def _get_values(self, lines: list[str]) -> dict[str, Any]:
        """Extracts value, description and variable name from MFILE lines.

        Parameters
        ----------
        lines : List[str]
            list of file lines to be parsed

        Returns
        -------
        Dict[str, Any]
            dictionary containing variable name, value and description
        """
        # Compile regex for converting underscores which are spaces into
        # a space character
        _space_re = r"(\_{5,})"
        _var_re = r"(\([a-z0-9\-\+\*\/\_\%\]\[]+\))"

        # Extract lines from the given line set that follow the variable
        # statement convention of 'desc_______(varname)________ value'
        _lines = [line for line in lines if re.findall(_var_re, line)]

        # Remove extra symbols such as quotation marks and split line into
        # the three required components using regex
        _lines = [
            [
                i.replace('"', "").replace("`", "").strip()
                for i in re.split(_space_re, line)
                if not (re.findall(_space_re, i)) and i.strip()
            ]
            for line in _lines
        ]

        # If there are not three components in a given line, try splitting
        # the components present by ' ' instead and append
        for i, line in enumerate(_lines):
            if len(line) != 3:
                _new_line = []
                for element in line:
                    if " " in element:
                        _new_line += element.split(" ")
                    else:
                        _new_line += element
                _lines[i] = _new_line[:3]

        # Use an ordered dictionary to match ordering in MFILE
        _vars_dict = OrderedDict()

        # Iterate through the resultant line sets and tidy them a little
        # finally creating a dictionary entry for each with the required
        # information
        for line in _lines:
            _var_key = line[1][1:-1]
            _var_key = _var_key.replace("%", ".")
            if not _var_key:
                continue
            _value = line[2]
            _desc = line[0].replace("_-_", "-").replace("_", " ")
            _desc = _desc.title().strip()
            _desc = _desc.replace('"', "")
            _desc = re.sub(r"\s{2,}", " ", _desc)
            if _var_key in _vars_dict:
                if not isinstance(_vars_dict[_var_key], list):
                    _vars_dict[_var_key]["value"] = [_vars_dict[_var_key]["value"]]
                _vars_dict[_var_key]["value"].append(
                    self._find_var_val_from_str(_value)
                )
            else:
                _vars_dict[_var_key] = {
                    "description": _desc,
                    "value": self._find_var_val_from_str(_value),
                }

        return _vars_dict

    def parse(self, mfile_addr: str) -> dict:
        """Parse an MFILE and extract output values.

        Parameters
        ----------
        mfile_addr : str
            address of MFILE to parse

        Returns
        -------
        Dict
            dictionary of output values

        Raises
        ------
        FileNotFoundError
            if the input file does not exist
        """
        if not os.path.exists(mfile_addr):
            raise FileNotFoundError(
                f"Could not open MFILE '{mfile_addr}', file does not exist."
            )

        self._logger.info("Parsing MFILE: %s", mfile_addr)

        with open(mfile_addr) as f:
            _lines = f.readlines()

        _end_of_output = self._line_string_search(_lines, MFILE_END)[0]

        self._logger.info("Extracting file headers")
        _header_indexes = [
            i for i, line in enumerate(_lines) if line.strip() and i < _end_of_output
        ]

        _header_indexes = [
            i
            for i in _header_indexes
            if _lines[i].strip()[0] == "#"
            and not any(k in _lines[i] for k in VETO_LIST)
        ]

        # Gets rid of multi-headers, taking the last one
        _header_indexes = [i for i in _header_indexes if i + 1 not in _header_indexes]

        self._logger.info("Retrieving output variable values")
        # Iterate through the file headers processing the "block" between them
        # extracting variable values. Where duplicate headers are found assume
        # that a parameter sweep is occuring and append values in lists
        for i in range(len(_header_indexes) - 1):
            _key = _lines[_header_indexes[i]].replace("#", "").strip()

            _new_vals = self._get_values(
                _lines[_header_indexes[i] + 1 : _header_indexes[i + 1]]
            )

            # The iscan variable is always present at start of sweep
            # no matter what the first header in an iteration is
            # need to move it into metadata
            _first_key = _lines[_header_indexes[0]]
            _check_iscan = self._mfile_data and "iscan" in _new_vals
            _check_iscan = _check_iscan and _key != _first_key
            if _check_iscan:
                _first_key = _first_key.replace("#", "").strip()
                _iscan_var = self._mfile_data[_first_key]["iscan"]["value"]
                if not isinstance(_iscan_var, list):
                    self._mfile_data[_first_key]["iscan"]["value"] = [_iscan_var]
                self._mfile_data[_first_key]["iscan"]["value"].append(
                    _new_vals["iscan"]["value"]
                )
                del _new_vals["iscan"]

            # Add header to dictionary of not present
            if _key not in self._mfile_data:
                self._mfile_data[_key] = _new_vals

            # If header already present, iterate through member parameters
            # appending the new values to each
            else:
                for param, var_dict in self._mfile_data[_key].items():
                    if param not in _new_vals:
                        self._logger.warning(
                            f"Expected parameter '{param}' in sweep, "
                            "but could not find entry"
                            " for this iteration"
                        )
                        continue
                    _value = _new_vals[param]["value"]
                    if not isinstance(var_dict["value"], list):
                        self._mfile_data[_key][param]["value"] = [
                            self._mfile_data[_key][param]["value"],
                            _value,
                        ]
                    else:
                        # Need to check if the find variables function
                        # returned a single value for the parameter or multiple
                        # and handle the cases
                        if not isinstance(_new_vals[param]["value"], list):
                            self._mfile_data[_key][param]["value"].append(_value)
                        else:
                            self._mfile_data[_key][param]["value"] += _value

        self._logger.info("Creating output dictionaries")
        # Remove any cases where there are no parameters under a given header
        self._mfile_data = {k: v for k, v in self._mfile_data.items() if v}

        # Use underscore keys and tidy them to be more computationally friendly
        def _key_update(key):
            _key = key.lower()
            _key = _key.replace(" ", "_")
            for sym in [":", "(", ")", "/"]:
                _key = _key.replace(sym, "")
            return _key.replace("__", "_")

        # Apply header mappings and tidy headers
        self._mfile_data = {
            _key_update(k) if k not in HEADER_MAPPING else HEADER_MAPPING[k]: v
            for k, v in self._mfile_data.items()
        }

        if not self._mfile_data or len(self._mfile_data) == 0:
            raise AssertionError("Failed to extract data from given MFILE")

        # Only run iscan check if iscan exists
        try:
            _first_key = next(iter(self._mfile_data.keys()))
            _second_key = list(self._mfile_data.keys())[1]
            _second_key_fp = list(self._mfile_data[_second_key])[8]
            _iscan_arr = self._mfile_data[_first_key]["iscan"]["value"]
            _test_param = self._mfile_data[_second_key][_second_key_fp]["value"]
            if len(_test_param) != _iscan_arr[-1]:
                print(_test_param)
                raise AssertionError(
                    "Failed to retrieve all parameter sweep values, "
                    f"expected {_iscan_arr[-1]} values for '{_second_key}:{_second_key_fp}' and got {len(_test_param)}"
                )
        except KeyError:
            pass

        self._logger.info("Extraction completed successfully")
        return self._mfile_data

    def write(self, output_filename: str) -> None:
        """Write output to file.

        Parameters
        ----------
        output_filename : str
            path of output file, file type is determined from the type and can
            be '.toml', '.yml', '.pckl', '.json'
        """
        self._logger.info("Writing to output file '%s'", output_filename)

        _suffix = os.path.splitext(output_filename)[1].lower()

        if _suffix == ".toml":
            self._logger.info("Output will be TOML file.")
            try:
                import tomlkit
            except ImportError:
                # If file suffix is TOML but TOMLKit is not installed
                import toml

                print(
                    "WARNING: Python module 'tomlkit' not found, "
                    "file comments will not be written to created TOML file"
                )
                with open(output_filename, "w") as file:
                    toml.dump(self._mfile_data, file)
                exit(0)

            # If TOMLKit is present, write TOML file as normal but add in
            # descriptions as docstrings instead and format
            _doc = tomlkit.document()
            _doc.add(tomlkit.comment("PROCESS Run Output"))
            for group_name, data in self._mfile_data.items():
                _new_dict = {}
                for var_name, var_data in data.items():
                    _new_dict[var_name] = var_data["value"]
                _header = group_name.replace("_", " ").title()
                _ls = int((75 - len(_header)) / 2)
                _rs = 75 - len(_header) - _ls
                _header = _ls * "-" + " " + _header + " " + _rs * "-"
                _doc.add(tomlkit.comment(_header))
                _doc.add(group_name, _new_dict)
                _doc.add(tomlkit.nl())
                _doc.add(tomlkit.nl())

            for group_name, data in self._mfile_data.items():
                for var_name in data:
                    _doc[group_name][var_name].comment(
                        self._mfile_data[group_name][var_name]["description"]
                    )

            with open(output_filename, "w") as f:
                f.write(tomlkit.dumps(_doc))
        elif _suffix == ".json":
            # If file suffix is JSON
            self._logger.info("Output will be JSON file.")
            import json

            with open(output_filename, "w") as file:
                json.dump(self._mfile_data, file)
        elif _suffix in [".yml", ".yaml"]:
            self._logger.info("Output will be YAML file.")
            # If file suffix is YAML
            import yaml

            with open(output_filename, "w") as file:
                yaml.dump(self._mfile_data, file)
        elif _suffix == ".pckl":
            self._logger.info("Output will be Pickle file.")
            # If file suffix is Pickle
            import pickle

            with open(output_filename, "wb") as file:
                pickle.dump(self._mfile_data, file)
        else:
            raise RuntimeError(f"Unrecognised file format '{_suffix}'")

        self._logger.info("File was written successfully.")


if __name__ in "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("input_mfile")
    parser.add_argument("output_file")

    args = parser.parse_args()

    parser = MFILEParser(args.input_mfile)
    parser.write(args.output_file)
