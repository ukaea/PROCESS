import logging
import os
import re
from collections import OrderedDict, abc
from pathlib import Path
from typing import Any

MFILE_END = "# Copy of PROCESS Input Follows #"
VETO_LIST = [" # PROCESS found a feasible solution #"]
HEADER_MAPPING = {"Power Reactor Optimisation Code": "metadata"}


class MFILEParser(abc.MutableMapping):
    """Parse an MFILE and extract output values."""

    def __init__(self, input_mfile: str = ""):
        self._input_file = input_mfile
        self._mfile_data: OrderedDict = OrderedDict()
        self._logger = logging.getLogger(self.__class__.__name__)
        if self._input_file:
            self.parse(self._input_file)

    def __iter__(self):
        for group in self._mfile_data:
            yield from self._mfile_data[group]

    def __len__(self):
        return sum(len(self._mfile_data[g]) for g in self._mfile_data)

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
        :
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
        lines :
            list of file lines to search
        search_str:
            search term to look for

        """
        return [i for i, line in enumerate(lines) if search_str in line]

    def _find_var_val_from_str(self, value_str: str) -> Any:
        """Convert a string variable to float, int etc.

        This function parsers values given within the MFILE removing other
        unneeded information such as the specific typing in PROCESS

        Parameters
        ----------
        value_str :
            value as a string
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
        lines :
            list of file lines to be parsed
        """
        # Compile regex for converting underscores which are spaces into
        # a space character
        space_re = r"(\_{5,})"
        var_re = r"(\([a-z0-9\-\+\*\/\_\%\]\[]+\))"
        # TODO remove underscores
        # Extract lines from the given line set that follow the variable
        # statement convention of 'desc_______(varname)________ value'
        lines_ = [line for line in lines if re.findall(var_re, line)]

        # Remove extra symbols such as quotation marks and split line into
        # the three required components using regex
        lines_ = [
            [
                i.replace('"', "").replace("`", "").strip()
                for i in re.split(space_re, line)
                if not (re.findall(space_re, i)) and i.strip()
            ]
            for line in lines_
        ]

        # If there are not three components in a given line, try splitting
        # the components present by ' ' instead and append
        for i, line in enumerate(lines_):
            if len(line) != 3:
                new_line = []
                for element in line:
                    if " " in element:
                        new_line += element.split(" ")
                    else:
                        new_line += element
                lines_[i] = new_line[:3]

        # Use an ordered dictionary to match ordering in MFILE
        vars_dict = OrderedDict()

        # Iterate through the resultant line sets and tidy them a little
        # finally creating a dictionary entry for each with the required
        # information
        for line in lines_:
            var_key = line[1][1:-1]
            var_key = var_key.replace("%", ".")
            if not var_key:
                continue
            value = line[2]
            desc = line[0].replace("_-_", "-").replace("_", " ")
            desc = desc.title().strip()
            desc = desc.replace('"', "")
            desc = re.sub(r"\s{2,}", " ", desc)
            if var_key in vars_dict:
                if not isinstance(vars_dict[var_key], list):
                    vars_dict[var_key]["value"] = [vars_dict[var_key]["value"]]
                vars_dict[var_key]["value"].append(self._find_var_val_from_str(value))
            else:
                vars_dict[var_key] = {
                    "description": desc,
                    "value": self._find_var_val_from_str(value),
                }

        return vars_dict

    def parse(self, mfile_addr: str) -> dict:
        """Parse an MFILE and extract output values.

        Parameters
        ----------
        mfile_addr :
            address of MFILE to parse

        Returns
        -------
        :
            dictionary of output values

        """
        if not os.path.exists(mfile_addr):
            raise FileNotFoundError(
                f"Could not open MFILE '{mfile_addr}', file does not exist."
            )

        self._logger.info("Parsing MFILE: %s", mfile_addr)

        with open(mfile_addr) as f:
            lines = f.readlines()

        end_of_output = self._line_string_search(lines, MFILE_END)[0]

        self._logger.info("Extracting file headers")
        header_indexes = [
            i for i, line in enumerate(lines) if line.strip() and i < end_of_output
        ]

        header_indexes = [
            i
            for i in header_indexes
            if lines[i].strip()[0] == "#" and not any(k in lines[i] for k in VETO_LIST)
        ]

        # Gets rid of multi-headers, taking the last one
        header_indexes = [i for i in header_indexes if i + 1 not in header_indexes]

        self._logger.info("Retrieving output variable values")
        # Iterate through the file headers processing the "block" between them
        # extracting variable values. Where duplicate headers are found assume
        # that a parameter sweep is occuring and append values in lists
        for i in range(len(header_indexes) - 1):
            key = lines[header_indexes[i]].replace("#", "").strip()

            new_vals = self._get_values(
                lines[header_indexes[i] + 1 : header_indexes[i + 1]]
            )

            # The iscan variable is always present at start of sweep
            # no matter what the first header in an iteration is
            # need to move it into metadata
            first_key = lines[header_indexes[0]]
            check_iscan = self._mfile_data and "iscan" in new_vals
            check_iscan = check_iscan and key != first_key
            if check_iscan:
                first_key = first_key.replace("#", "").strip()
                iscan_var = self._mfile_data[first_key]["iscan"]["value"]
                if not isinstance(iscan_var, list):
                    self._mfile_data[first_key]["iscan"]["value"] = [iscan_var]
                self._mfile_data[first_key]["iscan"]["value"].append(
                    new_vals["iscan"]["value"]
                )
                del new_vals["iscan"]

            # Add header to dictionary of not present
            if key not in self._mfile_data:
                self._mfile_data[key] = new_vals

            # If header already present, iterate through member parameters
            # appending the new values to each
            else:
                for param, var_dict in self._mfile_data[key].items():
                    if param not in new_vals:
                        self._logger.warning(
                            f"Expected parameter '{param}' in sweep, "
                            "but could not find entry"
                            " for this iteration"
                        )
                        continue
                    value = new_vals[param]["value"]
                    if not isinstance(var_dict["value"], list):
                        self._mfile_data[key][param]["value"] = [
                            self._mfile_data[key][param]["value"],
                            value,
                        ]
                    else:
                        # Need to check if the find variables function
                        # returned a single value for the parameter or multiple
                        # and handle the cases
                        if not isinstance(new_vals[param]["value"], list):
                            self._mfile_data[key][param]["value"].append(value)
                        else:
                            self._mfile_data[key][param]["value"] += value

        self._logger.info("Creating output dictionaries")
        # Remove any cases where there are no parameters under a given header
        self._mfile_data = {k: v for k, v in self._mfile_data.items() if v}

        # Use underscore keys and tidy them to be more computationally friendly
        def _key_update(key):
            key_ = key.lower()
            key_ = key_.replace(" ", "_")
            for sym in [":", "(", ")", "/"]:
                key_ = key_.replace(sym, "")
            return key_.replace("__", "_")

        # Apply header mappings and tidy headers
        self._mfile_data = {
            _key_update(k) if k not in HEADER_MAPPING else HEADER_MAPPING[k]: v
            for k, v in self._mfile_data.items()
        }

        if not self._mfile_data or len(self._mfile_data) == 0:
            raise AssertionError("Failed to extract data from given MFILE")

        # Only run iscan check if iscan exists
        try:
            first_key = next(iter(self._mfile_data.keys()))
            second_key = list(self._mfile_data.keys())[1]
            second_key_fp = list(self._mfile_data[second_key])[8]
            iscan_arr = self._mfile_data[first_key]["iscan"]["value"]
            test_param = self._mfile_data[second_key][second_key_fp]["value"]
            if len(test_param) != iscan_arr[-1]:
                print(test_param)
                raise AssertionError(
                    "Failed to retrieve all parameter sweep values, "
                    f"expected {iscan_arr[-1]} values for '{second_key}:{second_key_fp}' and got {len(test_param)}"
                )
        except KeyError:
            pass

        self._logger.info("Extraction completed successfully")
        return self._mfile_data

    def write(self, output_filename: str):
        """Write output to file.

        Parameters
        ----------
        output_filename : str
            path of output file, file type is determined from the type and can
            be '.toml', '.yml', '.pckl', '.json'
        """
        self._logger.info("Writing to output file '%s'", output_filename)

        suffix = os.path.splitext(output_filename)[1].lower()

        if suffix == ".toml":
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
            doc = tomlkit.document()
            doc.add(tomlkit.comment("PROCESS Run Output"))
            for group_name, data in self._mfile_data.items():
                new_dict = {}
                for var_name, var_data in data.items():
                    new_dict[var_name] = var_data["value"]
                header = group_name.replace("_", " ").title()
                ls = int((75 - len(header)) / 2)
                rs = 75 - len(header) - ls
                header = ls * "-" + " " + header + " " + rs * "-"
                doc.add(tomlkit.comment(header))
                doc.add(group_name, new_dict)
                doc.add(tomlkit.nl())
                doc.add(tomlkit.nl())

            for group_name, data in self._mfile_data.items():
                for var_name in data:
                    doc[group_name][var_name].comment(
                        self._mfile_data[group_name][var_name]["description"]
                    )

            Path(output_filename).write_text(tomlkit.dumps(doc))
        elif suffix == ".json":
            # If file suffix is JSON
            self._logger.info("Output will be JSON file.")
            import json

            with open(output_filename, "w") as file:
                json.dump(self._mfile_data, file)
        elif suffix in [".yml", ".yaml"]:
            self._logger.info("Output will be YAML file.")
            # If file suffix is YAML
            import yaml

            with open(output_filename, "w") as file:
                yaml.dump(self._mfile_data, file)
        elif suffix == ".pckl":
            self._logger.info("Output will be Pickle file.")
            # If file suffix is Pickle
            import pickle

            with open(output_filename, "wb") as file:
                pickle.dump(self._mfile_data, file)
        else:
            raise RuntimeError(f"Unrecognised file format '{suffix}'")

        self._logger.info("File was written successfully.")


if __name__ in "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("input_mfile")
    parser.add_argument("output_file")

    args = parser.parse_args()

    parser = MFILEParser(args.input_mfile)
    parser.write(args.output_file)
