"""

  PROCESS MFILE.DAT IO library

  process.io.mfile.

  James Morris
  CCFE

  Notes:
    + 12/03/2014: Initial version
    + 12/03/2014: Added MFILE variable class
    + 12/03/2014: Added MFILE class for containing all info from file.
    + 12/03/2014: Added ability to read MFILE.DAT into class
    + 12/03/2014: Added ability write MFILE.DAT from class
    + 12/05/2014: Fixed mfile issue with strings in MFILE.DAT with no scans
    + 16/05/2014: Cleaned up MFileVariable
    + 19/05/2014: Cleaned up MFile and put some functions outside class.
    + 12/06/2014: Fixed error handling for "variable not in MFILE" errors
    + 16/06/2014: Fixed library path error; fix in get_scans
    + 24/11/2021: Global dictionary variables moved within the functions
                  to avoid cyclic dependencies. This is because the dicts
                  generation script imports, and inspects, process.

  Compatible with PROCESS version 286

"""

from collections import OrderedDict
import logging
import json
from typing import List, Union

logger = logging.getLogger(__name__)


class MFileVariable(dict):
    """Class for containing a single mfile variable"""

    def __init__(
        self, var_name, var_description, var_unit=None, var_flag=None, *args, **kwargs
    ):
        """
        An object class to contain information (and data values) for a single
        variable from the PROCESS machine readable output file (MFILE.DAT).
        Each class can contain all the scan values for the variable.

        Init defines the variable name and description from arguments:
          var_name --> variable name
          var_description --> variable description/long name

        The class contains the following:
            set_scan    --> function to set a scan number to a given value
            get_scan    --> function to retrieve a given scan
            get_scans   --> function to retrieve all scans.
        """
        self.var_name = var_name
        self.var_description = var_description
        self.var_unit = var_unit
        self.var_flag = var_flag
        self.latest_scan = 0
        super().__init__(*args, **kwargs)
        logger.debug(
            "Initialising variable '{}': {}".format(self.var_name, self.var_description)
        )

    def __getattr__(self, name):
        result = self.get(name)
        # print("Trying to get({}) on {}, {}".format(name, self, id(self)))
        if result:
            return result
        else:
            raise AttributeError(
                "{} object has no attribute {}".format(self.__class__, name)
            )

    def set_scan(self, scan_number, scan_value):
        """Sets the class attribute self.scan# where # is scan number

        Arguments:
          scan_number --> scan number
          scan_value --> value of parameter for scan

        """
        self["scan{:02}".format(scan_number)] = scan_value
        if scan_number > self.latest_scan:
            self.latest_scan = scan_number
        logger.debug(
            "Scan {} for variable '{}' == {}".format(
                scan_number, self.var_name, scan_value
            )
        )

    def get_scan(self, scan_number):
        """Returns the value of a specific scan. For scan = -1 or None the last
        scan is given.

        Arguments:
          scan_number --> scan number to return [-1/None = last scan]

        Returns:
          [single scan requested]
        """

        try:
            if scan_number is None or scan_number == -1:
                return self["scan{:02}".format(self.latest_scan)]
            else:
                return self["scan{:02}".format(scan_number)]
        except KeyError:
            raise  # or substitute with any other exception type you want

    def get_scans(self):
        """Returns a list of scan values in order of scan number

        Returns:
          [List of all scans for variable]

        """
        return [
            v
            for k, v in sorted(
                filter(lambda x: True if "scan" in x[0] else False, self.items())
            )
        ]

    def get_number_of_scans(self):
        """Function to return the number of scans in the variable class"""
        # likely we can just use self.latest_scan, but not guaranteed, so
        # keeping this as it is for now...
        return len([key for key in self.keys() if "scan" in key])

    @property
    def exists(self):
        return True


class MFileErrorClass(object):
    """Error class for handling missing data from MFILE"""

    def __init__(self, item):
        self.item = item
        self.get_scan = self.get_error
        self.get_scans = self.get_error
        self.set_scan = self.get_error
        self.get_number_of_scans = self.get_error

    def get_error(self, *args, **kwargs):
        logger.error("Key '{}' not in MFILE. KeyError! Check MFILE".format(self.item))

        if self.item == "error_status":
            # Missing error_status key means Process exited prematurely, usually
            # due to a "STOP 1"
            raise KeyError(
                "error_status not found in MFILE. Process probably "
                "exited prematurely"
            )
        else:
            return 0

    @property
    def exists(self):
        return False


class MFileDataDictionary(OrderedDict):
    """Class object to act as a dictionary for the data."""

    def __getattr__(self, name):
        result = self.get(name)
        if result:
            return result
        else:
            raise AttributeError(
                "{} object has no attribute {}".format(self.__class__, name)
            )

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return MFileErrorClass(item)


class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if default_factory is not None and not isinstance(default_factory):
            raise TypeError("first argument must be callable")
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            return MFileErrorClass(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = (self.default_factory,)
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy

        return type(self)(self.default_factory, copy.deepcopy(self.items()))

    def __repr__(self):
        return "OrderedDefaultDict(%s, %s)" % (
            self.default_factory,
            OrderedDict.__repr__(self),
        )


class MFile(object):
    def __init__(self, filename="MFILE.DAT"):
        """Class object to store the MFile Objects"""
        logger.info("Creating MFile class for file '{}'".format(filename))
        self.filename = filename
        # self.data = MFileDataDictionary()
        # self.data = OrderedDict()
        self.data = DefaultOrderedDict()
        self.mfile_lines = list()
        self.mfile_modules = dict()
        self.des_name = list()
        self.mfile_modules["Misc"] = list()
        self.current_module = "Misc"
        if filename is not None:
            logger.info("Opening file '{}'".format(self.filename))
            self.open_mfile()
            logger.info("Parsing file '{}'".format(self.filename))
            self.parse_mfile()

    def open_mfile(self):
        """Function to open MFILE.DAT"""
        with open(self.filename, "r", encoding="utf-8") as mfile:
            self.mfile_lines = mfile.readlines()

        for i in range(len(self.mfile_lines)):
            if "*----" in self.mfile_lines[i] or "***" in self.mfile_lines[i]:
                self.mfile_lines = self.mfile_lines[:i]
                self.mfile_end = i
                return

    def parse_mfile(self):
        """Function to parse MFILE.DAT"""
        # for line in (c for c in (clean_line(l) for l in self.mfile_lines if '#' not in l[:2])
        for line in (
            c for c in (clean_line(lines) for lines in self.mfile_lines) if c != [""]
        ):
            self.add_line(line)

    def add_line(self, line):
        """Function to read the line from MFILE and add to the appropriate
        class or create a new class if it is the first instance of it.
        """
        if "#" in line[:2]:
            combined = " ".join(line[1:-1])
            EXCLUSIONS = [
                "Feasible",
                "feasible",
                "Errors",
                "Waveforms",
                "Power Reactor Optimisation Code",
            ]
            if any(exclusion not in combined for exclusion in EXCLUSIONS):
                self.current_module = combined
                self.mfile_modules[self.current_module] = list()

        else:

            var_des = line[0]
            extracted_var_name = sort_brackets(line[1])

            if extracted_var_name == "":
                var_name = var_des
                self.des_name.append(var_name)
            else:
                var_name = extracted_var_name

            if "runtitle" in var_name:
                var_value = " ".join(line[2:])
            else:
                # Pass all value "words"
                var_value = sort_value(line[2:])
            var_unit = get_unit(var_des)
            if len(line) >= 4:
                var_flag = line[3]
            else:
                var_flag = None

            self.mfile_modules[self.current_module].append(var_name)
            self.add_to_mfile_variable(var_des, var_name, var_value, var_unit, var_flag)

    def add_to_mfile_variable(self, des, name, value, unit, flag, scan=None):
        """Function to add value to MFile class for that name/description"""
        if name == "":
            var_key = des.lower().replace("_", " ")
        else:
            var_key = name.lower()

        if var_key in self.data.keys():
            scan_num = scan if scan else (self.data[var_key].get_number_of_scans() + 1)

            # Check for duplicate entries per scan point if there are scans and no scans
            a = len(self.data[var_key].get_scans())
            if "iscan" in self.data.keys():
                b = len(self.data["iscan"].get_scans())
            else:
                b = 1

            if var_key != "iscan":
                if a < b:
                    self.data[var_key].set_scan(scan_num, value)
            else:
                self.data[var_key].set_scan(scan_num, value)
        else:
            var = MFileVariable(
                name, des, unit, var_flag=flag, var_mod=self.current_module
            )
            self.data[var_key] = var
            self.data[var_key].set_scan(1, value)

    def write_to_json(self, keys_to_write=None, scan=-1, verbose=False):
        """Write MFILE object to JSON file"""

        if keys_to_write is None:
            keys_to_write = self.data.keys()

        filename = f"{self.filename}.json"

        dict_to_write = dict()

        if scan == 0:
            for i in range(self.data["rmajor"].get_number_of_scans()):
                sub_dict = {}
                for item in keys_to_write:
                    if self.data[item].get_number_of_scans() == 1:
                        dat_key = -1
                    else:
                        dat_key = i + 1
                    data = self.data[item].get_scan(dat_key)
                    des = self.data[item].var_description.replace("_", " ")
                    if verbose:
                        entry = {"value": data, "description": des}
                    else:
                        entry = data
                    sub_dict[item] = entry
                dict_to_write[f"scan-{i+1}"] = sub_dict
        else:
            for item in keys_to_write:
                # Initialize dat_key properly based on the number of scans
                if self.data[item].get_number_of_scans() == 1:
                    dat_key = -1
                else:
                    dat_key = (
                        scan if scan > 0 else 1
                    )  # Default to scan 1 if not specified
                data = self.data[item].get_scan(dat_key)
                des = self.data[item].var_description.replace("_", " ")
                if verbose:
                    entry = {"value": data, "description": des}
                else:
                    entry = data
                dict_to_write[item] = entry

        with open(filename, "w") as fp:
            json.dump(dict_to_write, fp, indent=4)


def sort_value(value_words: List[str]) -> Union[str, float]:
    """Parse value section of a line in MFILE.

    value_words is a list of strings, which is then parsed.
    :param value_words: value of var in MFILE as list of strings
    :type value_words: List[str]
    :return: string or float representation of value list
    :rtype: Union[str, float]
    """
    if '"' in value_words[0]:
        # First "word" begins with ": return words as single str
        return " ".join(value_words).strip().strip('"').strip()
    else:
        try:
            # Attempt float conversion of first word
            return float(value_words[0])
        except ValueError:
            # Log the exception with details
            logger.exception(f"Can't parse value in MFILE: {value_words}")
            # Return the original string as a fallback
            return " ".join(value_words).strip()


def sort_brackets(var):
    """Function to sort bracket madness on variable name."""
    if var != "":
        tmp_name = var.lstrip("(").split(")")
        if len(tmp_name) > 2:
            return tmp_name[0] + ")"
        else:
            return tmp_name[0]
    else:
        return ""


def clean_line(line):
    """Cleans an MFILE line into the three parts we care about"""
    cleaned_line = [item.strip("_ \n") for item in line.split(" ") if item != ""]
    return cleaned_line


def search_keys(dictionary, variable):
    """Searches the dictionary keys for matches to the variable name
    'variable' in arguments.

    Puts everything into lower case before searching.

    Argument:
      dictionary --> dictionary to search in
      variable --> variable name to search for

    Returns:
      matches --> List of matches to the searched for variable

    """
    matches = []
    for key in dictionary.keys():
        if variable.lower() in key.lower():
            matches.append(key)
    return matches


def search_des(dictionary, description):
    """Searches the dictionary descriptions for matches to the description
    'description' from the arguments.

    Puts everything into lower case before searching.

    Argument:
      dictionary --> dictionary to search in
      variable --> variable name to search for

    Returns:
      matches --> List of matches to the searched for description

    """
    descriptions = [
        dictionary[key].var_descript.lower() for key in dictionary.data.keys()
    ]
    matches = []
    for item in descriptions:
        if description.lower() in item.lower():
            matches.append(item)
    return matches


def get_unit(variable_desc):
    """Returns the unit from a variable description if possible, else None."""
    candidate = variable_desc.rsplit("_", 1)[-1]
    if candidate.startswith("(") and candidate.endswith(")"):
        return candidate[1:-1]
    else:
        return None


def is_number(val):
    """Check MFILE data entry"""
    try:
        float(val)
        return True
    except ValueError:
        pass

    return False


def test(f):
    """Testing

    :param f: file to test"""

    try:
        # print(m.data["rmajor"].get_scans())
        return True
    except Exception:
        return False
    return True


# if __name__ == "__main__":
# filename = sys.argv[1]
# m = MFile(filename)
# print(m.data["rmajor"].get_number_of_scans())
# print(m.data["rmajor"].get_scans())
# print(m.data["rmajor"].get_scan(2))
