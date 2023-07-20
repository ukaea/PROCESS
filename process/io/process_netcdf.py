"""
PROCESS I/O library NetCDF API.

Used for reading/writing NetCDF PROCESS data.
"""

import os
import re
import sys
from sys import stderr
import numpy as np
from netCDF4 import Dataset
from process.io.mfile import MFile, MFileErrorClass
from process.io.in_dat import InDat


NAME_MAPPINGS = {"/": "_slash_", "*": "_asterisk_", ".": "_dot_"}
METADATA = ("procver", "date", "time", "username", "isweep", "nsweep")


class NetCDFWriter(object):

    """Takes PROCESS data and writes it to a NetCDF file."""

    def __init__(self, netcdf_filename, append=True, overwrite=False):
        self.netcdf_filename = os.path.abspath(netcdf_filename)
        self._append = append
        self._overwrite = overwrite
        self.root = None

    def __enter__(self):
        self._open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._close()

    def _open(self):
        try:
            mode = (
                "a" if (os.path.exists(self.netcdf_filename) and self._append) else "w"
            )
            self.root = Dataset(self.netcdf_filename, mode, clobber=self._overwrite)
        except RuntimeError:
            raise OSError(
                "Cannot create {} - file may already "
                "exist".format(self.netcdf_filename)
            )

    def _close(self):
        """Correctly flush data to NetCDF file and close for writing."""
        try:
            self.root.close()
        except AttributeError:
            print("File not initially opened by NetCDFWriter")

    def _store_variable(self, var_name, value, var_group):

        if var_name == "bounds":
            boundfmt = "bound{}({})"
            var_type = "f8"
            for boundno, bound_dict in value.items():
                for boundtype, boundval in bound_dict.items():
                    bvar_name = boundfmt.format(boundtype, boundno)
                    bvar_val = boundval
                    stored_val = np.array(bvar_val, dtype=var_type)
                    stored_var = var_group.createVariable(bvar_name, var_type)
                    try:
                        stored_var[:] = stored_val
                    except RuntimeError as err:
                        print("Please debug this!", file=stderr)
                        print(var_name, value, var_group, file=stderr)
                        print(err, file=stderr)
        elif var_name in ["fimp", "zref"]:
            arrfmt = "{}({})"
            var_type = "f8"
            for i, a_value in enumerate(value):
                # fortran starts counting at 1!
                avar_name = arrfmt.format(var_name, i + 1)
                avar_val = a_value
                stored_val = np.array(avar_val, dtype=var_type)
                stored_var = var_group.createVariable(avar_name, var_type)
                stored_var[:] = stored_val

        else:
            try:
                var_type = "f8"
                stored_val = np.array(value, dtype=var_type)
            except ValueError:
                var_type = "str"
                stored_val = np.array(value, dtype=var_type)

            stored_var = var_group.createVariable(var_name, var_type)
            try:
                stored_var[:] = stored_val
            except RuntimeError as err:
                print("Please debug this!", file=stderr)
                print(var_name, value, var_group, file=stderr)
                print(err, file=stderr)
                exit()

    def handle_unknown_vars(self, save_vars, ignore_unknowns, source_data, source_type):
        """
        Common check for handling unspecified variables.
        save_vars       : string/list  : Variables to store
        ignore_unknowns : Bool         : Ignores unknown variables
        """
        if isinstance(save_vars, list):
            unknown_vars = set(save_vars) - set(source_data.keys())

            if any(unknown_vars) and ignore_unknowns:
                print(
                    "Cannot save these variables (not provided in {}"
                    "instance): {}".format(source_type, unknown_vars)
                )
            elif any(unknown_vars):
                raise KeyError(
                    "Cannot save these variables (not in provided "
                    "{} instance): {}".format(source_type, unknown_vars)
                )
            else:
                pass

    def handle_unknowns(
        self,
        save_vars,
        ignore_unknowns,
        source_data1,
        source_type1,
        source_data2,
        source_type2,
    ):
        """Common check for handling unspecified variables."""
        if isinstance(save_vars, list):
            save_vars1 = set(source_data1.keys()) & set(save_vars)
            unknown_vars = set(save_vars) - save_vars1
            save_vars2 = set(source_data2.keys()) & unknown_vars
            unknown_vars = unknown_vars - save_vars2

            if any(unknown_vars) and ignore_unknowns:
                print(
                    "Cannot save these variables (neither provided in {}"
                    "instance nor in {} instance): {}".format(
                        source_type1, source_type2, unknown_vars
                    )
                )
            elif any(unknown_vars):
                raise KeyError(
                    "Cannot save these variables (neither provided",
                    " in {} instance".format(source_type1),
                    " nor in {} instance): {}".format(source_type2, unknown_vars),
                )
            else:
                pass

            return list(save_vars1), list(save_vars2)

        else:
            return "all", "all"

    def write_in_data(self, in_dat, run_id, save_vars="all", ignore_unknowns=False):
        """
        Write the provided InDat instance out as an input group in the NetCDF
        associated with the given run_id.

        Raises KeyError if save_vars is a list of variables and any of which are
        not available for extraction from the provided InDat instance and
        ignore_unknowns is False.
        """

        in_data = in_dat.data
        indat_vars = self.root.createGroup("input_{}".format(run_id))

        # Don't need to include metadata here as it should be stored with the
        # output data group...
        keys = []

        # Stop/warn when there are requested to-save variables that don't exist
        self.handle_unknown_vars(save_vars, ignore_unknowns, in_data, "InDat")

        for k in in_data.keys():
            if k.endswith("."):
                continue

            # Swap illegal characters in key with substitutes in NAME_MAPPINGS
            rep_key = dict((re.escape(k), val) for k, val in NAME_MAPPINGS.items())
            pattern = re.compile("|".join(rep_key.keys()))
            # Swaps all illegal characters in one go
            if save_vars == "all" or k in save_vars:
                keys.append(pattern.sub(lambda m: rep_key[re.escape(m.group(0))], k))
            else:
                pass

        for key in keys:
            var_group = indat_vars.createGroup(key)
            val = in_data[key]
            var_group.v_type = val.v_type
            var_group.comment = val.comment

            self._store_variable(val.name, val.value, var_group)

    def write_mfile_data(self, mfile, run_id, save_vars="all", ignore_unknowns=False):
        """
        Write the provided MFile instance out as a run within the NetCDF.

        Raises KeyError if save_vars is a list of variables and any of which are
        not available for extraction from the provided MFile instance and
        ignore_unknowns is False.
        """

        mfile_data = mfile.data
        mfile_vars = self.root.createGroup("output_{}".format(run_id))

        # Make sure we include metadata
        keys = []
        keys += METADATA

        # Stop/warn when there are requested to-save variables that don't exist
        self.handle_unknown_vars(save_vars, ignore_unknowns, mfile_data, "MFile")

        # for k, val in mfile_data.items():
        for k in mfile_data.keys():
            if k.endswith("."):
                continue

            # Swap illegal characters in key with substitutes in NAME_MAPPINGS
            rep_key = dict((re.escape(k), val) for k, val in NAME_MAPPINGS.items())
            pattern = re.compile("|".join(rep_key.keys()))
            # Swaps all illegal characters in one go
            if save_vars == "all" or k in save_vars:
                keys.append(pattern.sub(lambda m: rep_key[re.escape(m.group(0))], k))
            else:
                pass

        for key in keys:
            var_group = mfile_vars.createGroup(key)
            var = mfile_data[key]
            if isinstance(var, MFileErrorClass):
                continue
            else:
                var_group.description = var["var_description"]
                var_group.group_name = var["var_name"]
                var_group.unit = (
                    var["var_unit"] if (var["var_unit"] is not None) else "None"
                )
            possible_scans = dict(
                (scan, scanval)
                for scan, scanval in mfile_data[key].items()
                if "scan" in scan
            )

            # uses only latest scan
            highest_scan = 0
            latest_scan = None
            for scan_k in possible_scans.keys():
                scan_num = int(scan_k.strip("scan"))
                if scan_num > highest_scan:
                    highest_scan = scan_num
                    latest_scan = scan_k
            self._store_variable(latest_scan, possible_scans[latest_scan], var_group)

    def write_data(self, mfile, in_dat, run_id, save_vars="all", ignore_unknowns=False):
        """
        Write the provided MFile and InDat instances out as into the NetCDF
        associated with the given run_id.
        mfile           : MFile object : input MFILE
        in_dat          : InDat object : input IN.DAT file
        run_id          : Integer      : Run ID
        save_vars       : string/list  : Variables to store
        ignore_unknowns : Bool         : Ignores unknown variables
        """

        mfile_data = mfile.data
        in_data = in_dat.data
        mfile_vars = self.root.createGroup("output_{}".format(run_id))
        indat_vars = self.root.createGroup("input_{}".format(run_id))

        # Make sure we include metadata from MFILE
        mfile_keys = []
        mfile_keys += METADATA
        indat_keys = []

        # Stop/warn when there are requested to-save variables that don't exist
        save_vars_mfile, save_vars_indat = self.handle_unknowns(
            save_vars, ignore_unknowns, mfile_data, "MFile", in_data, "InDat"
        )

        # MFILE
        for k in mfile_data.keys():
            if k.endswith("."):
                continue

            # Swap illegal characters in key with substitutes in NAME_MAPPINGS
            rep_key = dict((re.escape(k), val) for k, val in NAME_MAPPINGS.items())
            pattern = re.compile("|".join(rep_key.keys()))
            # Swaps all illegal characters in one go
            if save_vars_mfile == "all" or k in save_vars_mfile:
                mfile_keys.append(
                    pattern.sub(lambda m: rep_key[re.escape(m.group(0))], k)
                )
            else:
                pass

        for key in mfile_keys:
            var_group = mfile_vars.createGroup(key)
            var = mfile_data[key]

            if isinstance(var, MFileErrorClass):
                continue
            else:

                var_group.description = var.var_description
                var_group.group_name = var.var_name
                var_group.unit = var.var_unit if (var.var_unit is not None) else "None"

                possible_scans = dict(
                    (scan, scanval)
                    for scan, scanval in mfile_data[key].items()
                    if "scan" in scan
                )

            # using latest scan only
            highest_scan = 0
            latest_scan = None
            for scan_k in possible_scans.keys():
                scan_num = int(scan_k.strip("scan"))
                if scan_num > highest_scan:
                    highest_scan = scan_num
                    latest_scan = scan_k
            self._store_variable(latest_scan, possible_scans[latest_scan], var_group)

        # INDAT
        for k in in_data.keys():
            if k.endswith("."):
                continue

            # Swap illegal characters in key with substitutes in NAME_MAPPINGS
            rep_key = dict((re.escape(k), val) for k, val in NAME_MAPPINGS.items())
            pattern = re.compile("|".join(rep_key.keys()))
            # Swaps all illegal characters in one go
            if save_vars_indat == "all" or k in save_vars_indat:
                indat_keys.append(
                    pattern.sub(lambda m: rep_key[re.escape(m.group(0))], k)
                )
            else:
                pass

        for key in indat_keys:
            var_group = indat_vars.createGroup(key)
            val = in_data[key]
            var_group.v_type = val.v_type
            var_group.comment = val.comment

            self._store_variable(val.name, val.value, var_group)


class NetCDFReader(object):

    """Capable of reading NetCDF PROCESS data files and returning [...]."""

    def __init__(self, netcdf_filename):
        self.netcdf_filename = os.path.abspath(netcdf_filename)

    def __enter__(self):
        """Open NetCDF file and provide with statement usage."""
        self._open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close NetCDF file and provide with statement usage."""
        self._close()

    def _open(self):
        """Open the NetCDF file for reading."""
        try:
            self.root = Dataset(self.netcdf_filename, "r")
        except RuntimeError:
            raise FileNotFoundError("Cannot read" " {}".format(self.netcdf_filename))

    def _close(self):
        """Correctly close NetCDF file handle."""
        try:
            self.root.close()
        except AttributeError:
            print("File not initially opened by NetCDFReader")

    def _get_mfiledict(self, path):
        """Return a dict of variable, {datadict} from the output parameters."""
        mfile_dict = dict()
        try:
            mfile_data = self.root.groups[path]
        except KeyError:
            ("Cannot access {} in " "{}".format(path, self.netcdf_filename))
        for group in mfile_data.groups.values():
            for variable in group.variables.values():
                mfile_dict[group.group_name] = variable.getValue()

        return mfile_dict

    def _get_mfile(self, path):
        m_file = MFile(filename=None)
        try:
            mfile_data = self.root.groups[path]
        except KeyError:
            ("Cannot access {} in " "{}".format(path, self.netcdf_filename))
        for group in mfile_data.groups.values():
            for var_name, variable in group.variables.items():
                scan_num = None
                if "scan" in var_name:
                    scan_num = int(re.search(r"\d+", var_name).group())
                m_file.add_to_mfile_variable(
                    group.description,
                    group.group_name,
                    variable.getValue(),
                    group.unit,
                    scan_num,
                )
        return m_file

    def _get_indatdict(self, path):

        """Return a dict of variable, {datadict} from the input parameters."""

        indat_dict = dict()
        try:
            input_data = self.root.groups[path]
        except KeyError:
            ("Cannot access {} in " "{}".format(path, self.netcdf_filename))
        for group in input_data.groups.values():
            for var_name, variable in group.variables.items():
                indat_dict[var_name] = variable.getValue()

        return indat_dict

    def _get_indat(self, path):

        """Returns an indat instance, though exceptions are not treated
        and therefore ignored!"""

        indat = InDat(filename=None)
        try:
            input_data = self.root.groups[path]
        except KeyError:
            ("Cannot access {} in " "{}".format(path, self.netcdf_filename))
        for group in input_data.groups.values():
            for var_name, variable in group.variables.items():
                indat.add_parameter(var_name, variable.getValue())

        return indat

    def get_run_mfile(self, run_id=1):
        """Return the run_id data as an MFile instance.

        run_id is the ID number of the PROCESS output to retrieve. If it cannot
        be found, a KeyError is raised.
        """
        return self._get_mfile("output_{}".format(run_id))

    def get_run_indatdict(self, run_id=1):
        """Return the run_id data as indat_dict.

        run_id is the ID number of the PROCESS output to retrieve. If it cannot
        be found, a KeyError is raised.
        """
        return self._get_indatdict("input_{}".format(run_id))

    def get_run_indat(self, run_id=1):
        """Return the run_id data as indat instance.

        run_id is the ID number of the PROCESS output to retrieve. If it cannot
        be found, a KeyError is raised.
        """
        return self._get_indat("input_{}".format(run_id))

    def get_run_dicts(self, run_id=1):
        """returns dict with input and output data"""
        mfiledict = self._get_mfiledict("output_{}".format(run_id))
        indatdict = self._get_indatdict("input_{}".format(run_id))
        indatdict.update(mfiledict)  # in conflict takes mfile value!
        return indatdict

    def runs(self, start_run=1):
        """Generator to provide each run, starting from the ID start_run."""
        for run in self.root.groups.keys():
            try:
                run_id = int(run.split("output_")[1])
            except IndexError:
                continue
            if run_id >= start_run:
                yield self._get_mfile("output_{}".format(run_id))
            else:
                pass

    def run_dicts(self, start_run=1):
        """Generator to provide dicts of each run, starting from the ID start_run."""

        for run in self.root.groups.keys():
            try:
                run_id = int(run.split("output_")[1])
            except IndexError:  # for input_run_id to avoid double counting!
                continue
            if run_id >= start_run:
                mfiledict = self._get_mfiledict("output_{}".format(run_id))
                indatdict = self._get_indatdict("input_{}".format(run_id))
                indatdict.update(mfiledict)  # in conflict takes mfile value!
                yield indatdict
            else:
                pass


if __name__ == "__main__":

    mf = MFile("/home/edwardsj/example_MFILE.DAT")

    if sys.argv[1] == "write":
        # Writer example - now uses context manager (i.e. with statement)
        with NetCDFWriter(
            "/home/edwardsj/tester.nc", append=True, overwrite=False
        ) as ncdf_writer:
            ncdf_writer.write_mfile_data(mf, 3, save_vars="all", latest_scan_only=True)
    elif sys.argv[1] == "read":
        # Reader example - gives back MFile instances
        with NetCDFReader("/home/edwardsj/tester.nc") as ncdf_reader:
            # Get an individual run
            returned_mf = ncdf_reader.get_run_mfile(4)
            print(returned_mf.data)
            # Get multiple runs in a loop, starting at run 1
            # for mfile in ncdf_reader.runs(start_run=2):
            #     print(mfile)
