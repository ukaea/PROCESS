"""
Author: Steven Torrisi, University of Rochester

Data which is stored in the NCDF file from an MFIlE library,
with it's title in the NCDF file listed in parentheses
(Asterisk next to data which is optional and only stored if available):

Title (title)*
Author (author)*
Description (description)*

Dimensions (dimensions)
    -Name
    -Coordinates [From here you can infer # of steps, upper bound, lower bound]

Variables (variables)
    -The dimensions which they are plotted along [same as the ones above]
    -ifail is always included as one of these by default, so there is at least
 one variable of interest
"""

from getpass import getuser
from process.io.mfile import MFile
from numpy import array, transpose
from netCDF4 import Dataset
from datetime import datetime
from process.io.configuration import Config


class NCDFconverter(Config):
    """
    Contains all of the necessary methods to convert MFILEs from a scan
    to a NetCDF file.
    """

    def __init__(self, configfilename="ndscan.json", verbosity=1):
        """
        Specifies which configuration file to read.
        Arguments:
            configfile---> Specifies which configuration file to read.
            verbosity----> If set to 1, will output some information about
            the results. If set to 2, will output results as well as process.

        Class variables:
            currentstep------------->Array to keep track of the current step
                                     along the coordinates
            variablecollector-------> List of lists which collects the values
                                      for each variable of interest
            variablepointers--------> List of strings; the names of each
                                      variable of interest
            output_vars-------------> List of strings; the names of each
                                      variable of interest
            stepstring--------------> String which represents the current step
                                      along the coordinates.
            mfiledir----------------> Directory name in which the MFILES are
                                      stored.

        """
        super().__init__(configfilename)

        self.currentstep = []
        self.variablecollector = []
        self.variablepointers = []
        self.steptuple = ()
        self.axestuple = ()
        self.mfiledir = "MFILES"
        self.stepstring = ""
        self.verbosity = verbosity

        self.configfilename = configfilename
        self.output_vars = self.get("output_vars", default=[])
        if len(self.output_vars) == 0:
            print(
                "No output variables specified in config file.",
                self.configfilename,
                "\nNo NetCDF file created!",
            )
            exit()

        self.time = str(datetime.now())
        self.time = self.time.replace(" ", "-")
        self.time = self.time[0:16]

        self.title = self.get("title", default="NdscanOutput")
        self.author = self.get("author", default=getuser())
        # self.author = self.get("author", default="No author given")
        self.description = self.get("description", default="No description given")

        self.axes = self.get("axes", default=[])

    def convert_mfilelibrary(self):
        """
        The primary method which begins the conversion process on the config
        file named configname.
        This carries out the whole process.
        """

        # Determine roughly but dynamically what level of compression to use.
        scanvars = []
        scansteps = []
        for axis in self.axes:
            scanvars.append(axis["varname"])

            if type(axis["steps"]) is list:
                scansteps.append(len(axis["steps"]))
            else:
                scansteps.append(axis["steps"])

        self.steptuple = tuple(scansteps)
        self.axestuple = tuple(scanvars)

        compressionlevel = self.determine_compression()

        ncfile = self.create_ncfile()

        # A list of lists which will later be converted into n dimensional
        # arrays for each variable.
        # The first list corresponds to the first variable, and so on.
        self.create_dimensions_and_var_slots(ncfile)
        # The first case is done independent of the automation.

        self.do_first_extraction()

        self.start_automated_extraction()

        ncfile.scanvars = ",".join(self.axestuple)
        ncfile.scansteps = self.steptuple

        self.store_variables_from_collector(ncfile, compressionlevel)

        self.print_ncfile_variables(ncfile)

    def increment_current_step(self):
        """
        Increments the current step; compares it against self.steptuple, and
        when an individual axis as at the end of it's bound, rounds up the
        next axis.

        When an 'overflow' occurs along an axis, rounds up.
        For example- If 0th dimension has n steps, when the current step's 0th
        dimension reaches the (n+1)th step,
        increases the 1st dimension step by one, and take the 0th dimension
        back down to the 1st step.
        Then if the 1st dimension overflowed, take that back down to the 1st
        step and tick up the 2nd dimension-- etc...

        Modifies:

        self.currentstep

        WILL END PROGRAM IF:

        The currentstep counter overflows against steptuple. This shouldn't
        happen.

        """

        self.currentstep[0] += 1
        if self.verbosity > 1:
            print("Current step is now ", self.currentstep)

        # Checks to see if increasing the current step by 1 resulted in a
        # dimension overflow
        # if so- reset it to 0, and increment the next dimension by 1.
        for i in range(len(self.currentstep)):
            # Check to see if we are evaluating the last step. If so, an
            # overflow must be treated differently.
            if self.currentstep[-1] == self.steptuple[-1]:
                print(
                    "A serious error has occured, and the total counts\
 have overflowed."
                )
                print(
                    "Currentstep[-1] is",
                    self.currentstep[-1],
                    "and\
 steptuple[-1] is",
                    self.steptuple[-1],
                )
                exit()
            if self.currentstep[i] == self.steptuple[i]:
                self.currentstep[i] = 0
                self.currentstep[i + 1] += 1

    def get_next_mfile(self, currentstep):
        """
        Returns the next MFile in the sequence of mfiles specified by the
        config file.

        Accepts as input the tuple which describes the number of steps for
        each variable.
        Iterates the current step up by one, checks to make sure that the
        current step
        actually describes a coordinate.

        Arguments:

        currentstep----> Represents the current step the Mfile walker is on

        Returns:
        mfile----------> MFile class of the new MFILE to trawl for data

        """

        stepstring = "M."
        for i in range(len(currentstep)):
            stepstring += str(currentstep[i]) + "."
        stepstring += "DAT"

        mfile = MFile(self.mfiledir + "/" + stepstring)

        return mfile

    def get_mfile_variables(self, mfile):
        """
        Searches through the mfile for the name strings in varstofind, and
        returns their corresponding values.

        Arguments:
        mfile--------->MFile class describing the MFILE to be trawled

        Returns:
        varoutput----->List of values corresponding with the variables in
                       varstofind.
        """
        varoutput = [-1] * len(self.output_vars)

        error_status = mfile.data["error_status"].get_scan(-1)
        if error_status < 3:
            for ind, varname in enumerate(self.output_vars):
                varoutput[ind] = mfile.data[varname].get_scan(-1)
        # else:
        #    print('Unsuccessful PROCESS run, MFILE being ignored!')

        return varoutput

    def determine_compression(self):
        """
        From the current ndscan conf file, determine a level of compression
        based on the order of magnitude of the number of data points.

        A simple program which determines if the number of data pieces will
        be greater than 10^4, 10^6, or 10^7.

        Returns:
        compressionlevel-----> An integer which will be valued either 4, 6, 8,
                               or 9.

        """

        # Default level established first
        compressionlevel = 4

        totalsteps = 1
        dataarray = tuple(self.steptuple)

        for i in dataarray:
            totalsteps *= i

        totalpoints = totalsteps * len(self.output_vars)

        # In a rough way, dynamically set the compression based on how much
        # data points are going to be
        # in the file..

        if totalpoints >= 100000:
            compressionlevel = 6
            if totalpoints >= 1000000:
                compressionlevel = 8
                if totalpoints >= 10000000:
                    compressionlevel = 9
        return compressionlevel

    def create_ncfile(self):
        """
        Constructs a NetCDF file from the descriptive information contained
        in an ndscan.conf file.

        The data pieces it searches for are title, output directory,
        description, and author.
        """

        ncfile = Dataset(self.title + ".nc", "w")
        ncfile.title = self.title
        ncfile.time = self.time
        ncfile.description = self.description
        ncfile.author = self.author
        ncfile.ndim = len(self.axes)

        return ncfile

    def create_dimensions_and_var_slots(self, ncfile):
        """
        Creates a dimension for each Axis in the config file, and a spot in
        memory for each variable of interest.

        The variable of interest collector lists generated will act as 'trays'
        to store the variables of interest
        extracted as the MFILEs have their data extracted. They will each be
        coupled with the dimensions specified
        in the scan.

        Arguments:

        ncfile-----> A netCDF file class which will have dimensions and
        variables added to it.

        """
        self.stepstring = "M."
        for i in range(len(self.axes)):

            if type(self.axes[i]["steps"]) is list:
                axissize = len(self.axes[i]["steps"])
            else:
                axissize = int(self.axes[i]["steps"])

            if self.verbosity > 1:
                print(
                    "Creating dimension:",
                    self.axes[i]["varname"],
                    "with size",
                    axissize,
                )

            ncfile.createDimension(dimname=self.axes[i]["varname"], size=axissize)

            self.currentstep.append(0)
            self.stepstring += "0."

        for var in self.output_vars:
            if self.verbosity > 1:
                print("Creating variable", var)
            self.variablecollector.append([])
            self.variablepointers.append(0)

        if "ifail" not in self.output_vars:
            if self.verbosity > 1:
                print("Adding ifail in by default...")

            self.variablecollector.append([])
            self.variablepointers.append(0)
            self.output_vars.append("ifail")

        self.stepstring += "DAT"

    def do_first_extraction(self):
        """
        The first MFILE extraction features more error checks than the
        others for efficiency, and is handled as a seperate function
        before the automated extraction begins.
        """
        mstring = self.mfiledir + "/" + self.stepstring

        try:
            currentmfile = MFile(mstring)
        except FileNotFoundError:
            print("Error: The first Mfile could not be found.")
            print("\t", mstring)
            print("Stopped creation of NetCDF file!")
            exit()

        # Scan the MFILE for the variables of interest
        recentoutput = self.get_mfile_variables(currentmfile)

        # Store output in the collector trays for each variable of interest
        for collectnumber in range(len(recentoutput)):
            self.variablecollector[collectnumber].append(recentoutput[collectnumber])

    def start_automated_extraction(self):
        """
        Initiates the conversion tool's automated trawl across the MFILEs
        for the variables of interest.
        Stores the results in variable collector's respective list for each
        variable.

        Calls:  increment_current_step(steptuple)
                get_next_mfile
                get_mfile_variables

        Modifies:
        self.variablecollector

        """

        if self.verbosity > 1:
            print("Steptuple:", self.steptuple)

        totalsteps = 1
        for i in self.steptuple:
            totalsteps *= i

        for step in range(totalsteps - 1):

            self.increment_current_step()
            currentmfile = self.get_next_mfile(self.currentstep)
            recentoutput = self.get_mfile_variables(currentmfile)

            for collectnumber in range(len(recentoutput)):
                self.variablecollector[collectnumber].append(
                    recentoutput[collectnumber]
                )

    def store_variables_from_collector(self, ncfile, compressionlevel=4):
        """
        Once the data has been extracted from the MFILEs, stores it in the
        NetCDF file.

        Arguments:

        ncfile------------> The NetCDF file to write to
        compressionlevel--> The level of compression to be used on the
        variables

        Modifies:
        self.variablepointers

        """
        steptuple = self.steptuple
        elputpets = steptuple[::-1]  # Reverses the steptuple because the
        # format NCfiles store in is different from Python.

        for i in range(len(self.variablecollector)):
            variablearray = array(self.variablecollector[i], order="C")
            # IMPORTANT- I use column major ordering because I think it
            # makes more sense
            # in the context of the process ndscan suite.
            if len(elputpets) > 1:
                variablearray.shape = elputpets
            self.variablepointers[i] = ncfile.createVariable(
                self.output_vars[i],
                "f4",
                tuple(self.axestuple),
                zlib=True,
                complevel=compressionlevel,
                least_significant_digit=4,
            )

            self.variablepointers[i][:] = transpose(variablearray[:])

        scanvarpointers = []

        for ind, axis in enumerate(self.axes):
            scanvarname = "SCAN" + axis["varname"]
            scanvarpointers.append(0)
            scanvarpointers[ind] = ncfile.createVariable(
                scanvarname, "f4", axis["varname"]
            )
            if type(axis["steps"]) is list:
                scanvarpointers[ind][:] = array(axis["steps"])
            else:
                coords = [
                    axis["lowerbound"]
                    + j
                    * (axis["upperbound"] - axis["lowerbound"])
                    / (axis["steps"] - 1)
                    for j in range(int(axis["steps"]))
                ]
                scanvarpointers[ind][:] = array(coords)

    def print_ncfile_variables(self, ncfile):
        """
        Prints the variables which are currently stored in ncfile. By default,
        runs at the conclusion of each conversion.

        Arguments:
        ncfile-------> Points to the NetCDF file to have variables printed.
        """
        if self.verbosity > 0:
            print("Results:")
            for var in ncfile.variables:
                print("-----------------------")
                print(var)
                print(ncfile.variables[var][:])
