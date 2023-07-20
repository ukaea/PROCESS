""" Visualisation Class for Ndscan utility """


from numpy import meshgrid, amax, amin, std, sum
import matplotlib.pyplot as pl
from netCDF4 import Dataset


class VisUtility(object):

    """Visualisation class containing Netcdf files from ndscan tool"""

    def __init__(self, ncfilename="Demonstrationdata.nc"):
        """
        Opens the ncfile at ncfilename and initializes its' variables
        into memory.

        """
        try:
            self.ncfile = Dataset(ncfilename, "r")
        except RuntimeError:
            print("Error: NetCDF file does not exist, specify with -f!")
            exit()
        scanvars = self.ncfile.scanvars
        self.scanvars = scanvars.split(",")
        self.currentslice = None
        self.stepslist = []

    def make_contour_plot(self, xdim, ydim, zvar, masterconfigs):
        """
        Generates a contour plot using matplotlib and the arguments given.

        Arguments:
        xdim-------------> String representing the X scan dimension
        ydim-------------> String representing the Y scan dimension
        zvar-------------> String which represents the Z variable to plot.
        masterconfigs----> A configuration dict with some details.


        NOTE! XVAR AND Y VAR SHOULD BE A STRING REPRESENTING SCANNING VARIABLES,
        ZVAR SHOULD BE A STRING REPRESENTING A VARIABLE OF INTEREST.
        """

        print("xdim is", xdim, "ydim is", ydim)
        if xdim[0:4] != "SCAN":
            xdim = "SCAN" + xdim
        if ydim[0:4] != "SCAN":
            ydim = "SCAN" + ydim
        xcoords = self.ncfile.variables[xdim][:]
        ycoords = self.ncfile.variables[ydim][:]
        xarr, yarr = meshgrid(xcoords, ycoords)
        if self.currentslice is not None:
            zarr = self.ncfile.variables[zvar][self.currentslice]
        else:
            zarr = self.ncfile.variables[zvar][:]
        zarr = zarr.transpose()

        conp = pl.contour(xarr, yarr, zarr, 10)
        pl.xlabel(masterconfigs["xname"])
        pl.ylabel(masterconfigs["yname"])
        pl.title(masterconfigs["zname"])
        pl.clabel(conp, inline=1, fontsize=10)
        pl.show()

    def make_scatter_plot(self, xdim, ydim, zvar, masterconfigs):
        """
        Generates a contour plot using matplotlib and the arguments given.

        Arguments:
        xdim-------------> String representing the X scan dimension
        ydim-------------> String representing the Y scan dimension
        zvar-------------> String which represents the Z variable to plot.
        masterconfigs----> A configuration dicts which is not implemented.

        NOTE! XVAR AND Y VAR SHOULD BE A NUMBER REPRESENTING SCANNING VARIABLES,
        ZVAR SHOULD BE A NUMBER REPRESENTING A VARIABLE OF INTEREST.
        ."""

        print("xdim is", xdim, "ydim is", ydim)
        if xdim[0:4] != "SCAN":
            xdim = "SCAN" + xdim
        if ydim[0:4] != "SCAN":
            ydim = "SCAN" + ydim
        xcoords = self.ncfile.variables[xdim]
        ycoords = self.ncfile.variables[ydim]
        xarr, yarr = meshgrid(xcoords, ycoords)
        if self.currentslice is not None:
            zarr = self.ncfile.variables[zvar][self.currentslice]
        else:
            zarr = self.ncfile.variables[zvar][:]

        zarr = zarr.transpose()

        pl.scatter(xarr, yarr, c=zarr, marker="s")
        pl.xlabel(masterconfigs["xname"])
        pl.ylabel(masterconfigs["yname"])
        pl.title(masterconfigs["zname"])
        pl.colorbar()
        pl.show()

    def slicescanner(self, xvarnum, yvarnum, zvar):

        """
        Figures out what the 2 dimensional slices of a >=3 dimensional space
        looks like, then asks the user
        to pick one to use for later 2d plots.

        How it works: A dimension tuple looks like
         (x,y,z,p,q) === (0,1,2,3,4)
         From there we identify which two tuples we are scanning along, say x
         and p or 0 and 3
         Then we iterate through them piece by piece so extract the variable
         from the ncscan designated by zvar,
         then do zvar[x,0,0,p,0] from netcdf, then that's a 2d array!
         and from there do analytics, store the results, and grab the next 2d
         array.
         And so on.

        Modifies:
        stepslist


        """
        print("Activating slice scanner for xvarnum", xvarnum, "and yvarnum", yvarnum)
        superlatives = {
            "ZenithSlice": None,
            "NadirSlice": None,
            "MaximalSlice": None,
            "MinimalSlice": None,
            "SteepestSlice": None,
            "FlattestSlice": None,
        }
        sliceindexlist = []
        maxlist = []
        minlist = []
        sumlist = []
        stdlist = []

        steptuple = tuple(self.ncfile.scansteps)
        stepslist = []

        for i in range(len(steptuple)):
            if i != xvarnum and i != yvarnum:
                stepslist.append(0)
            else:
                stepslist.append(list(range(steptuple[i])))
        self.stepslist = stepslist
        # We just reset it so that it's [x][0][0][p][0]...[0] for our
        # dimensions.
        # Now we need this to iterate by incrementing the leftmost dimension
        # until it is equal to it's steptuple
        # value-- then we reduce it back to zero and add 1 to the one on the
        # left. once stepslist's non-list elements
        # are equal to the steptuple, then we are good.

        def keep_going(steptuple, stepslist):

            for index in range(len(stepslist)):
                if type(stepslist[index]) is list:

                    continue
                else:
                    stepslist[index] += 1

                if stepslist[index] < steptuple[index]:
                    return stepslist

                if stepslist[index] == steptuple[index] and index != len(stepslist) - 1:

                    stepslist[index] = 0
                    for countahead in range(len(stepslist))[index + 1 : :]:

                        if type(stepslist[countahead]) is int:
                            stepslist[countahead] += 1

                            countback = -1
                            if type(stepslist[-1]) is list:
                                while type(stepslist[countback]) is not int:
                                    countback -= 1

                            if stepslist[countback] == steptuple[countback]:
                                print("Done!")
                                return False

                            else:
                                return stepslist
                        else:
                            continue
                    print("Done!")

                    return False

                print("Done!")
                return False

        # A little hack to contain all of the machinery of analyzing slices
        # within the loop..
        while stepslist:

            sliceindexlist.append(list(stepslist))
            temp = self.ncfile.variables[zvar][stepslist]
            maxlist.append(amax(temp))
            minlist.append(amin(temp))
            sumlist.append(sum(temp))
            stdlist.append(std(temp))
            print(
                "Analytics: Max:",
                amax(temp),
                "Min:",
                amin(temp),
                "Sum:",
                sum(temp),
                "Std:",
                std(temp),
            )
            print(temp)
            stepslist = keep_going(steptuple, stepslist)

        print("Your options:")
        for i in range(len(sliceindexlist)):
            print(i, ":", sliceindexlist[i])

        superlatives["ZenithSlice"] = sliceindexlist[maxlist.index(max(maxlist))]
        superlatives["NadirSlice"] = sliceindexlist[minlist.index(min(minlist))]

        superlatives["MaximalSlice"] = sliceindexlist[sumlist.index(max(sumlist))]
        superlatives["MinimalSlice"] = sliceindexlist[sumlist.index(min(sumlist))]

        superlatives["SteepestSlice"] = sliceindexlist[stdlist.index(max(stdlist))]
        superlatives["FlattestSlice"] = sliceindexlist[stdlist.index(min(stdlist))]

        for keys in superlatives.keys():
            print(keys, superlatives[keys])

        choice = input("Choose a slice by number, please:")
        print()
        self.currentslice = sliceindexlist[int(choice)]
        return sliceindexlist[int(choice)]

    def slice_manager(self):
        """
        Helps the user to determine what 2 dimensional slice to make.

        Arguments:

        ncfile------>NetCDF file to read from

        Returns:

        If everything works, a list describing the slice as a list containing
        both integers and lists, with lists
        indicating what dimensions the slices contain and the integers
        representing where the cross-sections are made.
        """
        print(
            "You have chosen to make a 2d plot in ",
            len(self.ncfile.dimensions),
            " space.",
        )
        print(
            "You must summarily choose a 2d slice for the purposes of 2d\
 plotting."
        )
        print("Your scan dimensions are:")
        print(self.scanvars)
        print("Choose which two are your dimensions of choice.")
        xaxis = input("X:")
        print()
        yaxis = input("Y:")
        print()

        while (
            xaxis not in self.ncfile.dimensions or yaxis not in self.ncfile.dimensions
        ):
            print("Incorrect x and y axes: Please try again.")
            print(self.scanvars)
            print("Choose which two are your dimensions of choice.")
            xaxis = input("X:")
            print()

            yaxis = input("Y:")
            print()

        xvarnum = int(self.scanvars.index(xaxis))
        yvarnum = int(self.scanvars.index(yaxis))

        print("What is the third variable that you are interested in?")
        for var in self.ncfile.variables:
            print(var)

        zaxis = input("Z:")
        print()

        print("Generating information about slices...")

        return self.slicescanner(xvarnum, yvarnum, zaxis)

    def main(self):
        """
        The main loop which manages the main menu, and from which the
        other functions are called.

        Calls:

        slice_manager
        make_scatter_plot
        make_contour_plot

        Dependencies:
        Dataset module
        self.currentslice

        """

        print("``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``")
        print("``'-.,_,.-'`VISUALIZATION UTILITY`'-.,_,.='``")
        print("``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``")

        masterconfigs = {"xname": "", "yname": "", "zname": ""}

        while 1 == 1:

            print("Hello, welcome to the main menu!!")
            print("1) Choose a plotting routine.")
            print("2) Use the slice manager (for n>2 dimensional scans)")
            print("3) Exit.")
            userinput = input(
                "What would you like to do? Please enter a\
 number."
            )
            print()

            if userinput == "1":
                if len(self.ncfile.dimensions) > 2 and self.currentslice is None:
                    print("Warning-- You have not gone to the slice manager yet.")
                    print(
                        "Suggest that you run the slice manager first then\
 return to 2d plots"
                    )
                print("You have the following options.")
                print("2-dimensional scan, or higher dimensions?")
                print("1) 2d.")
                print("2) 3d or more.")
                scaninput = input("Pick one:")
                print()

                if scaninput == "1":
                    print("Choose your x and y axes:")
                    print(self.scanvars)
                    xaxis = input("x:")
                    print()
                    while xaxis not in self.scanvars:
                        print("Sorry- invalid x axis! Please try again...")
                        xaxis = input("x:")
                        print()
                    masterconfigs["xname"] = xaxis
                    yaxis = input("y:")
                    print()
                    while yaxis not in self.scanvars:
                        print("Sorry- invalid y axis! Please try again...")
                        yaxis = input("y:")
                        print()
                    masterconfigs["yname"] = yaxis
                    print("Please choose a Z:")
                    for varbs in self.ncfile.variables:
                        if varbs[0:4] != "SCAN":
                            print(varbs)
                    zaxis = input("z:")
                    print()
                    while zaxis not in self.ncfile.variables:
                        print("Sorry- invalid z axis! Please try again...")
                        zaxis = input("z:")
                        print()
                    masterconfigs["zname"] = zaxis
                    print("And what plotting routine would you like?")
                    print("1) Contour")
                    print("2) Scatter w/ Z Coloring")

                    plotinput = input("Please choose!:")
                    print()
                    if plotinput == "1":
                        self.make_contour_plot(xaxis, yaxis, zaxis, masterconfigs)

                    elif plotinput == "2":
                        self.make_scatter_plot(xaxis, yaxis, zaxis, masterconfigs)

                elif scaninput == "2":
                    print("Okay- activating slice manager!")
                    self.currentslice = self.slice_manager()
                    # TODO: I don't think this actually leads to any plots!
            elif userinput == "2":
                print("Okay- activating slice manager!")
                self.currentslice = self.slice_manager()
            else:
                print("Exiting.")
                break
