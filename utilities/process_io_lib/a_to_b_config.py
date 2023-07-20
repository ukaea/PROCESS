import os
from os.path import abspath


class AToBConfig:
    """Class to read in and store configuration parameters for a_to_b

    Parameters:
     wdir --> Working directory in which the .DAT files for process runs
              will be created
     keep_output --> Switch to keep output for intermediate steps
     outdir --> Directory where the .DAT files for each step will be copied
                after completion
     a_filename --> Name of IN.DAT file for A
     b_filename --> Name of IN.DAT file for B
     path_to_process --> Location of process binary
     vary_niter --> Number of vary_iteration_variables to run when a step
                    can't be completed
     nsteps --> Number of steps to travel from A to B. Setting to 0 will
                mean only the IN.DAT file for A is run
     factor --> When vary_iteration_variables is called, variables will
                be randomly varied between current_value * factor and
                current_value / factor
     bound_gap --> If a variable is set as iteration variable in A but not
                   in B, the targets for the upper and lower bounds of that
                   variable will be set as B_value * bound_gap and
                   B_vale / bound_gap

    """

    wdir = abspath("wdir")
    keep_output = True
    outdir = abspath("steps")
    a = "A.DAT"
    b = "B.DAT"
    path_to_process = "/home/PROCESS/develop/process.exe"
    vary_niter = 20
    nsteps = 10
    factor = 1.2
    bound_gap = 1.001

    def __init__(self, configfile="a_to_b.conf"):
        if not os.path.isfile(configfile):
            return
        with open(configfile) as myfile:
            for line in myfile:
                if line[0] == "*" or line.strip() == "":
                    continue
                elif "=" in line:
                    t_lhs, t_rhs = line.split("=")
                    lhs = t_lhs.strip().lower()
                    rhs = t_rhs.strip()
                    if lhs == "wdir":
                        self.wdir = abspath(rhs)
                        continue
                    elif lhs == "keep_output":
                        if rhs.lower() == "true":
                            self.keep_output = True
                            continue
                        elif rhs.lower() == "false":
                            self.keep_output = False
                            continue
                    elif lhs == "outdir":
                        self.outdir = abspath(rhs)
                        continue
                    elif lhs == "a_filename":
                        self.a = abspath(rhs)
                        continue
                    elif lhs == "b_filename":
                        self.b = abspath(rhs)
                        continue
                    elif lhs == "path_to_process":
                        self.path_to_process = abspath(rhs)
                        continue
                    elif lhs == "vary_niter":
                        if int(rhs) < 0:
                            print("vary_niter must be positive")
                            exit()
                        self.vary_niter = int(rhs)
                        continue
                    elif lhs == "nsteps":
                        if int(rhs) < 0:
                            print("nsteps must be positive")
                            exit()
                        self.nsteps = int(rhs)
                        continue
                    elif lhs == "factor":
                        if float(rhs) < 0:
                            print("factor must be positive")
                            exit()
                        self.factor = float(rhs)
                        continue
                    elif lhs == "bound_gap":
                        if float(rhs) < 0:
                            print("bound_gap must be positive")
                            exit()
                        self.bound_gap = float(rhs)
                        continue

                print("Unrecognised line in {}:".format(configfile))
                print(line)
                break
