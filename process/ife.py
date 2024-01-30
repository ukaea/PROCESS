from process.fortran import constants
from process.fortran import ife_module as ife


# currently the ife module is only partially wrapped
# to unblock the wrapping of availability
class IFE:
    """Module containing Inertial Fusion Energy device routines
    author: P J Knight, CCFE, Culham Science Centre
    N/A
    This module contains routines for calculating the
    parameters of an Inertial Fusion Energy power plant.
    AEA FUS 251: A User's Guide to the PROCESS Systems Code

    NOTE: currently the IFE module is only partially wrapped to unblock the wrapping of availability
    """

    def __init__(self, availability, costs) -> None:
        """Initialises the IFE module's variables

        :param availability: a pointer to the availability model, allowing use of availability's variables/methods
        :type availability: process.availability.Availability
        :param costs: a pointer to the costs model, allowing the use of costs' variables/methods
        :type costs: process.costs.Costs
        """

        self.outfile: int = constants.nout
        self.availability = availability
        self.costs = costs

    def run(self, output: bool):
        """Routine to output the physics and engineering information
        relevant to inertial fusion energy power plants
        author: P J Knight, CCFE, Culham Science Centre

        This routine outputs the physics and engineering information
        relevant to inertial fusion energy power plants.
        F/MI/PJK/LOGBOOK12, p.66
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        # write to output file
        if output:
            # Costs
            self.costs.costs(output=True)

            # Plant availability
            self.availability.avail(output=True)

            # IFE physics
            ife.ifephy(self.outfile, 1)

            # Device build
            ife.ifebld(self.outfile, 1)

            # First wall, blanket and shield
            ife.ifefbs(self.outfile, 1)

            # Device structure
            ife.ifestr()

            # Target data
            ife.ifetgt()

            # Primary thermal power
            ife.ifepw1()

            # Vacuum system
            ife.ifevac()

            # Buildings
            ife.ifebdg(self.outfile, 1)

            # AC power requirements
            ife.ifeacp(self.outfile, 1)

            # Secondary thermal power
            ife.ifepw2(self.outfile, 1)

            return

        # Device build
        ife.ifebld(constants.nout, 0)

        # IFE physics
        ife.ifephy(constants.nout, 0)

        # Device structure
        ife.ifestr()

        # Target data
        ife.ifetgt()

        # First wall, blanket and shield
        ife.ifefbs(constants.nout, 0)

        # Primary thermal power
        ife.ifepw1()

        # Vacuum system
        ife.ifevac()

        # Buildings
        ife.ifebdg(constants.nout, 0)

        # AC power requirements
        ife.ifeacp(constants.nout, 0)

        # Secondary thermal power
        ife.ifepw2(constants.nout, 0)

        # Plant availability
        # TODO: should availability.run be called
        # rather than availability.avail?
        self.availability.avail(output=False)

        # Costs
        self.costs.costs(output=False)
