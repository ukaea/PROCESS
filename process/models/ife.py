"""Module containing Inertial Fusion Energy device routines


This module contains routines for calculating the
parameters of an Inertial Fusion Energy power plant.
"""

import numpy as np

from process.core import constants, process_output
from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.data_structure import physics_variables
from process.data_structure.ife_variables import MAXMAT

MATERIALS = [
    "void",
    "steel",
    "carbon",
    "FLiBe",
    "Li2O",
    "concrete",
    "helium",
    "zenon",
    "lithium",
]


class IFE(Model):
    """Module containing Inertial Fusion Energy device routines

    N/A
    This module contains routines for calculating the
    parameters of an Inertial Fusion Energy power plant.
    """

    def __init__(self, availability, costs):
        """Initialises the IFE module's variables

        Parameters
        ----------
        availability :
            A pointer to the availability model, allowing use of availability's
            variables/methods
        costs :
            A pointer to the costs model, allowing the use of costs'
            variables/methods
        """
        self.outfile: int = constants.NOUT
        self.availability = availability
        self.costs = costs

    def output(self):
        self.run(output=True)

    def run(self, output: bool = False):
        """Routine to output the physics and engineering information
        relevant to inertial fusion energy power plants


        This routine outputs the physics and engineering information
        relevant to inertial fusion energy power plants.
        F/MI/PJK/LOGBOOK12, p.66

        Parameters
        ----------
        output:
            indicate whether output should be written to the output file, or not

        """
        # Device build
        self.ifebld(output=output)

        # IFE physics
        self.ifephy(output=output)

        # Device structure
        self.ifestr()

        # Target data
        self.ifetgt()

        # First wall, blanket and shield
        self.ifefbs(output=output)

        # Primary thermal power
        self.ifepw1()

        # Vacuum system
        self.ifevac()

        # Buildings
        self.ifebdg(output=output)

        # AC power requirements
        self.ifeacp(output=output)

        # Secondary thermal power
        self.ifepw2(output=output)

        # Plant availability
        # TODO: should availability.run be called
        # rather than availability.avail?
        self.availability.avail(output=output)

        # Costs
        self.costs.run()

        if output:
            self.costs.output()

    def ifebld(self, output: bool = False):
        """Routine to create the build of an inertial fusion energy device
        and to calculate the material volumes for the device core


        This routine constructs the build of an inertial fusion energy device
        and calculates the material volumes for the device core.

        Parameters
        ----------
        output:
            boolean to control writing of output to outfile/mfile
        """
        match self.data.ife.ifetyp:
            case 1:
                self.osibld()
            case 2:
                self.sombld()
            case 3:
                self.hylbld()
            case 4:
                self.bld2019()
            case _:
                self.genbld()

        if not output:
            return

        radial_build_data = [
            ["Device centreline", None, 0.0, 0.0],
            ["Chamber", "chrad", self.data.ife.chrad, self.data.ife.r1],
            ["First Wall", "fwdr", self.data.ife.fwdr, self.data.ife.r2],
            ["Void 1", "v1dr", self.data.ife.v1dr, self.data.ife.r3],
            ["Blanket", "bldr", self.data.ife.bldr, self.data.ife.r4],
            ["Void 2", "v2dr", self.data.ife.v2dr, self.data.ife.r5],
            ["Shield", "shdr", self.data.ife.shdr, self.data.ife.r6],
            ["Void 3", "v3dr", self.data.ife.v3dr, self.data.ife.r7],
        ]

        process_output.oheadr(self.outfile, "Radial build")
        process_output.write(
            self.outfile, "\t" * 20 + "Thickness (m)" + "\t" * 2 + "Height (m)"
        )
        for title, _name, thickness, radius in radial_build_data:
            process_output.obuild(self.outfile, title, thickness, radius)

        for title, name, thickness, _radius in radial_build_data:
            if name is None:
                continue
            process_output.ovarre(self.outfile, f"{title} (m)", f"({name})", thickness)

        vertical_build_data = [
            ["Base of device", None, 0.0, -self.data.ife.zl7],
            ["Void 3 lower", "v3dzl", self.data.ife.v3dzl, -self.data.ife.zl6],
            ["Shield lower", "shdzl", self.data.ife.shdzl, -self.data.ife.zl5],
            ["Void 2 lower", "v2dzl", self.data.ife.v2dzl, -self.data.ife.zl4],
            ["Blanket lower", "bldzl", self.data.ife.bldzl, -self.data.ife.zl3],
            ["Void 1 lower", "v1dzl", self.data.ife.v1dzl, -self.data.ife.zl2],
            ["First wall lower", "fwdzl", self.data.ife.fwdzl, -self.data.ife.zl1],
            ["Chamber lower", "chdzl", self.data.ife.chdzl, 0.0],
            ["Chamber upper", "chdzu", self.data.ife.chdzu, self.data.ife.zu1],
            ["First wall upper", "fwdzu", self.data.ife.fwdzu, self.data.ife.zu2],
            ["Void 1 upper", "v1dzu", self.data.ife.v1dzu, self.data.ife.zu3],
            ["Blanket upper", "bldzu", self.data.ife.bldzu, self.data.ife.zu4],
            ["Void 2 upper", "v2dzu", self.data.ife.v2dzu, self.data.ife.zu5],
            ["Shield upper", "shdzu", self.data.ife.shdzu, self.data.ife.zu6],
            ["Void 3 upper", "v3dzu", self.data.ife.v3dzu, self.data.ife.zu7],
        ]

        if self.data.ife.ifetyp == 4:
            process_output.oheadr(self.outfile, "Vertical build - Midplane")
            process_output.write(
                self.outfile, "\t" * 20 + "Thickness (m)" + "\t" * 3 + "Radius (m)"
            )

            for title, _name, thickness, radius in vertical_build_data[:11]:
                process_output.obuild(self.outfile, title, thickness, radius)

            process_output.obuild(
                self.outfile,
                "Blanket upper",
                self.data.ife.bldzu - self.data.ife.bldzu,
                self.data.ife.zu4 - self.data.ife.bldzu,
            )
            process_output.obuild(
                self.outfile,
                "Void 2 upper",
                self.data.ife.v2dzu,
                self.data.ife.zu5 - self.data.ife.bldzu,
            )
            process_output.obuild(
                self.outfile,
                "Shield upper",
                self.data.ife.shdzu,
                self.data.ife.zu6 - self.data.ife.bldzu,
            )
            process_output.obuild(
                self.outfile,
                "Void 3 upper",
                self.data.ife.v3dzu + self.data.ife.bldzu,
                self.data.ife.zu7,
            )

            process_output.oheadr(self.outfile, "Vertical build - Edge")
            process_output.write(
                self.outfile, "\t" * 20 + "Thickness (m)" + "\t" * 3 + "Height (m)"
            )
            for title, _name, thickness, radius in vertical_build_data:
                process_output.obuild(self.outfile, title, thickness, radius)
        else:
            process_output.oheadr(self.outfile, "Vertical build")
            process_output.write(
                self.outfile, "\t" * 20 + "Thickness (m)" + "\t" * 3 + "Height (m)"
            )
            for title, _name, thickness, radius in vertical_build_data:
                process_output.obuild(self.outfile, title, thickness, radius)
            for title, name, thickness, _radius in vertical_build_data:
                if name is None:
                    continue
                process_output.ovarre(
                    self.outfile, f"{title} (m)", f"({name})", thickness
                )

        process_output.write(
            self.outfile,
            "\t" * 6
            + "Chamber".ljust(10)
            + "1st wall".ljust(10)
            + "Void 1".ljust(10)
            + "Blanket".ljust(10)
            + "Void 2".ljust(10)
            + "Shield".ljust(10)
            + "Void 3".ljust(10),
        )
        material_row_gen = _material_string_generator(
            self.data.ife.chmatv,
            self.data.ife.fwmatv,
            self.data.ife.v1matv,
            self.data.ife.blmatv,
            self.data.ife.v2matv,
            self.data.ife.shmatv,
            self.data.ife.v3matv,
        )

        for i in range(len(MATERIALS)):
            process_output.write(self.outfile, material_row_gen(i))

    def osibld(self):
        """Routine to create the build of an inertial fusion energy
        device, based on the design of the OSIRIS study,
        and to calculate the material volumes for the device core
        This routine constructs the build of an inertial fusion energy
        device, based on the design of the OSIRIS study, and to calculate
        the material volumes for the device core.
        """
        # Careful choice of thicknesses, and assuming that the FLiBe
        # inlet radius is small, allows the generic build calculation
        # to be roughly applicable.
        self.genbld()

        # First wall area: no true first wall at bottom of chamber
        self.data.first_wall.a_fw_total = (
            2.0 * np.pi * self.data.ife.r1 * (self.data.ife.zu1 + self.data.ife.zl1)
            + np.pi * self.data.ife.r1 * self.data.ife.r1
        )

    def sombld(self):
        """Routine to create the build of an inertial fusion energy
        device, based on the design of the SOMBRERO study,
        and to calculate the material volumes for the device core
        This routine constructs the build of an inertial fusion energy
        device, based on the design of the SOMBRERO study, and to calculate
        the material volumes for the device core.
        Sviatoslavsky et al, Fusion Technology vol.21 (1992) 1470
        """
        # Radial build
        self.data.ife.r1 = self.data.ife.chrad
        self.data.ife.r2 = self.data.ife.r1 + self.data.ife.fwdr
        self.data.ife.r3 = self.data.ife.r2 + self.data.ife.v1dr
        self.data.ife.r4 = self.data.ife.r3 + self.data.ife.bldr
        self.data.ife.r5 = self.data.ife.r4 + self.data.ife.v2dr
        self.data.ife.r6 = self.data.ife.r5 + self.data.ife.shdr
        self.data.ife.r7 = self.data.ife.r6 + self.data.ife.v3dr

        # Vertical build (below midplane)
        self.data.ife.zl1 = self.data.ife.chdzl
        self.data.ife.zl2 = self.data.ife.zl1 + self.data.ife.fwdzl
        self.data.ife.zl3 = self.data.ife.zl2 + self.data.ife.v1dzl
        self.data.ife.zl4 = self.data.ife.zl3 + self.data.ife.bldzl
        self.data.ife.zl5 = self.data.ife.zl4 + self.data.ife.v2dzl
        self.data.ife.zl6 = self.data.ife.zl5 + self.data.ife.shdzl
        self.data.ife.zl7 = self.data.ife.zl6 + self.data.ife.v3dzl

        # Vertical build (above midplane)
        self.data.ife.zu1 = self.data.ife.chdzu
        self.data.ife.zu2 = self.data.ife.zu1 + self.data.ife.fwdzu
        self.data.ife.zu3 = self.data.ife.zu2 + self.data.ife.v1dzu
        self.data.ife.zu4 = self.data.ife.zu3 + self.data.ife.bldzu
        self.data.ife.zu5 = self.data.ife.zu4 + self.data.ife.v2dzu
        self.data.ife.zu6 = self.data.ife.zu5 + self.data.ife.shdzu
        self.data.ife.zu7 = self.data.ife.zu6 + self.data.ife.v3dzu

        # The SOMBRERO chamber is made up of a cylindrical first wall/
        # blanket, with conical regions above and below. Outside this is
        # a cylindrical shield.

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part

        # Chamber : CHCYLH is the height of the cylindrical part
        chcylh = self.data.ife.chdzu + self.data.ife.chdzl - 2.0 * self.data.ife.chrad

        self.data.ife.chvol = (
            np.pi
            * self.data.ife.r1
            * self.data.ife.r1
            * (chcylh + (2.0 / 3.0) * self.data.ife.chrad)
        )

        # First wall
        self.data.ife.fwvol[0] = (
            np.pi
            * (self.data.ife.r2 * self.data.ife.r2 - self.data.ife.r1 * self.data.ife.r1)
            * chcylh
        )
        self.data.ife.fwvol[1] = (
            (1.0 / 3.0)
            * np.pi
            * (
                self.data.ife.r2
                * self.data.ife.r2
                * (self.data.ife.chrad + self.data.ife.fwdzu)
                - self.data.ife.r1 * self.data.ife.r1 * self.data.ife.chrad
            )
        )
        self.data.ife.fwvol[2] = (
            (1.0 / 3.0)
            * np.pi
            * (
                self.data.ife.r2
                * self.data.ife.r2
                * (self.data.ife.chrad + self.data.ife.fwdzl)
                - self.data.ife.r1 * self.data.ife.r1 * self.data.ife.chrad
            )
        )

        # First void
        self.data.ife.v1vol[0] = (
            np.pi
            * (self.data.ife.r3 * self.data.ife.r3 - self.data.ife.r2 * self.data.ife.r2)
            * chcylh
        )
        self.data.ife.v1vol[1] = (
            (1.0 / 3.0)
            * np.pi
            * (
                self.data.ife.r3
                * self.data.ife.r3
                * (self.data.ife.chrad + self.data.ife.fwdzu + self.data.ife.v1dzu)
                - self.data.ife.r2
                * self.data.ife.r2
                * (self.data.ife.chrad + self.data.ife.fwdzu)
            )
        )
        self.data.ife.v1vol[2] = (
            (1.0 / 3.0)
            * np.pi
            * (
                self.data.ife.r3
                * self.data.ife.r3
                * (self.data.ife.chrad + self.data.ife.fwdzl + self.data.ife.v1dzl)
                - self.data.ife.r2
                * self.data.ife.r2
                * (self.data.ife.chrad + self.data.ife.fwdzl)
            )
        )

        # Blanket:  SOMTDR and SOMBDR are the radii of the cylindrical
        # sections at the top/bottom of the blanket
        # DDZ = Height of top cylindrical section (by similar triangles)
        # DVOL = Volume of top cylindrical section, less the internal cone

        self.data.ife.blvol[0] = (
            np.pi
            * (self.data.ife.r4 * self.data.ife.r4 - self.data.ife.r3 * self.data.ife.r3)
            * chcylh
        )

        self.data.ife.blvol[1] = (
            (1.0 / 3.0)
            * np.pi
            * (
                self.data.ife.r4
                * self.data.ife.r4
                * (
                    self.data.ife.chrad
                    + self.data.ife.fwdzu
                    + self.data.ife.v1dzu
                    + self.data.ife.bldzu
                )
                - self.data.ife.r3
                * self.data.ife.r3
                * (self.data.ife.chrad + self.data.ife.fwdzu + self.data.ife.v1dzu)
            )
        )
        ddz = (
            (
                self.data.ife.chrad
                + self.data.ife.fwdzu
                + self.data.ife.v1dzu
                + self.data.ife.bldzu
            )
            / (
                self.data.ife.chrad
                + self.data.ife.fwdr
                + self.data.ife.v1dr
                + self.data.ife.bldr
            )
            * self.data.ife.somtdr
        )
        dvol = (
            2.0 * (1.0 / 3.0) * np.pi * self.data.ife.somtdr * self.data.ife.somtdr * ddz
        )

        self.data.ife.blvol[1] += dvol

        # Ditto for bottom region...

        self.data.ife.blvol[2] = (
            (1.0 / 3.0)
            * np.pi
            * (
                self.data.ife.r4
                * self.data.ife.r4
                * (
                    self.data.ife.chrad
                    + self.data.ife.fwdzl
                    + self.data.ife.v1dzl
                    + self.data.ife.bldzl
                )
                - self.data.ife.r3
                * self.data.ife.r3
                * (self.data.ife.chrad + self.data.ife.fwdzl + self.data.ife.v1dzl)
            )
        )
        ddz = (
            (
                self.data.ife.chrad
                + self.data.ife.fwdzl
                + self.data.ife.v1dzl
                + self.data.ife.bldzl
            )
            / (
                self.data.ife.chrad
                + self.data.ife.fwdr
                + self.data.ife.v1dr
                + self.data.ife.bldr
            )
            * self.data.ife.sombdr
        )
        dvol = (
            2.0 * (1.0 / 3.0) * np.pi * self.data.ife.sombdr * self.data.ife.sombdr * ddz
        )

        self.data.ife.blvol[2] += dvol

        # Second void
        self.data.ife.v2vol[0] = (
            np.pi
            * (self.data.ife.r5 * self.data.ife.r5 - self.data.ife.r4 * self.data.ife.r4)
            * chcylh
        )
        self.data.ife.v2vol[1] = np.pi * self.data.ife.r5 * self.data.ife.r5 * (
            self.data.ife.zu5 - self.data.ife.chdzu + self.data.ife.chrad
        ) - (
            self.data.ife.fwvol[1]
            + self.data.ife.v1vol[1]
            + self.data.ife.blvol[1]
            + (
                (1.0 / 3.0)
                * np.pi
                * self.data.ife.r1
                * self.data.ife.r1
                * self.data.ife.chrad
            )
        )
        self.data.ife.v2vol[2] = np.pi * self.data.ife.r5 * self.data.ife.r5 * (
            self.data.ife.zl5 - self.data.ife.chdzl + self.data.ife.chrad
        ) - (
            self.data.ife.fwvol[2]
            + self.data.ife.v1vol[2]
            + self.data.ife.blvol[2]
            + (
                (1.0 / 3.0)
                * np.pi
                * self.data.ife.r1
                * self.data.ife.r1
                * self.data.ife.chrad
            )
        )

        # Shield
        self.data.ife.shvol[0] = (
            np.pi
            * (self.data.ife.r6 * self.data.ife.r6 - self.data.ife.r5 * self.data.ife.r5)
            * (self.data.ife.zu6 + self.data.ife.zl6)
        )
        self.data.ife.shvol[1] = (
            np.pi
            * self.data.ife.r5
            * self.data.ife.r5
            * (self.data.ife.zu6 - self.data.ife.zu5)
        )
        self.data.ife.shvol[2] = (
            np.pi
            * self.data.ife.r5
            * self.data.ife.r5
            * (self.data.ife.zl6 - self.data.ife.zl5)
        )

        # Third void
        self.data.ife.v3vol[0] = (
            np.pi
            * (self.data.ife.r7 * self.data.ife.r7 - self.data.ife.r6 * self.data.ife.r6)
            * (self.data.ife.zu7 + self.data.ife.zl7)
        )
        self.data.ife.v3vol[1] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zu7 - self.data.ife.zu6)
        )
        self.data.ife.v3vol[2] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zl7 - self.data.ife.zl6)
        )

        # Material volumes
        for i in range(MAXMAT):
            self.data.ife.chmatv[i] = max(
                0.0, self.data.ife.chvol * self.data.ife.chmatf[i]
            )
            for j in range(3):
                self.data.ife.fwmatv[j, i] = max(
                    0.0, self.data.ife.fwvol[j] * self.data.ife.fwmatf[j, i]
                )
                self.data.ife.v1matv[j, i] = max(
                    0.0, self.data.ife.v1vol[j] * self.data.ife.v1matf[j, i]
                )
                self.data.ife.blmatv[j, i] = max(
                    0.0, self.data.ife.blvol[j] * self.data.ife.blmatf[j, i]
                )
                self.data.ife.v2matv[j, i] = max(
                    0.0, self.data.ife.v2vol[j] * self.data.ife.v2matf[j, i]
                )
                self.data.ife.shmatv[j, i] = max(
                    0.0, self.data.ife.shvol[j] * self.data.ife.shmatf[j, i]
                )
                self.data.ife.v3matv[j, i] = max(
                    0.0, self.data.ife.v3vol[j] * self.data.ife.v3matf[j, i]
                )

        # First wall area
        self.data.first_wall.a_fw_total = (
            2.0
            * np.pi
            * self.data.ife.r1
            * ((self.data.ife.zu1 + self.data.ife.zl1) + self.data.ife.r1 * np.sqrt(2.0))
        )

    def hylbld(self):
        # Radial build
        self.data.ife.r1 = self.data.ife.chrad
        self.data.ife.r2 = self.data.ife.r1 + self.data.ife.fwdr
        self.data.ife.r3 = self.data.ife.r2 + self.data.ife.v1dr
        self.data.ife.r4 = self.data.ife.r3 + self.data.ife.bldr
        self.data.ife.r5 = self.data.ife.r4 + self.data.ife.v2dr
        self.data.ife.r6 = self.data.ife.r5 + self.data.ife.shdr
        self.data.ife.r7 = self.data.ife.r6 + self.data.ife.v3dr

        # Vertical build (below midplane)
        self.data.ife.zl1 = self.data.ife.chdzl
        self.data.ife.zl2 = self.data.ife.zl1 + self.data.ife.fwdzl
        self.data.ife.zl3 = self.data.ife.zl2 + self.data.ife.v1dzl
        self.data.ife.zl4 = self.data.ife.zl3 + self.data.ife.bldzl
        self.data.ife.zl5 = self.data.ife.zl4 + self.data.ife.v2dzl
        self.data.ife.zl6 = self.data.ife.zl5 + self.data.ife.shdzl
        self.data.ife.zl7 = self.data.ife.zl6 + self.data.ife.v3dzl

        # Vertical build (above midplane)
        self.data.ife.zu1 = self.data.ife.chdzu
        self.data.ife.zu2 = self.data.ife.zu1 + self.data.ife.fwdzu
        self.data.ife.zu3 = self.data.ife.zu2 + self.data.ife.v1dzu
        self.data.ife.zu4 = self.data.ife.zu3 + self.data.ife.bldzu
        self.data.ife.zu5 = self.data.ife.zu4 + self.data.ife.v2dzu
        self.data.ife.zu6 = self.data.ife.zu5 + self.data.ife.shdzu
        self.data.ife.zu7 = self.data.ife.zu6 + self.data.ife.v3dzu

        # The HYLIFE-II chamber is assumed to be mostly cylindrical, but
        # with a conical region below the midplane that causes the Flibe
        # to flow downwards and outwards towards the outlet.

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part
        self.data.ife.chvol = (
            np.pi
            * self.data.ife.r1
            * self.data.ife.r1
            * (
                (self.data.ife.zu1 + self.data.ife.zl5)
                - (1 / 3) * (self.data.ife.zl5 - self.data.ife.zl1)
            )
        )

        # First wall
        # FLIRAD is the radius of the Flibe inlet
        self.data.ife.fwvol[0] = (
            np.pi
            * (self.data.ife.r2 * self.data.ife.r2 - self.data.ife.r1 * self.data.ife.r1)
            * (self.data.ife.zu2 + self.data.ife.zl5)
        )
        self.data.ife.fwvol[1] = (
            np.pi
            * (
                self.data.ife.r1 * self.data.ife.r1
                - self.data.ife.flirad * self.data.ife.flirad
            )
            * (self.data.ife.zu2 - self.data.ife.zu1)
        )
        self.data.ife.fwvol[2] = (
            (1 / 3)
            * np.pi
            * (
                self.data.ife.r2
                * self.data.ife.r2
                * (self.data.ife.zl5 - self.data.ife.zl1)
                - self.data.ife.r1
                * self.data.ife.r1
                * (self.data.ife.zl5 - self.data.ife.zl2)
            )
        )

        # First void
        self.data.ife.v1vol[0] = (
            np.pi
            * (self.data.ife.r3 * self.data.ife.r3 - self.data.ife.r2 * self.data.ife.r2)
            * (self.data.ife.zu2 + self.data.ife.zl3)
        )
        self.data.ife.v1vol[1] = (
            np.pi
            * (
                self.data.ife.r4 * self.data.ife.r4
                - self.data.ife.flirad * self.data.ife.flirad
            )
            * (self.data.ife.zu3 - self.data.ife.zu2)
        )
        self.data.ife.v1vol[2] = (
            (1 / 3)
            * np.pi
            * self.data.ife.r1
            * self.data.ife.r1
            * (self.data.ife.zl3 - self.data.ife.zl2)
        )

        # Blanket
        self.data.ife.blvol[0] = (
            np.pi
            * (self.data.ife.r4 * self.data.ife.r4 - self.data.ife.r3 * self.data.ife.r3)
            * (self.data.ife.zu2 + self.data.ife.zl3)
        )
        self.data.ife.blvol[1] = (
            np.pi
            * (
                self.data.ife.r4 * self.data.ife.r4
                - self.data.ife.flirad * self.data.ife.flirad
            )
            * (self.data.ife.zu4 - self.data.ife.zu3)
        )
        self.data.ife.blvol[2] = (
            np.pi
            * self.data.ife.r4
            * self.data.ife.r4
            * (self.data.ife.zl4 - self.data.ife.zl3)
        )

        # Second void
        self.data.ife.v2vol[0] = (
            np.pi
            * (self.data.ife.r5 * self.data.ife.r5 - self.data.ife.r4 * self.data.ife.r4)
            * (self.data.ife.zu4 + self.data.ife.zl4)
        )
        self.data.ife.v2vol[1] = (
            np.pi
            * (
                self.data.ife.r5 * self.data.ife.r5
                - self.data.ife.flirad * self.data.ife.flirad
            )
            * (self.data.ife.zu5 - self.data.ife.zu4)
        )
        self.data.ife.v2vol[2] = (
            np.pi
            * self.data.ife.r5
            * self.data.ife.r5
            * (self.data.ife.zl5 - self.data.ife.zl4)
        )

        # Shield
        self.data.ife.shvol[0] = (
            np.pi
            * (self.data.ife.r6 * self.data.ife.r6 - self.data.ife.r5 * self.data.ife.r5)
            * (self.data.ife.zu5 + self.data.ife.zl5)
        )
        self.data.ife.shvol[1] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zu6 - self.data.ife.zu5)
        )
        self.data.ife.shvol[2] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zl6 - self.data.ife.zl5)
        )

        # Third void
        self.data.ife.v3vol[0] = (
            np.pi
            * (self.data.ife.r7 * self.data.ife.r7 - self.data.ife.r6 * self.data.ife.r6)
            * (self.data.ife.zu6 + self.data.ife.zl6)
        )
        self.data.ife.v3vol[1] = (
            np.pi
            * self.data.ife.r7
            * self.data.ife.r7
            * (self.data.ife.zu7 - self.data.ife.zu6)
        )
        self.data.ife.v3vol[2] = (
            np.pi
            * self.data.ife.r7
            * self.data.ife.r7
            * (self.data.ife.zl7 - self.data.ife.zl6)
        )

        # Material volumes
        for i in range(MAXMAT + 1):
            self.data.ife.chmatv[i] = max(
                0.0, self.data.ife.chvol * self.data.ife.chmatf[i]
            )
            for j in range(3):
                self.data.ife.fwmatv[j, i] = max(
                    0.0, self.data.ife.fwvol[j] * self.data.ife.fwmatf[j, i]
                )
                self.data.ife.v1matv[j, i] = max(
                    0.0, self.data.ife.v1vol[j] * self.data.ife.v1matf[j, i]
                )
                self.data.ife.blmatv[j, i] = max(
                    0.0, self.data.ife.blvol[j] * self.data.ife.blmatf[j, i]
                )
                self.data.ife.v2matv[j, i] = max(
                    0.0, self.data.ife.v2vol[j] * self.data.ife.v2matf[j, i]
                )
                self.data.ife.shmatv[j, i] = max(
                    0.0, self.data.ife.shvol[j] * self.data.ife.shmatf[j, i]
                )
                self.data.ife.v3matv[j, i] = max(
                    0.0, self.data.ife.v3vol[j] * self.data.ife.v3matf[j, i]
                )

        # First wall area
        self.data.first_wall.a_fw_total = (
            2.0 * np.pi * self.data.ife.r1 * (self.data.ife.zu1 + self.data.ife.zl5)
        )
        self.data.first_wall.a_fw_total += np.pi * (
            self.data.ife.r1 * self.data.ife.r1
            - self.data.ife.flirad * self.data.ife.flirad
        )
        self.data.first_wall.a_fw_total += (
            np.pi
            * self.data.ife.r1
            * np.sqrt(
                self.data.ife.r1 * self.data.ife.r1
                + (self.data.ife.zl3 - self.data.ife.zl1) ** 2
            )
        )

    def bld2019(self):
        """Routine to create the build of a 2019 inertial fusion energy
        device, and to calculate the material volumes for the device core
        This routine constructs the build of a modern inertial fusion energy
        device, assumed to be cylindrically-symmetric, with a pool at bottom
        and top corners and with a lower shield at the centre.  See diagram
        attached to Issue #907.
        Issue #907
        """
        # Check input
        if self.data.ife.fwdr > 0 or self.data.ife.v1dr > 0:
            raise ProcessValueError("fwdr and v1dr should be zero for 2019 IFE build")
        if self.data.ife.fwdzu > 0 or self.data.ife.v1dzu > 0 or self.data.ife.v2dzu > 0:
            raise ProcessValueError(
                "fwdzu, v1dzu and v2dzu should be zero for 2019 IFE build"
            )
        if self.data.ife.fwdzl > 0 or self.data.ife.v1dzl > 0 or self.data.ife.v2dzu > 0:
            raise ProcessValueError(
                "fwdzl, v1dzl and v2dzl should be zero for 2019 IFE build"
            )

        # Lithium Pump
        # Velocity
        vel = np.sqrt(
            2.0
            * constants.ACCELERATION_GRAVITY
            * (self.data.ife.chdzu + self.data.ife.bldzu)
        )

        # Lithium Fraction
        self.data.ife.blmatf[0, 8] = 0.91 * np.sqrt(
            self.data.ife.bldzu / (self.data.ife.chdzu + self.data.ife.bldzu)
        )
        self.data.ife.blmatf[0, 0] = 1.0 - self.data.ife.blmatf[0, 8]

        # Spatial Thickness
        self.data.ife.bldr = self.data.ife.bldrc / self.data.ife.blmatf[0, 8]

        # Area
        acurt = np.pi * (
            (self.data.ife.chrad + self.data.ife.bldr) ** 2.0 - self.data.ife.chrad**2.0
        )

        # Mass Flow
        mdot = 512.0 * vel * self.data.ife.blmatf[0, 8] * acurt

        # Pump Power (MW)
        self.data.ife.lipmw = (
            1e-6
            * mdot
            * constants.ACCELERATION_GRAVITY
            * (
                self.data.ife.chdzl
                + self.data.ife.chdzu
                + self.data.ife.bldzu
                + self.data.ife.bldzl
            )
            / self.data.ife.etali
        )

        # Fall Time
        self.data.ife.taufall = (
            2.0 * (self.data.ife.chdzl + self.data.ife.chdzu + self.data.ife.bldzu) / vel
        )

        self.data.ife.rrmax = 1.0 / self.data.ife.taufall

        # TBR and Emult model was for spherical lithium
        # Remove reactor head
        phi = np.arctan(self.data.ife.chrad / self.data.ife.chdzu)
        sang = 1.0 - np.cos(phi)
        li_frac = 1.0 - 0.5 * sang

        # TBR
        self.data.fwbs.tbr = (
            3.7418
            * (1.0 / (1.0 + np.exp(-2.6366 * self.data.ife.bldrc)) - 0.5)
            * li_frac
        )

        # Energy Multiplication
        self.data.fwbs.f_p_blkt_multiplication = (
            2.2414
            * (1.0 / (1.0 + np.exp(-3.0038 * self.data.ife.bldrc)) - 0.5)
            * li_frac
        )

        # Radial build

        self.data.ife.r1 = self.data.ife.chrad
        self.data.ife.r2 = self.data.ife.r1 + self.data.ife.fwdr
        self.data.ife.r3 = self.data.ife.r2 + self.data.ife.v1dr
        self.data.ife.r4 = self.data.ife.r3 + self.data.ife.bldr
        self.data.ife.r5 = self.data.ife.r4 + self.data.ife.v2dr
        self.data.ife.r6 = self.data.ife.r5 + self.data.ife.shdr
        self.data.ife.r7 = self.data.ife.r6 + self.data.ife.v3dr

        # Vertical build (below midplane)

        self.data.ife.zl1 = self.data.ife.chdzl
        self.data.ife.zl2 = self.data.ife.zl1 + self.data.ife.fwdzl
        self.data.ife.zl3 = self.data.ife.zl2 + self.data.ife.v1dzl
        self.data.ife.zl4 = self.data.ife.zl3 + self.data.ife.bldzl
        self.data.ife.zl5 = self.data.ife.zl4 + self.data.ife.v2dzl
        self.data.ife.zl6 = self.data.ife.zl5 + self.data.ife.shdzl
        self.data.ife.zl7 = self.data.ife.zl6 + self.data.ife.v3dzl

        # Vertical build (above midplane)

        self.data.ife.zu1 = self.data.ife.chdzu
        self.data.ife.zu2 = self.data.ife.zu1 + self.data.ife.fwdzu
        self.data.ife.zu3 = self.data.ife.zu2 + self.data.ife.v1dzu
        self.data.ife.zu4 = self.data.ife.zu3 + self.data.ife.bldzu
        self.data.ife.zu5 = self.data.ife.zu4 + self.data.ife.v2dzu
        self.data.ife.zu6 = self.data.ife.zu5 + self.data.ife.shdzu

        self.data.ife.v3dzu = (
            (self.data.ife.zu6 + self.data.ife.zl6)
            + self.data.buildings.trcl
            + self.data.buildings.stcl
            + 5.1
            + 9.41e-6 * 1.0e5
        )

        self.data.ife.zu7 = self.data.ife.zu6 + self.data.ife.v3dzu

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part

        # Chamber

        chvol = (
            np.pi
            * self.data.ife.r1
            * self.data.ife.r1
            * (self.data.ife.zu1 + self.data.ife.zl1)
        )

        # First wall

        self.data.ife.fwvol[0] = (
            np.pi
            * (self.data.ife.r2 * self.data.ife.r2 - self.data.ife.r1 * self.data.ife.r1)
            * (self.data.ife.zu1 + self.data.ife.zl1)
        )
        self.data.ife.fwvol[1] = (
            np.pi
            * self.data.ife.r2
            * self.data.ife.r2
            * (self.data.ife.zu2 - self.data.ife.zu1)
        )
        self.data.ife.fwvol[2] = (
            np.pi
            * self.data.ife.r2
            * self.data.ife.r2
            * (self.data.ife.zl2 - self.data.ife.zl1)
        )

        # First void

        self.data.ife.v1vol[0] = (
            np.pi
            * (self.data.ife.r3 * self.data.ife.r3 - self.data.ife.r2 * self.data.ife.r2)
            * (self.data.ife.zu2 + self.data.ife.zl2)
        )
        self.data.ife.v1vol[1] = (
            np.pi
            * self.data.ife.r3
            * self.data.ife.r3
            * (self.data.ife.zu3 - self.data.ife.zu2)
        )
        self.data.ife.v1vol[2] = (
            np.pi
            * self.data.ife.r3
            * self.data.ife.r3
            * (self.data.ife.zl3 - self.data.ife.zl2)
        )

        # Blanket
        # Radial Blanket - between void 2 and chamber
        self.data.ife.blvol[0] = (
            np.pi
            * (self.data.ife.r4 * self.data.ife.r4 - self.data.ife.r3 * self.data.ife.r3)
            * (self.data.ife.zu3 + self.data.ife.zl3)
        )
        # Upper Blanket - Pool radially between shield and
        # chamber of input height.
        self.data.ife.blvol[1] = (
            np.pi
            * (self.data.ife.r5 * self.data.ife.r5 - self.data.ife.r3 * self.data.ife.r3)
            * self.data.ife.bldzu
        )
        # Lower Blanket - Pool filling base of device
        self.data.ife.blvol[2] = (
            np.pi
            * self.data.ife.r5
            * self.data.ife.r5
            * (self.data.ife.zl4 - self.data.ife.zl3)
        )

        # Second void

        self.data.ife.v2vol[0] = (
            np.pi
            * (self.data.ife.r5 * self.data.ife.r5 - self.data.ife.r4 * self.data.ife.r4)
            * (self.data.ife.chdzl + self.data.ife.chdzu)
        )
        self.data.ife.v2vol[1] = 0.0
        self.data.ife.v2vol[2] = 0.0

        # Shield
        self.data.ife.shvol[0] = (
            np.pi
            * (self.data.ife.r6 * self.data.ife.r6 - self.data.ife.r5 * self.data.ife.r5)
            * (self.data.ife.zu5 + self.data.ife.zl5)
        )
        # Top Section is in three parts to account for the dip at
        # the centre.  The first is the horizontal top, the second is the
        # horizontal
        self.data.ife.shvol[1] = np.pi * (
            (
                (
                    self.data.ife.r6 * self.data.ife.r6
                    - (self.data.ife.chrad - self.data.ife.shdr)
                    * (self.data.ife.chrad - self.data.ife.shdr)
                )
                * self.data.ife.shdzu
            )
            + (
                (
                    self.data.ife.r1 * self.data.ife.r1
                    - self.data.ife.flirad * self.data.ife.flirad
                )
                * self.data.ife.shdzu
            )
            + (
                (
                    self.data.ife.r1 * self.data.ife.r1
                    - (self.data.ife.r1 - self.data.ife.shdzu)
                    * (self.data.ife.r1 - self.data.ife.shdzu)
                )
                * (self.data.ife.bldzu - self.data.ife.shdzu)
            )
        )
        self.data.ife.shvol[2] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zl6 - self.data.ife.zl5)
        )

        # Third void

        self.data.ife.v3vol[0] = (
            np.pi
            * (self.data.ife.r7 * self.data.ife.r7 - self.data.ife.r6 * self.data.ife.r6)
            * (self.data.ife.zu6 + self.data.ife.zl6)
        )
        self.data.ife.v3vol[1] = (
            np.pi
            * self.data.ife.r7
            * self.data.ife.r7
            * (self.data.ife.zu7 - self.data.ife.zu6)
            + np.pi
            * (
                (self.data.ife.r1 - self.data.ife.shdzu)
                * (self.data.ife.r1 - self.data.ife.shdzu)
                - self.data.ife.flirad * self.data.ife.flirad
            )
            * self.data.ife.bldzu
        )
        self.data.ife.v3vol[2] = (
            np.pi
            * self.data.ife.r7
            * self.data.ife.r7
            * (self.data.ife.zl7 - self.data.ife.zl6)
        )

        # Material volumes

        for i in range(MAXMAT + 1):
            self.data.ife.chmatv[i] = max(0.0, chvol * self.data.ife.chmatf[i])
            for j in range(3):
                self.data.ife.fwmatv[j, i] = max(
                    0.0, self.data.ife.fwvol[j] * self.data.ife.fwmatf[j, i]
                )
                self.data.ife.v1matv[j, i] = max(
                    0.0, self.data.ife.v1vol[j] * self.data.ife.v1matf[j, i]
                )
                self.data.ife.blmatv[j, i] = max(
                    0.0, self.data.ife.blvol[j] * self.data.ife.blmatf[j, i]
                )
                self.data.ife.v2matv[j, i] = max(
                    0.0, self.data.ife.v2vol[j] * self.data.ife.v2matf[j, i]
                )
                self.data.ife.shmatv[j, i] = max(
                    0.0, self.data.ife.shvol[j] * self.data.ife.shmatf[j, i]
                )
                self.data.ife.v3matv[j, i] = max(
                    0.0, self.data.ife.v3vol[j] * self.data.ife.v3matf[j, i]
                )

        # First wall area
        # The chamber is surrounded by liquid on three sides
        # with only the top being solid.  This is considered part
        # of the shield. There is a target injector tube at the
        # centre of this area.
        self.data.first_wall.a_fw_total = np.pi * (
            self.data.ife.r1 * self.data.ife.r1
            - self.data.ife.flirad * self.data.ife.flirad
        )

    def genbld(self):
        """Routine to create the build of a generic inertial fusion energy
        device, and to calculate the material volumes for the device core
        This routine constructs the build of a generic inertial fusion energy
        device, assumed to be cylindrically-symmetric, and to calculate
        the material volumes for the device core.
        F/MI/PJK/LOGBOOK12, p.52
        """
        # Radial build

        self.data.ife.r1 = self.data.ife.chrad
        self.data.ife.r2 = self.data.ife.r1 + self.data.ife.fwdr
        self.data.ife.r3 = self.data.ife.r2 + self.data.ife.v1dr
        self.data.ife.r4 = self.data.ife.r3 + self.data.ife.bldr
        self.data.ife.r5 = self.data.ife.r4 + self.data.ife.v2dr
        self.data.ife.r6 = self.data.ife.r5 + self.data.ife.shdr
        self.data.ife.r7 = self.data.ife.r6 + self.data.ife.v3dr

        # Vertical build (below midplane)

        self.data.ife.zl1 = self.data.ife.chdzl
        self.data.ife.zl2 = self.data.ife.zl1 + self.data.ife.fwdzl
        self.data.ife.zl3 = self.data.ife.zl2 + self.data.ife.v1dzl
        self.data.ife.zl4 = self.data.ife.zl3 + self.data.ife.bldzl
        self.data.ife.zl5 = self.data.ife.zl4 + self.data.ife.v2dzl
        self.data.ife.zl6 = self.data.ife.zl5 + self.data.ife.shdzl
        self.data.ife.zl7 = self.data.ife.zl6 + self.data.ife.v3dzl

        # Vertical build (above midplane)

        self.data.ife.zu1 = self.data.ife.chdzu
        self.data.ife.zu2 = self.data.ife.zu1 + self.data.ife.fwdzu
        self.data.ife.zu3 = self.data.ife.zu2 + self.data.ife.v1dzu
        self.data.ife.zu4 = self.data.ife.zu3 + self.data.ife.bldzu
        self.data.ife.zu5 = self.data.ife.zu4 + self.data.ife.v2dzu
        self.data.ife.zu6 = self.data.ife.zu5 + self.data.ife.shdzu
        self.data.ife.zu7 = self.data.ife.zu6 + self.data.ife.v3dzu

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part

        # Chamber

        chvol = (
            np.pi
            * self.data.ife.r1
            * self.data.ife.r1
            * (self.data.ife.zu1 + self.data.ife.zl1)
        )

        # First wall

        self.data.ife.fwvol[0] = (
            np.pi
            * (self.data.ife.r2 * self.data.ife.r2 - self.data.ife.r1 * self.data.ife.r1)
            * (self.data.ife.zu1 + self.data.ife.zl1)
        )
        self.data.ife.fwvol[1] = (
            np.pi
            * self.data.ife.r2
            * self.data.ife.r2
            * (self.data.ife.zu2 - self.data.ife.zu1)
        )
        self.data.ife.fwvol[2] = (
            np.pi
            * self.data.ife.r2
            * self.data.ife.r2
            * (self.data.ife.zl2 - self.data.ife.zl1)
        )

        # First void

        self.data.ife.v1vol[0] = (
            np.pi
            * (self.data.ife.r3 * self.data.ife.r3 - self.data.ife.r2 * self.data.ife.r2)
            * (self.data.ife.zu2 + self.data.ife.zl2)
        )
        self.data.ife.v1vol[1] = (
            np.pi
            * self.data.ife.r3
            * self.data.ife.r3
            * (self.data.ife.zu3 - self.data.ife.zu2)
        )
        self.data.ife.v1vol[2] = (
            np.pi
            * self.data.ife.r3
            * self.data.ife.r3
            * (self.data.ife.zl3 - self.data.ife.zl2)
        )

        # Blanket

        self.data.ife.blvol[0] = (
            np.pi
            * (self.data.ife.r4 * self.data.ife.r4 - self.data.ife.r3 * self.data.ife.r3)
            * (self.data.ife.zu3 + self.data.ife.zl3)
        )
        self.data.ife.blvol[1] = (
            np.pi
            * self.data.ife.r4
            * self.data.ife.r4
            * (self.data.ife.zu4 - self.data.ife.zu3)
        )
        self.data.ife.blvol[2] = (
            np.pi
            * self.data.ife.r4
            * self.data.ife.r4
            * (self.data.ife.zl4 - self.data.ife.zl3)
        )

        # Second void

        self.data.ife.v2vol[0] = (
            np.pi
            * (self.data.ife.r5 * self.data.ife.r5 - self.data.ife.r4 * self.data.ife.r4)
            * (self.data.ife.zu4 + self.data.ife.zl4)
        )
        self.data.ife.v2vol[1] = (
            np.pi
            * self.data.ife.r5
            * self.data.ife.r5
            * (self.data.ife.zu5 - self.data.ife.zu4)
        )
        self.data.ife.v2vol[2] = (
            np.pi
            * self.data.ife.r5
            * self.data.ife.r5
            * (self.data.ife.zl5 - self.data.ife.zl4)
        )

        # Shield

        self.data.ife.shvol[0] = (
            np.pi
            * (self.data.ife.r6 * self.data.ife.r6 - self.data.ife.r5 * self.data.ife.r5)
            * (self.data.ife.zu5 + self.data.ife.zl5)
        )
        self.data.ife.shvol[1] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zu6 - self.data.ife.zu5)
        )
        self.data.ife.shvol[2] = (
            np.pi
            * self.data.ife.r6
            * self.data.ife.r6
            * (self.data.ife.zl6 - self.data.ife.zl5)
        )

        # Third void

        self.data.ife.v3vol[0] = (
            np.pi
            * (self.data.ife.r7 * self.data.ife.r7 - self.data.ife.r6 * self.data.ife.r6)
            * (self.data.ife.zu6 + self.data.ife.zl6)
        )
        self.data.ife.v3vol[1] = (
            np.pi
            * self.data.ife.r7
            * self.data.ife.r7
            * (self.data.ife.zu7 - self.data.ife.zu6)
        )
        self.data.ife.v3vol[2] = (
            np.pi
            * self.data.ife.r7
            * self.data.ife.r7
            * (self.data.ife.zl7 - self.data.ife.zl6)
        )

        # Material volumes

        for i in range(MAXMAT + 1):
            self.data.ife.chmatv[i] = max(0.0, chvol * self.data.ife.chmatf[i])
            for j in range(3):
                self.data.ife.fwmatv[j, i] = max(
                    0.0, self.data.ife.fwvol[j] * self.data.ife.fwmatf[j, i]
                )
                self.data.ife.v1matv[j, i] = max(
                    0.0, self.data.ife.v1vol[j] * self.data.ife.v1matf[j, i]
                )
                self.data.ife.blmatv[j, i] = max(
                    0.0, self.data.ife.blvol[j] * self.data.ife.blmatf[j, i]
                )
                self.data.ife.v2matv[j, i] = max(
                    0.0, self.data.ife.v2vol[j] * self.data.ife.v2matf[j, i]
                )
                self.data.ife.shmatv[j, i] = max(
                    0.0, self.data.ife.shvol[j] * self.data.ife.shmatf[j, i]
                )
                self.data.ife.v3matv[j, i] = max(
                    0.0, self.data.ife.v3vol[j] * self.data.ife.v3matf[j, i]
                )

        # First wall area

        self.data.first_wall.a_fw_total = (
            2.0
            * np.pi
            * self.data.ife.r1
            * ((self.data.ife.zu1 + self.data.ife.zl1) + self.data.ife.r1)
        )

    def ifephy(self, output: bool = False):
        """Routine to calculate the physics parameters of an Inertial Fusion
        Energy power plant


        This routine calculates the physics parameters of an Inertial Fusion
        Energy power plant.
        F/MI/PJK/LOGBOOK12, pp.68,85

        Parameters
        ----------
        output: bool
             (Default value = False)
        """
        match self.data.ife.ifedrv:
            case -1:
                # Target gain and driver efficiency dependencies on
                # driver energy are input
                self.data.ife.gain, self.data.ife.etadrv = self.driver(
                    self.data.ife.edrive, self.data.ife.gainve, self.data.ife.etave
                )
            case 0:  # Target gain and driver efficiency are input
                self.data.ife.gain = self.data.ife.tgain
                self.data.ife.etadrv = self.data.ife.drveff
            case 1:  # Laser driver based on SOMBRERO design
                self.data.ife.gain, self.data.ife.etadrv = self.lasdrv(
                    self.data.ife.edrive
                )
            case 2:  # Heavy-ion beam driver based on OSIRIS design
                self.data.ife.gain, self.data.ife.etadrv = self.iondrv(
                    self.data.ife.edrive
                )
            case 3:
                self.data.ife.etadrv = self.data.ife.drveff
            case _:
                raise ProcessValueError(
                    f"ifedrv={self.data.ife.ifedrv} is an invalid option"
                )

        if self.data.ife.ifedrv != 3:
            # Repetition rate (Hz)
            self.data.ife.reprat = self.data.ife.pdrive / self.data.ife.edrive
            # Fusion power (MW)
            physics_variables.p_fusion_total_mw = (
                1.0e-6 * self.data.ife.pdrive * self.data.ife.gain
            )
        else:
            # Driver Power
            self.data.ife.reprat = self.data.ife.rrin
            self.data.ife.pdrive = self.data.ife.reprat * self.data.ife.edrive
            # Gain
            physics_variables.p_fusion_total_mw = self.data.ife.pfusife
            self.data.ife.gain = physics_variables.p_fusion_total_mw / (
                1.0e-6 * self.data.ife.pdrive
            )

        # Wall load (assume total fusion power applies)

        if self.data.ife.ifetyp == 1:
            # OSIRIS-type build: First wall subtends a solid angle of 2 pi * SANG

            phi = 0.5 * np.pi + np.arctan(self.data.ife.zl1 / self.data.ife.r1)
            sang = 1.0 - np.cos(phi)
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_fusion_total_mw
                * 0.5
                * sang
                / self.data.first_wall.a_fw_total
            )

        elif self.data.ife.ifetyp == 4:
            # 2019 build only has first wall at the top which has a tube at
            # its centre.  This calculates solid angle and removes tube.

            phi = np.arctan(self.data.ife.r1 / self.data.ife.zu1)
            sang = 1.0 - np.cos(phi)
            phi = np.arctan(self.data.ife.flirad / self.data.ife.zu1)
            sang -= 1.0 - np.cos(phi)
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_fusion_total_mw
                * 0.5
                * sang
                / self.data.first_wall.a_fw_total
            )

        else:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_fusion_total_mw / self.data.first_wall.a_fw_total
            )

        if not output:
            return

        process_output.oheadr(self.outfile, "Physics / Driver Issues")

        match self.data.ife.ifedrv:
            case -1 | 0:
                process_output.ocmmnt(self.outfile, "Driver type : generic")
            case 1:
                process_output.ocmmnt(self.outfile, "Driver type : laser")
            case 2:
                process_output.ocmmnt(self.outfile, "Driver type : heavy ion beam")

        process_output.oblnkl(self.outfile)

        process_output.ovarre(
            self.outfile, "Driver energy (J)", "(edrive)", self.data.ife.edrive
        )
        process_output.ovarre(
            self.outfile, "Driver efficiency", "(etadrv)", self.data.ife.etadrv
        )
        process_output.ovarre(
            self.outfile,
            "Driver power reaching target (W)",
            "(pdrive)",
            self.data.ife.pdrive,
        )
        process_output.ovarre(
            self.outfile,
            "Driver repetition rate (Hz)",
            "(reprat)",
            self.data.ife.reprat,
        )
        process_output.ovarre(self.outfile, "Target gain", "(gain)", self.data.ife.gain)
        process_output.ovarre(
            self.outfile,
            "Fusion power (MW)",
            "(p_fusion_total_mw)",
            physics_variables.p_fusion_total_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Neutron wall load (MW/m2)",
            "(pflux_fw_neutron_mw)",
            physics_variables.pflux_fw_neutron_mw,
        )

    def driver(self, edrive, gainve, etave):
        """Routine to calculate parameters of a generic driver
        suitable for inertial fusion energy

        Parameters
        ----------
        edrive:
            Driver energy (J)
        gainve(10):
            Gain vs energy data
        etave(10):
            Driver efficiency vs energy data

        Returns
        -------
        gain   :
            Target gain
        etadrv :
            Driver efficiency

        Notes
        -----
        This routine calculates the parameters of a generic driver
        suitable for inertial fusion energy.
        Gain and driver efficiency data are interpolated from input data.
        F/MI/PJK/LOGBOOK12, p.85
        """
        # The arrays contain data points for EDRIVE = 1MJ, 2MJ, ... , 10MJ

        e = 1e-6 * edrive
        ie = int(e)
        de = e - ie

        # Assume linear interpolations and extrapolations

        if ie <= 1:
            gain = gainve[1] - 1.0e-6 * (edrive - 2.0e6) * (gainve[0] - gainve[1])
            etadrv = etave[1] - 1.0e-6 * (edrive - 2.0e6) * (etave[0] - etave[1])
        elif ie >= 9:
            gain = gainve[8] + 1.0e-6 * (edrive - 9.0e6) * (gainve[9] - gainve[8])
            etadrv = etave[8] + 1.0e-6 * (edrive - 9.0e6) * (etave[9] - etave[8])
        else:
            gain = gainve[ie - 1] + de * (gainve[ie] - gainve[ie - 1])
            etadrv = etave[ie - 1] + de * (etave[ie] - etave[ie - 1])

        # Ensure sensible values

        gain = max(0.01, gain)
        etadrv = max(0.01, etadrv)

        return gain, etadrv

    def lasdrv(self, edrive):
        """Routine to calculate parameters of a laser driver
        suitable for inertial fusion energy

        Parameters
        ----------
        edrive:
            Driver energy (J)

        Returns
        -------
        gain   : output
            Target gain
        etadrv : output
            Driver efficiency

        Notes
        -----
        This routine calculates the parameters of a laser driver
        suitable for inertial fusion energy.
        Gain and driver efficiency data are taken from Figures 1 and 2 of
        Meier and Rosenberg.
        Meier and Rosenberg, Fusion Technology vol.21 (1992) p.1552
        F/MI/PJK/LOGBOOK12, p.86
        """
        # GVE(K): target gain at EDRIVE = K MegaJoules
        gve = [63.0, 95.0, 112.0, 125.0, 136.0, 144.0, 151.0, 157.0, 162.0, 166.0]

        # EVE(K): driver efficiency at EDRIVE = K MegaJoules
        eve = [0.082, 0.079, 0.076, 0.072, 0.069, 0.064, 0.059, 0.054, 0.048, 0.042]

        e = 1.0e-6 * edrive
        ie = int(e)
        de = e - ie

        # Assume linear interpolations and extrapolations
        # Would be better to prevent extrapolation

        if ie <= 1:
            gain = gve[1] - 1.0e-6 * (edrive - 2.0e6) * (gve[0] - gve[1])
            etadrv = eve[1] - 1.0e-6 * (edrive - 2.0e6) * (eve[0] - eve[1])

        elif ie >= 9:
            gain = gve[8] + 1.0e-6 * (edrive - 9.0e6) * (gve[9] - gve[8])
            etadrv = eve[8] + 1.0e-6 * (edrive - 9.0e6) * (eve[9] - eve[8])

        else:
            gain = gve[ie - 1] + de * (gve[ie] - gve[ie - 1])
            etadrv = eve[ie - 1] + de * (eve[ie] - eve[ie - 1])

        # Ensure sensible values

        gain = max(0.01e0, gain)
        etadrv = max(0.01e0, etadrv)

        return gain, etadrv

    def iondrv(self, edrive):
        """Routine to calculate parameters of a heavy ion driver
        suitable for inertial fusion energy

        Parameters
        ----------
        edrive:
            Driver energy (J)

        Returns
        -------
        gain   : output
            Target gain
        etadrv : output
            Driver efficiency

        Notes
        -----
        This routine calculates the parameters of a heavy ion driver
        suitable for inertial fusion energy.

        Heavy-ion Driver Design and Scaling, R. Bieri et al.,
        Fusion Technology, vol.21 (1992) 1583
        Meier and Bieri, Fusion Technology, vol.21 (1992) 1547
        """
        # NOTE: TN and SM agree that the "incomplete" IONDRV model should
        # be removed. It is impossible to use without modifications to the
        # source code indicating is probably should not be used.
        # It can be found in the Git history before commit
        # cc4dd5091b2967673d3b4c7047cd54dd5a016433
        # (30/10/2024)

        # gve(k): target gain at edrive = k MegaJoules
        gve = [25.0, 44.0, 62.0, 76.0, 87.0, 97.0, 107.0, 115.0, 125.0, 132.0]

        # eve(k): driver efficiency at edrive = k MegaJoules
        eve = [0.232, 0.256, 0.269, 0.276, 0.282, 0.286, 0.29, 0.292, 0.294, 0.296]

        e = 1.0e-6 * edrive
        ie = int(e)
        de = e - ie

        # Assume linear interpolations and extrapolations
        # Would be better to prevent extrapolation

        if ie <= 1:
            gain = gve[1] - 1.0e-6 * (edrive - 2.0e6) * (gve[0] - gve[1])
            etadrv = eve[1] - 1.0e-6 * (edrive - 2.0e6) * (eve[0] - eve[1])

        elif ie >= 9:
            gain = gve[8] + 1.0e-6 * (edrive - 9.0e6) * (gve[9] - gve[8])
            etadrv = eve[8] + 1.0e-6 * (edrive - 9.0e6) * (eve[9] - eve[8])

        else:
            gain = gve[ie - 1] + de * (gve[ie] - gve[ie - 1])
            etadrv = eve[ie - 1] + de * (eve[ie] - eve[ie - 1])

        # Ensure sensible values
        gain = max(0.01, gain)
        etadrv = max(0.01, etadrv)

        return gain, etadrv

    def ifestr(self):
        """Routine to calculate the support structural masses for the core of
        an Inertial Fusion Energy power plant

        This routine calculates the support structural masses for the core of
        an Inertial Fusion Energy power plant.

        The outputs are all, trivially, 0 as they are magnetic fusion specific.
        """
        self.data.structure.aintmass = 0.0
        self.data.structure.clgsmass = 0.0
        self.data.structure.coldmass = 0.0
        self.data.structure.fncmass = 0.0
        self.data.structure.gsmass = 0.0

    def ifetgt(self):
        """Routine to calculate the power requirements of the target
        delivery system and the target factory

        This routine calculates the power requirements of the target
        delivery system and the target factory, for an Inertial
        Fusion Energy power plant.
        F/MI/PJK/LOGBOOK12, pp.87-88
        """
        # Target delivery system power (MWe) - effectively negligible

        # Target factory power (MWe)
        # Assumed to scale with repetition rate (not quite linearly)
        self.data.ife.tfacmw = self.data.ife.ptargf * (self.data.ife.reprat / 6.0) ** 0.7

    def ifefbs(self, output: bool = False):
        """Routine to calculate the first wall, blanket and shield volumes,
        masses and other parameters, for an Inertial Fusion Energy device

        This routine calculates the first wall, blanket and shield volumes,
        masses and other parameters, for an Inertial Fusion Energy device.
        F/MI/PJK/LOGBOOK12, p.86
        Moir et al., Fusion Technology, vol.25 (1994) p.5

        Parameters
        ----------
        output: bool
             (Default value = False)
        """
        # Material densities
        # 0 = void
        # 1 = steel
        # 2 = carbon
        # 3 = FLiBe (inferred from Moir et al)
        # 4 = Li2O
        # 5 = concrete
        # 6 = helium (at typical coolant temperatures)
        # 7 = xenon (taken as ten times the normal tabulated value)
        # 8 = lithium (liquid, source Wikipedia)

        matden = [
            0.0,
            self.data.fwbs.den_steel,
            2300.0,
            2020.0,
            2010.0,
            2400.0,
            1.517,
            55.0,
            512.0,
        ]

        # Material masses
        for i in range(MAXMAT + 1):
            den = matden[i]
            self.data.ife.chmatm[i] = self.data.ife.chmatv[i] * den
            for j in range(3):
                self.data.ife.fwmatm[j, i] = self.data.ife.fwmatv[j, i] * den
                self.data.ife.v1matm[j, i] = self.data.ife.v1matv[j, i] * den
                self.data.ife.blmatm[j, i] = self.data.ife.blmatv[j, i] * den
                self.data.ife.v2matm[j, i] = self.data.ife.v2matv[j, i] * den
                self.data.ife.shmatm[j, i] = self.data.ife.shmatv[j, i] * den
                self.data.ife.v3matm[j, i] = self.data.ife.v3matv[j, i] * den

        # Total masses of components (excluding coolant)
        self.data.fwbs.m_fw_total = 0.0
        self.data.fwbs.m_blkt_total = 0.0
        self.data.fwbs.whtshld = 0.0
        for i in range(5):
            for j in range(3):
                self.data.fwbs.m_fw_total += self.data.ife.fwmatm[j, i]
                self.data.fwbs.m_blkt_total += self.data.ife.blmatm[j, i]
                self.data.fwbs.whtshld += self.data.ife.shmatm[j, i]

        # Other masses
        self.data.fwbs.m_blkt_beryllium = 0.0
        self.data.fwbs.m_blkt_vanadium = 0.0
        self.data.fwbs.m_blkt_steel_total = 0.0
        self.data.fwbs.m_blkt_li2o = 0.0
        self.data.fwbs.m_blkt_lithium = 0.0

        for j in range(3):
            self.data.fwbs.m_blkt_steel_total += self.data.ife.blmatm[j, 1]
            self.data.fwbs.m_blkt_li2o += self.data.ife.blmatm[j, 4]
            self.data.fwbs.m_blkt_lithium += self.data.ife.blmatm[j, 8]

        # Total mass of FLiBe
        self.data.ife.mflibe = self.data.ife.chmatm[3]
        for j in range(3):
            self.data.ife.mflibe = (
                self.data.ife.mflibe
                + self.data.ife.fwmatm[j, 3]
                + self.data.ife.v1matm[j, 3]
                + self.data.ife.blmatm[j, 3]
                + self.data.ife.v2matm[j, 3]
                + self.data.ife.shmatm[j, 3]
                + self.data.ife.v3matm[j, 3]
            )

        # A fraction FBREED of the total breeder inventory is outside the
        # core region, i.e. is in the rest of the heat transport system
        if (self.data.ife.fbreed < 0.0) or (self.data.ife.fbreed > 0.999):
            raise ProcessValueError("Illegal fbreed value", fbreed=self.data.ife.fbreed)

        #  Following assumes that use of FLiBe and Li2O are
        # mutually exclusive
        self.data.ife.mflibe /= 1.0 - self.data.ife.fbreed
        self.data.fwbs.m_blkt_li2o /= 1.0 - self.data.ife.fbreed
        self.data.fwbs.m_blkt_lithium /= 1.0 - self.data.ife.fbreed

        # Blanket and first wall lifetimes (HYLIFE-II: = plant life)
        if self.data.ife.ifetyp in {3, 4}:
            life = self.data.costs.life_plant
        else:
            life = min(
                self.data.costs.life_plant,
                self.data.costs.abktflnc
                / (
                    physics_variables.pflux_fw_neutron_mw
                    * self.data.costs.f_t_plant_available
                ),
            )

        self.data.fwbs.life_blkt_fpy = life
        self.data.fwbs.life_fw_fpy = life

        if not output:
            return

        process_output.oheadr(self.outfile, "First Wall, Blanket, Shield")
        process_output.ovarre(
            self.outfile,
            "First wall area (m2)",
            "(a_fw_total)",
            self.data.first_wall.a_fw_total,
        )
        process_output.ovarre(
            self.outfile,
            "First wall mass (kg)",
            "(m_fw_total)",
            self.data.fwbs.m_fw_total,
        )
        process_output.ovarre(
            self.outfile,
            "Blanket mass (kg)",
            "(m_blkt_total)",
            self.data.fwbs.m_blkt_total,
        )
        process_output.ovarre(
            self.outfile,
            "Blanket lithium mass (kg)",
            "(m_blkt_lithium)",
            self.data.fwbs.m_blkt_lithium,
        )
        process_output.ovarre(
            self.outfile,
            "Total mass of FLiBe (kg)",
            "(mflibe)",
            self.data.ife.mflibe,
        )
        process_output.ovarre(
            self.outfile, "Shield mass (kg)", "(whtshld)", self.data.fwbs.whtshld
        )

    def ifepw1(self):
        """Routine to calculate the first part of the heat transport
        and plant power balance constituents, for an IFE power plant
        This routine calculates the first part of the heat transport
        and plant power balance constituents, for an IFE power plant.
        F/MI/PJK/LOGBOOK12, pp.67,89
        Bourque et al., Fusion Technology vol.21 (1992) 1465
        """
        # Driver power reaching target (MW)

        pdrvmw = 1.0e-6 * self.data.ife.pdrive

        # Primary nuclear heating (MW)
        # Total thermal power removed from fusion core

        self.data.heat_transport.priheat = (
            self.data.fwbs.f_p_blkt_multiplication * physics_variables.p_fusion_total_mw
        )

        # Useful (high-grade) thermal power (MW)

        self.data.heat_transport.p_plant_primary_heat_mw = (
            self.data.heat_transport.priheat * (1.0 - self.data.fwbs.fhole)
        )

        # Assume 0.24 of thermal power is intercepted by the first wall
        # (Bourque et al)
        # HYLIFE-II case: Assume FLiBe flows intercept all fusion power
        # and provide the energy multiplication as though it were a
        # conventional blanket

        if self.data.ife.ifetyp not in {3, 4}:
            self.data.heat_transport.p_fw_div_heat_deposited_mw = (
                0.24 * self.data.heat_transport.p_plant_primary_heat_mw
            )
            self.data.fwbs.p_blkt_nuclear_heat_total_mw = (
                self.data.heat_transport.p_plant_primary_heat_mw
                - self.data.heat_transport.p_fw_div_heat_deposited_mw
            )
        else:
            self.data.heat_transport.p_fw_div_heat_deposited_mw = 0.0
            self.data.fwbs.p_blkt_nuclear_heat_total_mw = (
                self.data.heat_transport.p_plant_primary_heat_mw
            )

        self.data.fwbs.p_shld_nuclear_heat_mw = 0.0

        # Lost fusion power (MW)

        self.data.fwbs.pnucloss = (
            self.data.heat_transport.priheat
            - self.data.heat_transport.p_plant_primary_heat_mw
        )  # = priheat*fhole

        # Number of primary heat exchangers

        self.data.heat_transport.n_primary_heat_exchangers = np.ceil(
            self.data.heat_transport.p_plant_primary_heat_mw / 1000.0
        )

        # Secondary heat (some of it... rest calculated in IFEPW2)

        # Wall plug driver power (MW)

        self.data.heat_transport.p_hcd_electric_total_mw = pdrvmw / self.data.ife.etadrv

        # Waste driver power (MW)

        self.data.heat_transport.p_hcd_electric_loss_mw = (
            self.data.heat_transport.p_hcd_electric_total_mw - pdrvmw
        )

        # Cryogenic power (MW)
        # Cryogenic temperature is assumed to be 4.5K

        self.data.heat_transport.p_cryo_plant_electric_mw = self.data.ife.pifecr
        self.data.heat_transport.helpow = (
            1.0e6
            * self.data.heat_transport.p_cryo_plant_electric_mw
            * (0.13 * 4.5)
            / (constants.TEMP_ROOM - 4.5)
        )

    def ifepw2(self, output: bool = False):
        """Routine to calculate the rest of the IFE heat transport
        and plant power balance constituents, not already calculated in
        IFEPW1 or IFEACP

        This routine calculates the rest of the IFE heat transport
        and plant power balance constituents, not already calculated in
        routines <A HREF="ifepw1.html">IFEPW1</A> or
        <A HREF="ifeacp.html">IFEACP</A>.
        F/MI/PJK/LOGBOOK12, p.67

        Parameters
        ----------
        output: bool
             (Default value = False)
        """
        # Facility heat removal (p_plant_electric_base_total_mw calculated in IFEACP)
        self.data.heat_transport.fachtmw = (
            self.data.heat_transport.p_plant_electric_base_total_mw
        )

        # Total secondary heat
        self.data.heat_transport.p_plant_secondary_heat_mw = (
            self.data.heat_transport.p_hcd_electric_loss_mw
            + self.data.fwbs.pnucloss
            + self.data.heat_transport.fachtmw
            + self.data.heat_transport.vachtmw
            + self.data.heat_transport.p_tritium_plant_electric_mw
            + self.data.ife.tdspmw
            + self.data.ife.tfacmw
            + self.data.heat_transport.p_cryo_plant_electric_mw
            + self.data.ife.htpmw_ife
        )

        # Calculate powers relevant to a power-producing plant
        if self.data.costs.ireactor == 1:
            # Gross electric power
            self.data.heat_transport.p_plant_electric_gross_mw = (
                self.data.heat_transport.p_plant_primary_heat_mw
                * self.data.heat_transport.eta_turbine
            )

            # Balance of plant recirculating power fraction
            self.data.heat_transport.fgrosbop = min(
                0.5,
                (
                    self.data.ife.fauxbop
                    / (self.data.heat_transport.p_plant_electric_gross_mw / 1000.0)
                    ** 0.6
                ),
            )

            # Total recirculating power
            self.data.heat_transport.p_plant_electric_recirc_mw = (
                self.data.heat_transport.fgrosbop
                * self.data.heat_transport.p_plant_electric_gross_mw
            ) + self.data.heat_transport.pacpmw

            # Net electric power
            self.data.heat_transport.p_plant_electric_net_mw = (
                self.data.heat_transport.p_plant_electric_gross_mw
                - self.data.heat_transport.p_plant_electric_recirc_mw
            )

            if not output:
                return

            process_output.oheadr(self.outfile, "Power / Heat Transport")
            process_output.ovarre(
                self.outfile,
                "Fusion power escaping via holes (MW)",
                "(pnucloss)",
                self.data.fwbs.pnucloss,
            )
            process_output.ovarre(
                self.outfile,
                "Power multiplication factor",
                "(f_p_blkt_multiplication)",
                self.data.fwbs.f_p_blkt_multiplication,
            )
            if self.data.ife.ifetyp == 4:
                process_output.ovarre(
                    self.outfile, "Tritium Breeding Ratio", "(tbr)", self.data.fwbs.tbr
                )
                process_output.ovarre(
                    self.outfile,
                    "Lithium Fall Time (s)",
                    "(taufall)",
                    self.data.ife.taufall,
                )

            process_output.ovarre(
                self.outfile,
                "Driver wall plug power (MW)",
                "(p_hcd_electric_total_mw)",
                self.data.heat_transport.p_hcd_electric_total_mw,
            )
            process_output.ovarre(
                self.outfile,
                "First wall nuclear heating (MW)",
                "(p_fw_div_heat_deposited_mw)",
                self.data.heat_transport.p_fw_div_heat_deposited_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Blanket nuclear heating (MW)",
                "(p_blkt_nuclear_heat_total_mw)",
                self.data.fwbs.p_blkt_nuclear_heat_total_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Primary heat (MW)",
                "(p_plant_primary_heat_mw)",
                self.data.heat_transport.p_plant_primary_heat_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Secondary heat (MW)",
                "(p_plant_secondary_heat_mw)",
                self.data.heat_transport.p_plant_secondary_heat_mw,
            )
            process_output.oblnkl(self.outfile)
            process_output.ovarre(
                self.outfile,
                "Heat removal from driver power (MW)",
                "(p_hcd_electric_loss_mw)",
                self.data.heat_transport.p_hcd_electric_loss_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from cryogenic plant (MW)",
                "(p_cryo_plant_electric_mw)",
                self.data.heat_transport.p_cryo_plant_electric_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from vacuum pumps (MW)",
                "(vachtmw)",
                self.data.heat_transport.vachtmw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from target factory (MW)",
                "(tfacmw)",
                self.data.ife.tfacmw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from delivery system (MW)",
                "(tdspmw)",
                self.data.ife.tdspmw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from tritium plant (MW)",
                "(p_tritium_plant_electric_mw)",
                self.data.heat_transport.p_tritium_plant_electric_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from facilities (MW)",
                "(fachtmw)",
                self.data.heat_transport.fachtmw,
            )
            process_output.ovarin(
                self.outfile,
                "Number of primary heat exchangers",
                "(n_primary_heat_exchangers)",
                self.data.heat_transport.n_primary_heat_exchangers,
            )

            if self.data.costs.ireactor == 1:
                process_output.osubhd(self.outfile, "Reactor powers :")
                process_output.ovarre(
                    self.outfile,
                    "Gross electric power (MW)",
                    "(p_plant_electric_gross_mw)",
                    self.data.heat_transport.p_plant_electric_gross_mw,
                )
                process_output.ovarre(
                    self.outfile,
                    "Net electric power (MW)",
                    "(p_plant_electric_net_mw)",
                    self.data.heat_transport.p_plant_electric_net_mw,
                )
                process_output.ovarre(
                    self.outfile,
                    "Balance of plant aux. power fraction",
                    "(fgrosbop)",
                    self.data.heat_transport.fgrosbop,
                )

    def ifeacp(self, output: bool = False):
        """Routine to calculate AC power requirements for an IFE power plant

        This routine calculates the AC power requirements for an IFE power plant.
        F/MI/PJK/LOGBOOK12, p.68

        Parameters
        ----------
        output: bool
             (Default value = False)
        """
        # Facility base load, MW (loads not dependent on floor area)

        basemw = self.data.heat_transport.p_plant_electric_base * 1e-6

        # Power needed per floor area, MW/m2

        pmwpm2 = self.data.heat_transport.pflux_plant_floor_electric * 1e-6

        # Total pulsed power system load, MW

        self.data.heat_transport.pacpmw = (
            self.data.heat_transport.p_cryo_plant_electric_mw
            + self.data.heat_transport.vachtmw
            + self.data.ife.tdspmw
            + self.data.ife.tfacmw
            + (self.data.ife.htpmw_ife * self.data.ife.reprat / 6.0)
            + self.data.heat_transport.p_tritium_plant_electric_mw
            + self.data.heat_transport.p_hcd_electric_total_mw
            + basemw
            + (self.data.buildings.a_plant_floor_effective * pmwpm2)
            + self.data.ife.lipmw
        )

        # Total baseline power to facility loads, MW

        self.data.heat_transport.p_plant_electric_base_total_mw = basemw + (
            self.data.buildings.a_plant_floor_effective * pmwpm2
        )

        # Estimate of the total low voltage power, MW

        self.data.heat_transport.tlvpmw = (
            self.data.heat_transport.p_plant_electric_base_total_mw
            + self.data.heat_transport.p_tritium_plant_electric_mw
            + (self.data.ife.htpmw_ife * self.data.ife.reprat / 6.0)
            + self.data.heat_transport.vachtmw
            + 0.5 * self.data.heat_transport.p_cryo_plant_electric_mw
            + self.data.ife.tfacmw
        )

        if not output:
            return

        process_output.oheadr(self.outfile, "AC Power")

        process_output.ovarre(
            self.outfile, "Facility base load (MW)", "(basemw)", basemw
        )
        process_output.ovarre(
            self.outfile,
            "Total floor space (m2)",
            "(a_plant_floor_effective)",
            self.data.buildings.a_plant_floor_effective,
        )
        process_output.ovarre(
            self.outfile, "Power/floor area (MW/m2)", "(pmwpm2)", pmwpm2
        )
        process_output.ovarre(
            self.outfile,
            "Driver power supplies (MW)",
            "(p_hcd_electric_total_mw)",
            self.data.heat_transport.p_hcd_electric_total_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Target delivery system (MW)",
            "(tdspmw)",
            self.data.ife.tdspmw,
        )
        process_output.ovarre(
            self.outfile, "Target factory (MW)", "(tfacmw)", self.data.ife.tfacmw
        )
        process_output.ovarre(
            self.outfile,
            "Tritium processing plant (MW)",
            "(p_tritium_plant_electric_mw)",
            self.data.heat_transport.p_tritium_plant_electric_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Vacuum pump motors (MW)",
            "(vachtmw)",
            self.data.heat_transport.vachtmw,
        )
        process_output.ovarre(
            self.outfile,
            "Cryogenic comp motors (MW)",
            "(p_cryo_plant_electric_mw)",
            self.data.heat_transport.p_cryo_plant_electric_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Heat transport system pump motors (MW)",
            "(htpmw_ife*reprat/6)",
            self.data.ife.htpmw_ife * self.data.ife.reprat / 6.0,
        )
        if self.data.ife.ifetyp == 4:
            process_output.ovarre(
                self.outfile,
                "Lithium Pump Power (MW)",
                "(lipmw)",
                self.data.ife.lipmw,
            )
        process_output.oblnkl(self.outfile)
        process_output.ovarre(
            self.outfile,
            "Total pulsed power (MW)",
            "(pacpmw)",
            self.data.heat_transport.pacpmw,
        )
        process_output.ovarre(
            self.outfile,
            "Total base power reqd at all times (MW)",
            "(p_plant_electric_base_total_mw)",
            self.data.heat_transport.p_plant_electric_base_total_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Total low voltage power (MW)",
            "(tlvpmw)",
            self.data.heat_transport.tlvpmw,
        )

    def ifebdg(self, output: bool = False):
        """Routine to calculate the volumes of the buildings required for
        an Inertial Fusion Energy power plant

        This routine calculates the volumes of the buildings required for
        an Inertial Fusion Energy power plant. The method is based
        closely on that for tokamaks etc. in routine

        Parameters
        ----------
        output: bool
             (Default value = False)
        """
        # Reactor building
        # ================

        # Total internal height
        hrbi = self.data.ife.zl7 + self.data.ife.zu7

        # Distance from centre of device to wall
        self.data.buildings.wrbi = self.data.ife.r7

        # Internal volume (square floor)
        vrci = (2.0 * self.data.buildings.wrbi) ** 2 * hrbi

        # External dimensions
        # RBWT = wall thickness
        # RBRT = roof thickness
        # FNDT = foundation thickness
        rbw = 2.0 * (self.data.ife.r7 + self.data.buildings.rbwt)
        rbl = rbw
        rbh = hrbi + self.data.buildings.rbrt + self.data.buildings.fndt

        # External volume
        rbv = rbw * rbl * rbh

        # Maintenance building
        # ====================

        # The reactor maintenance building includes the hot cells, the
        # decontamination chamber, the transfer corridors, and the waste
        # treatment building.  The dimensions of these areas are scaled
        # from a reference (tokamak) design based on the shield sector size.

        # Shield height

        shh = self.data.ife.zl6 + self.data.ife.zu6

        # Transport corridor size

        tcw = self.data.ife.r6 + 4.0 * self.data.buildings.trcl
        tcl = 5.0 * tcw + 2.0 * self.data.buildings.hcwt

        # Decontamination cell size

        dcw = 2.0 * tcw + 1.0

        # Hot cell size

        hcw = self.data.ife.r6 + 3.0 * self.data.buildings.hccl + 2.0
        hcl = 3.0 * self.data.ife.r6 + 4.0 * self.data.buildings.hccl + tcw

        # Radioactive waste treatment
        # rww = dcw
        # rwl = hcl - dcl - self.data.buildings.hcwt

        # Maintenance building dimensions

        rmbw = hcw + dcw + 3.0 * self.data.buildings.hcwt
        rmbl = hcl + 2.0 * self.data.buildings.hcwt

        # Height

        if self.data.buildings.wgt2 > 1.0:
            wgts = self.data.buildings.wgt2
        else:
            wgts = self.data.fwbs.whtshld

        cran = 9.41e-6 * wgts + 5.1
        rmbh = (
            10.0
            + (self.data.ife.zl6 + self.data.ife.zu6)
            + self.data.buildings.trcl
            + cran
            + 5.1
            + self.data.buildings.stcl
            + self.data.buildings.fndt
        )
        tch = shh + self.data.buildings.stcl + self.data.buildings.fndt

        # Volume

        fac2 = 2.8
        rmbv = fac2 * rmbw * rmbl * rmbh + tcw * tcl * tch

        # Warm shop and hot cell gallery

        wsa = (rmbw + 7.0) * 20.0 + rmbl * 7.0
        fac3 = 1.9
        wsv = fac3 * wsa * rmbh

        # Cryogenic building volume

        cryv = 55.0 * np.sqrt(self.data.heat_transport.helpow)

        # Electrical building volume
        # (set equal to power injection (i.e. driver) building volume)

        elev = self.data.buildings.pibv

        # Calculate effective floor area for ac power module

        self.data.buildings.a_plant_floor_effective = (
            rbv
            + rmbv
            + wsv
            + self.data.buildings.triv
            + elev
            + self.data.buildings.conv
            + cryv
            + self.data.buildings.admv
            + self.data.buildings.shov
        ) / 6.0

        # Convert local into global variables

        self.data.buildings.admvol = self.data.buildings.admv
        self.data.buildings.convol = self.data.buildings.conv
        self.data.buildings.elevol = elev
        self.data.buildings.rbvol = rbv
        self.data.buildings.rmbvol = rmbv
        self.data.buildings.shovol = self.data.buildings.shov
        self.data.buildings.volrci = vrci
        self.data.buildings.wsvol = wsv

        # Total volume of nuclear buildings

        self.data.buildings.volnucb = vrci + rmbv + wsv + self.data.buildings.triv + cryv

        if not output:
            return

        process_output.oheadr(self.outfile, "Plant Buildings System")
        process_output.ovarre(
            self.outfile, "Internal volume of reactor building (m3)", "(vrci)", vrci
        )
        process_output.ovarre(
            self.outfile,
            "Dist from device centre to bldg wall (m)",
            "(wrbi)",
            self.data.buildings.wrbi,
        )
        process_output.ovarre(
            self.outfile,
            "Effective floor area (m2)",
            "(a_plant_floor_effective)",
            self.data.buildings.a_plant_floor_effective,
        )
        process_output.ovarre(self.outfile, "Reactor building volume (m3)", "(rbv)", rbv)
        process_output.ovarre(
            self.outfile, "Reactor maintenance building volume (m3)", "(rmbv)", rmbv
        )
        process_output.ovarre(self.outfile, "Warmshop volume (m3)", "(wsv)", wsv)
        process_output.ovarre(
            self.outfile,
            "Tritium building volume (m3)",
            "(triv)",
            self.data.buildings.triv,
        )
        process_output.ovarre(
            self.outfile, "Electrical building volume (m3)", "(elev)", elev
        )
        process_output.ovarre(
            self.outfile,
            "Control building volume (m3)",
            "(conv)",
            self.data.buildings.conv,
        )
        process_output.ovarre(
            self.outfile, "Cryogenics building volume (m3)", "(cryv)", cryv
        )
        process_output.ovarre(
            self.outfile,
            "Administration building volume (m3)",
            "(admv)",
            self.data.buildings.admv,
        )
        process_output.ovarre(
            self.outfile, "Shops volume (m3)", "(shov)", self.data.buildings.shov
        )
        process_output.ovarre(
            self.outfile,
            "Total volume of nuclear buildings (m3)",
            "(volnucb)",
            self.data.buildings.volnucb,
        )

    def ifevac(self):
        """Routine to calculate parameters of the vacuum system for an
        Inertial Fusion Energy power plant

        This routine calculates the parameters of the vacuum system for an
        Inertial Fusion Energy power plant.
        <P>The calculated values are hard-wired; they are based loosely
        on those for a tokamak of 6m major radius. F/MI/PJK/LOGBOOK12, p.87
        """
        self.data.vacuum.dlscal = 2.0
        self.data.vacuum.n_vv_vacuum_ducts = 16
        self.data.vacuum.m_vv_vacuum_duct_shield = 0.0
        self.data.vacuum.dia_vv_vacuum_ducts = 0.3
        self.data.vacuum.n_vac_pumps_high = 32


def _material_string_generator(chmatv, fwmatv, v1matv, blmatv, v2matv, shmatv, v3matv):
    """Return a function that generates the table rows for the volumes of
    each material in each section of the IFE device.

    `chmatv` provides the volume of each MATERIAL in the chamber.

    All other variables (`fwmatv`, `v1matv`, `blmatv`, `v2matv`, shmatv`, `v3matv`) are
    shaped [3, MATERIAL] and provide the volume of each MATERIAL in the three sections.
    """

    def _material_string(material_index: int):
        fwmatv_sum = fwmatv[:, material_index].sum()
        v1matv_sum = v1matv[:, material_index].sum()
        blmatv_sum = blmatv[:, material_index].sum()
        v2matv_sum = v2matv[:, material_index].sum()
        shmatv_sum = shmatv[:, material_index].sum()
        v3matv_sum = v3matv[:, material_index].sum()

        return (
            MATERIALS[material_index].ljust(11)
            + f"{chmatv[material_index]:.2E}".ljust(10)
            + f"{fwmatv_sum:.2E}".ljust(10)
            + f"{v1matv_sum:.2E}".ljust(10)
            + f"{blmatv_sum:.2E}".ljust(10)
            + f"{v2matv_sum:.2E}".ljust(10)
            + f"{shmatv_sum:.2E}".ljust(10)
            + f"{v3matv_sum:.2E}".ljust(10)
        )

    return _material_string
