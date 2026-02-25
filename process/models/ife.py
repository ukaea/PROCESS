"""Module containing Inertial Fusion Energy device routines


This module contains routines for calculating the
parameters of an Inertial Fusion Energy power plant.
"""

import numpy as np

from process import process_output
from process.core import constants
from process.core.exceptions import ProcessValueError
from process.data_structure import (
    buildings_variables,
    cost_variables,
    first_wall_variables,
    fwbs_variables,
    heat_transport_variables,
    ife_variables,
    physics_variables,
    structure_variables,
    vacuum_variables,
)

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


class IFE:
    """Module containing Inertial Fusion Energy device routines

    N/A
    This module contains routines for calculating the
    parameters of an Inertial Fusion Energy power plant.
    """

    def __init__(self, availability, costs):
        """Initialises the IFE module's variables

        :param availability: a pointer to the availability model, allowing use of availability's variables/methods
        :type availability: process.availability.Availability
        :param costs: a pointer to the costs model, allowing the use of costs' variables/methods
        :type costs: process.costs.Costs
        """

        self.outfile: int = constants.NOUT
        self.availability = availability
        self.costs = costs

    def run(self, output: bool):
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
        match ife_variables.ifetyp:
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
            ["Chamber", "chrad", ife_variables.chrad, ife_variables.r1],
            ["First Wall", "fwdr", ife_variables.fwdr, ife_variables.r2],
            ["Void 1", "v1dr", ife_variables.v1dr, ife_variables.r3],
            ["Blanket", "bldr", ife_variables.bldr, ife_variables.r4],
            ["Void 2", "v2dr", ife_variables.v2dr, ife_variables.r5],
            ["Shield", "shdr", ife_variables.shdr, ife_variables.r6],
            ["Void 3", "v3dr", ife_variables.v3dr, ife_variables.r7],
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
            ["Base of device", None, 0.0, -ife_variables.zl7],
            ["Void 3 lower", "v3dzl", ife_variables.v3dzl, -ife_variables.zl6],
            ["Shield lower", "shdzl", ife_variables.shdzl, -ife_variables.zl5],
            ["Void 2 lower", "v2dzl", ife_variables.v2dzl, -ife_variables.zl4],
            ["Blanket lower", "bldzl", ife_variables.bldzl, -ife_variables.zl3],
            ["Void 1 lower", "v1dzl", ife_variables.v1dzl, -ife_variables.zl2],
            ["First wall lower", "fwdzl", ife_variables.fwdzl, -ife_variables.zl1],
            ["Chamber lower", "chdzl", ife_variables.chdzl, 0.0],
            ["Chamber upper", "chdzu", ife_variables.chdzu, ife_variables.zu1],
            ["First wall upper", "fwdzu", ife_variables.fwdzu, ife_variables.zu2],
            ["Void 1 upper", "v1dzu", ife_variables.v1dzu, ife_variables.zu3],
            ["Blanket upper", "bldzu", ife_variables.bldzu, ife_variables.zu4],
            ["Void 2 upper", "v2dzu", ife_variables.v2dzu, ife_variables.zu5],
            ["Shield upper", "shdzu", ife_variables.shdzu, ife_variables.zu6],
            ["Void 3 upper", "v3dzu", ife_variables.v3dzu, ife_variables.zu7],
        ]

        if ife_variables.ifetyp == 4:
            process_output.oheadr(self.outfile, "Vertical build - Midplane")
            process_output.write(
                self.outfile, "\t" * 20 + "Thickness (m)" + "\t" * 3 + "Radius (m)"
            )

            for title, _name, thickness, radius in vertical_build_data[:11]:
                process_output.obuild(self.outfile, title, thickness, radius)

            process_output.obuild(
                self.outfile,
                "Blanket upper",
                ife_variables.bldzu - ife_variables.bldzu,
                ife_variables.zu4 - ife_variables.bldzu,
            )
            process_output.obuild(
                self.outfile,
                "Void 2 upper",
                ife_variables.v2dzu,
                ife_variables.zu5 - ife_variables.bldzu,
            )
            process_output.obuild(
                self.outfile,
                "Shield upper",
                ife_variables.shdzu,
                ife_variables.zu6 - ife_variables.bldzu,
            )
            process_output.obuild(
                self.outfile,
                "Void 3 upper",
                ife_variables.v3dzu + ife_variables.bldzu,
                ife_variables.zu7,
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
            ife_variables.chmatv,
            ife_variables.fwmatv,
            ife_variables.v1matv,
            ife_variables.blmatv,
            ife_variables.v2matv,
            ife_variables.shmatv,
            ife_variables.v3matv,
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
        first_wall_variables.a_fw_total = (
            2.0 * np.pi * ife_variables.r1 * (ife_variables.zu1 + ife_variables.zl1)
            + np.pi * ife_variables.r1 * ife_variables.r1
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
        ife_variables.r1 = ife_variables.chrad
        ife_variables.r2 = ife_variables.r1 + ife_variables.fwdr
        ife_variables.r3 = ife_variables.r2 + ife_variables.v1dr
        ife_variables.r4 = ife_variables.r3 + ife_variables.bldr
        ife_variables.r5 = ife_variables.r4 + ife_variables.v2dr
        ife_variables.r6 = ife_variables.r5 + ife_variables.shdr
        ife_variables.r7 = ife_variables.r6 + ife_variables.v3dr

        # Vertical build (below midplane)
        ife_variables.zl1 = ife_variables.chdzl
        ife_variables.zl2 = ife_variables.zl1 + ife_variables.fwdzl
        ife_variables.zl3 = ife_variables.zl2 + ife_variables.v1dzl
        ife_variables.zl4 = ife_variables.zl3 + ife_variables.bldzl
        ife_variables.zl5 = ife_variables.zl4 + ife_variables.v2dzl
        ife_variables.zl6 = ife_variables.zl5 + ife_variables.shdzl
        ife_variables.zl7 = ife_variables.zl6 + ife_variables.v3dzl

        # Vertical build (above midplane)
        ife_variables.zu1 = ife_variables.chdzu
        ife_variables.zu2 = ife_variables.zu1 + ife_variables.fwdzu
        ife_variables.zu3 = ife_variables.zu2 + ife_variables.v1dzu
        ife_variables.zu4 = ife_variables.zu3 + ife_variables.bldzu
        ife_variables.zu5 = ife_variables.zu4 + ife_variables.v2dzu
        ife_variables.zu6 = ife_variables.zu5 + ife_variables.shdzu
        ife_variables.zu7 = ife_variables.zu6 + ife_variables.v3dzu

        # The SOMBRERO chamber is made up of a cylindrical first wall/
        # blanket, with conical regions above and below. Outside this is
        # a cylindrical shield.

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part

        # Chamber : CHCYLH is the height of the cylindrical part
        chcylh = ife_variables.chdzu + ife_variables.chdzl - 2.0 * ife_variables.chrad

        ife_variables.chvol = (
            np.pi
            * ife_variables.r1
            * ife_variables.r1
            * (chcylh + (2.0 / 3.0) * ife_variables.chrad)
        )

        # First wall
        ife_variables.fwvol[0] = (
            np.pi
            * (ife_variables.r2 * ife_variables.r2 - ife_variables.r1 * ife_variables.r1)
            * chcylh
        )
        ife_variables.fwvol[1] = (
            (1.0 / 3.0)
            * np.pi
            * (
                ife_variables.r2
                * ife_variables.r2
                * (ife_variables.chrad + ife_variables.fwdzu)
                - ife_variables.r1 * ife_variables.r1 * ife_variables.chrad
            )
        )
        ife_variables.fwvol[2] = (
            (1.0 / 3.0)
            * np.pi
            * (
                ife_variables.r2
                * ife_variables.r2
                * (ife_variables.chrad + ife_variables.fwdzl)
                - ife_variables.r1 * ife_variables.r1 * ife_variables.chrad
            )
        )

        # First void
        ife_variables.v1vol[0] = (
            np.pi
            * (ife_variables.r3 * ife_variables.r3 - ife_variables.r2 * ife_variables.r2)
            * chcylh
        )
        ife_variables.v1vol[1] = (
            (1.0 / 3.0)
            * np.pi
            * (
                ife_variables.r3
                * ife_variables.r3
                * (ife_variables.chrad + ife_variables.fwdzu + ife_variables.v1dzu)
                - ife_variables.r2
                * ife_variables.r2
                * (ife_variables.chrad + ife_variables.fwdzu)
            )
        )
        ife_variables.v1vol[2] = (
            (1.0 / 3.0)
            * np.pi
            * (
                ife_variables.r3
                * ife_variables.r3
                * (ife_variables.chrad + ife_variables.fwdzl + ife_variables.v1dzl)
                - ife_variables.r2
                * ife_variables.r2
                * (ife_variables.chrad + ife_variables.fwdzl)
            )
        )

        # Blanket:  SOMTDR and SOMBDR are the radii of the cylindrical
        # sections at the top/bottom of the blanket
        # DDZ = Height of top cylindrical section (by similar triangles)
        # DVOL = Volume of top cylindrical section, less the internal cone

        ife_variables.blvol[0] = (
            np.pi
            * (ife_variables.r4 * ife_variables.r4 - ife_variables.r3 * ife_variables.r3)
            * chcylh
        )

        ife_variables.blvol[1] = (
            (1.0 / 3.0)
            * np.pi
            * (
                ife_variables.r4
                * ife_variables.r4
                * (
                    ife_variables.chrad
                    + ife_variables.fwdzu
                    + ife_variables.v1dzu
                    + ife_variables.bldzu
                )
                - ife_variables.r3
                * ife_variables.r3
                * (ife_variables.chrad + ife_variables.fwdzu + ife_variables.v1dzu)
            )
        )
        ddz = (
            (
                ife_variables.chrad
                + ife_variables.fwdzu
                + ife_variables.v1dzu
                + ife_variables.bldzu
            )
            / (
                ife_variables.chrad
                + ife_variables.fwdr
                + ife_variables.v1dr
                + ife_variables.bldr
            )
            * ife_variables.somtdr
        )
        dvol = (
            2.0 * (1.0 / 3.0) * np.pi * ife_variables.somtdr * ife_variables.somtdr * ddz
        )

        ife_variables.blvol[1] = ife_variables.blvol[1] + dvol

        # Ditto for bottom region...

        ife_variables.blvol[2] = (
            (1.0 / 3.0)
            * np.pi
            * (
                ife_variables.r4
                * ife_variables.r4
                * (
                    ife_variables.chrad
                    + ife_variables.fwdzl
                    + ife_variables.v1dzl
                    + ife_variables.bldzl
                )
                - ife_variables.r3
                * ife_variables.r3
                * (ife_variables.chrad + ife_variables.fwdzl + ife_variables.v1dzl)
            )
        )
        ddz = (
            (
                ife_variables.chrad
                + ife_variables.fwdzl
                + ife_variables.v1dzl
                + ife_variables.bldzl
            )
            / (
                ife_variables.chrad
                + ife_variables.fwdr
                + ife_variables.v1dr
                + ife_variables.bldr
            )
            * ife_variables.sombdr
        )
        dvol = (
            2.0 * (1.0 / 3.0) * np.pi * ife_variables.sombdr * ife_variables.sombdr * ddz
        )

        ife_variables.blvol[2] = ife_variables.blvol[2] + dvol

        # Second void
        ife_variables.v2vol[0] = (
            np.pi
            * (ife_variables.r5 * ife_variables.r5 - ife_variables.r4 * ife_variables.r4)
            * chcylh
        )
        ife_variables.v2vol[1] = np.pi * ife_variables.r5 * ife_variables.r5 * (
            ife_variables.zu5 - ife_variables.chdzu + ife_variables.chrad
        ) - (
            ife_variables.fwvol[1]
            + ife_variables.v1vol[1]
            + ife_variables.blvol[1]
            + (
                (1.0 / 3.0)
                * np.pi
                * ife_variables.r1
                * ife_variables.r1
                * ife_variables.chrad
            )
        )
        ife_variables.v2vol[2] = np.pi * ife_variables.r5 * ife_variables.r5 * (
            ife_variables.zl5 - ife_variables.chdzl + ife_variables.chrad
        ) - (
            ife_variables.fwvol[2]
            + ife_variables.v1vol[2]
            + ife_variables.blvol[2]
            + (
                (1.0 / 3.0)
                * np.pi
                * ife_variables.r1
                * ife_variables.r1
                * ife_variables.chrad
            )
        )

        # Shield
        ife_variables.shvol[0] = (
            np.pi
            * (ife_variables.r6 * ife_variables.r6 - ife_variables.r5 * ife_variables.r5)
            * (ife_variables.zu6 + ife_variables.zl6)
        )
        ife_variables.shvol[1] = (
            np.pi
            * ife_variables.r5
            * ife_variables.r5
            * (ife_variables.zu6 - ife_variables.zu5)
        )
        ife_variables.shvol[2] = (
            np.pi
            * ife_variables.r5
            * ife_variables.r5
            * (ife_variables.zl6 - ife_variables.zl5)
        )

        # Third void
        ife_variables.v3vol[0] = (
            np.pi
            * (ife_variables.r7 * ife_variables.r7 - ife_variables.r6 * ife_variables.r6)
            * (ife_variables.zu7 + ife_variables.zl7)
        )
        ife_variables.v3vol[1] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zu7 - ife_variables.zu6)
        )
        ife_variables.v3vol[2] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zl7 - ife_variables.zl6)
        )

        # Material volumes
        for i in range(ife_variables.MAXMAT):
            ife_variables.chmatv[i] = max(
                0.0, ife_variables.chvol * ife_variables.chmatf[i]
            )
            for j in range(3):
                ife_variables.fwmatv[j, i] = max(
                    0.0, ife_variables.fwvol[j] * ife_variables.fwmatf[j, i]
                )
                ife_variables.v1matv[j, i] = max(
                    0.0, ife_variables.v1vol[j] * ife_variables.v1matf[j, i]
                )
                ife_variables.blmatv[j, i] = max(
                    0.0, ife_variables.blvol[j] * ife_variables.blmatf[j, i]
                )
                ife_variables.v2matv[j, i] = max(
                    0.0, ife_variables.v2vol[j] * ife_variables.v2matf[j, i]
                )
                ife_variables.shmatv[j, i] = max(
                    0.0, ife_variables.shvol[j] * ife_variables.shmatf[j, i]
                )
                ife_variables.v3matv[j, i] = max(
                    0.0, ife_variables.v3vol[j] * ife_variables.v3matf[j, i]
                )

        # First wall area
        first_wall_variables.a_fw_total = (
            2.0
            * np.pi
            * ife_variables.r1
            * ((ife_variables.zu1 + ife_variables.zl1) + ife_variables.r1 * np.sqrt(2.0))
        )

    def hylbld(self):
        # Radial build
        ife_variables.r1 = ife_variables.chrad
        ife_variables.r2 = ife_variables.r1 + ife_variables.fwdr
        ife_variables.r3 = ife_variables.r2 + ife_variables.v1dr
        ife_variables.r4 = ife_variables.r3 + ife_variables.bldr
        ife_variables.r5 = ife_variables.r4 + ife_variables.v2dr
        ife_variables.r6 = ife_variables.r5 + ife_variables.shdr
        ife_variables.r7 = ife_variables.r6 + ife_variables.v3dr

        # Vertical build (below midplane)
        ife_variables.zl1 = ife_variables.chdzl
        ife_variables.zl2 = ife_variables.zl1 + ife_variables.fwdzl
        ife_variables.zl3 = ife_variables.zl2 + ife_variables.v1dzl
        ife_variables.zl4 = ife_variables.zl3 + ife_variables.bldzl
        ife_variables.zl5 = ife_variables.zl4 + ife_variables.v2dzl
        ife_variables.zl6 = ife_variables.zl5 + ife_variables.shdzl
        ife_variables.zl7 = ife_variables.zl6 + ife_variables.v3dzl

        # Vertical build (above midplane)
        ife_variables.zu1 = ife_variables.chdzu
        ife_variables.zu2 = ife_variables.zu1 + ife_variables.fwdzu
        ife_variables.zu3 = ife_variables.zu2 + ife_variables.v1dzu
        ife_variables.zu4 = ife_variables.zu3 + ife_variables.bldzu
        ife_variables.zu5 = ife_variables.zu4 + ife_variables.v2dzu
        ife_variables.zu6 = ife_variables.zu5 + ife_variables.shdzu
        ife_variables.zu7 = ife_variables.zu6 + ife_variables.v3dzu

        # The HYLIFE-II chamber is assumed to be mostly cylindrical, but
        # with a conical region below the midplane that causes the Flibe
        # to flow downwards and outwards towards the outlet.

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part
        ife_variables.chvol = (
            np.pi
            * ife_variables.r1
            * ife_variables.r1
            * (
                (ife_variables.zu1 + ife_variables.zl5)
                - (1 / 3) * (ife_variables.zl5 - ife_variables.zl1)
            )
        )

        # First wall
        # FLIRAD is the radius of the Flibe inlet
        ife_variables.fwvol[0] = (
            np.pi
            * (ife_variables.r2 * ife_variables.r2 - ife_variables.r1 * ife_variables.r1)
            * (ife_variables.zu2 + ife_variables.zl5)
        )
        ife_variables.fwvol[1] = (
            np.pi
            * (
                ife_variables.r1 * ife_variables.r1
                - ife_variables.flirad * ife_variables.flirad
            )
            * (ife_variables.zu2 - ife_variables.zu1)
        )
        ife_variables.fwvol[2] = (
            (1 / 3)
            * np.pi
            * (
                ife_variables.r2
                * ife_variables.r2
                * (ife_variables.zl5 - ife_variables.zl1)
                - ife_variables.r1
                * ife_variables.r1
                * (ife_variables.zl5 - ife_variables.zl2)
            )
        )

        # First void
        ife_variables.v1vol[0] = (
            np.pi
            * (ife_variables.r3 * ife_variables.r3 - ife_variables.r2 * ife_variables.r2)
            * (ife_variables.zu2 + ife_variables.zl3)
        )
        ife_variables.v1vol[1] = (
            np.pi
            * (
                ife_variables.r4 * ife_variables.r4
                - ife_variables.flirad * ife_variables.flirad
            )
            * (ife_variables.zu3 - ife_variables.zu2)
        )
        ife_variables.v1vol[2] = (
            (1 / 3)
            * np.pi
            * ife_variables.r1
            * ife_variables.r1
            * (ife_variables.zl3 - ife_variables.zl2)
        )

        # Blanket
        ife_variables.blvol[0] = (
            np.pi
            * (ife_variables.r4 * ife_variables.r4 - ife_variables.r3 * ife_variables.r3)
            * (ife_variables.zu2 + ife_variables.zl3)
        )
        ife_variables.blvol[1] = (
            np.pi
            * (
                ife_variables.r4 * ife_variables.r4
                - ife_variables.flirad * ife_variables.flirad
            )
            * (ife_variables.zu4 - ife_variables.zu3)
        )
        ife_variables.blvol[2] = (
            np.pi
            * ife_variables.r4
            * ife_variables.r4
            * (ife_variables.zl4 - ife_variables.zl3)
        )

        # Second void
        ife_variables.v2vol[0] = (
            np.pi
            * (ife_variables.r5 * ife_variables.r5 - ife_variables.r4 * ife_variables.r4)
            * (ife_variables.zu4 + ife_variables.zl4)
        )
        ife_variables.v2vol[1] = (
            np.pi
            * (
                ife_variables.r5 * ife_variables.r5
                - ife_variables.flirad * ife_variables.flirad
            )
            * (ife_variables.zu5 - ife_variables.zu4)
        )
        ife_variables.v2vol[2] = (
            np.pi
            * ife_variables.r5
            * ife_variables.r5
            * (ife_variables.zl5 - ife_variables.zl4)
        )

        # Shield
        ife_variables.shvol[0] = (
            np.pi
            * (ife_variables.r6 * ife_variables.r6 - ife_variables.r5 * ife_variables.r5)
            * (ife_variables.zu5 + ife_variables.zl5)
        )
        ife_variables.shvol[1] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zu6 - ife_variables.zu5)
        )
        ife_variables.shvol[2] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zl6 - ife_variables.zl5)
        )

        # Third void
        ife_variables.v3vol[0] = (
            np.pi
            * (ife_variables.r7 * ife_variables.r7 - ife_variables.r6 * ife_variables.r6)
            * (ife_variables.zu6 + ife_variables.zl6)
        )
        ife_variables.v3vol[1] = (
            np.pi
            * ife_variables.r7
            * ife_variables.r7
            * (ife_variables.zu7 - ife_variables.zu6)
        )
        ife_variables.v3vol[2] = (
            np.pi
            * ife_variables.r7
            * ife_variables.r7
            * (ife_variables.zl7 - ife_variables.zl6)
        )

        # Material volumes
        for i in range(ife_variables.MAXMAT + 1):
            ife_variables.chmatv[i] = max(
                0.0, ife_variables.chvol * ife_variables.chmatf[i]
            )
            for j in range(3):
                ife_variables.fwmatv[j, i] = max(
                    0.0, ife_variables.fwvol[j] * ife_variables.fwmatf[j, i]
                )
                ife_variables.v1matv[j, i] = max(
                    0.0, ife_variables.v1vol[j] * ife_variables.v1matf[j, i]
                )
                ife_variables.blmatv[j, i] = max(
                    0.0, ife_variables.blvol[j] * ife_variables.blmatf[j, i]
                )
                ife_variables.v2matv[j, i] = max(
                    0.0, ife_variables.v2vol[j] * ife_variables.v2matf[j, i]
                )
                ife_variables.shmatv[j, i] = max(
                    0.0, ife_variables.shvol[j] * ife_variables.shmatf[j, i]
                )
                ife_variables.v3matv[j, i] = max(
                    0.0, ife_variables.v3vol[j] * ife_variables.v3matf[j, i]
                )

        # First wall area
        first_wall_variables.a_fw_total = (
            2.0 * np.pi * ife_variables.r1 * (ife_variables.zu1 + ife_variables.zl5)
        )
        first_wall_variables.a_fw_total = first_wall_variables.a_fw_total + np.pi * (
            ife_variables.r1 * ife_variables.r1
            - ife_variables.flirad * ife_variables.flirad
        )
        first_wall_variables.a_fw_total = (
            first_wall_variables.a_fw_total
            + np.pi
            * ife_variables.r1
            * np.sqrt(
                ife_variables.r1 * ife_variables.r1
                + (ife_variables.zl3 - ife_variables.zl1) ** 2
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
        if ife_variables.fwdr > 0 or ife_variables.v1dr > 0:
            raise ProcessValueError("fwdr and v1dr should be zero for 2019 IFE build")
        if ife_variables.fwdzu > 0 or ife_variables.v1dzu > 0 or ife_variables.v2dzu > 0:
            raise ProcessValueError(
                "fwdzu, v1dzu and v2dzu should be zero for 2019 IFE build"
            )
        if ife_variables.fwdzl > 0 or ife_variables.v1dzl > 0 or ife_variables.v2dzu > 0:
            raise ProcessValueError(
                "fwdzl, v1dzl and v2dzl should be zero for 2019 IFE build"
            )

        # Lithium Pump
        # Velocity
        vel = np.sqrt(
            2.0
            * constants.ACCELERATION_GRAVITY
            * (ife_variables.chdzu + ife_variables.bldzu)
        )

        # Lithium Fraction
        ife_variables.blmatf[0, 8] = 0.91 * np.sqrt(
            ife_variables.bldzu / (ife_variables.chdzu + ife_variables.bldzu)
        )
        ife_variables.blmatf[0, 0] = 1.0 - ife_variables.blmatf[0, 8]

        # Spatial Thickness
        ife_variables.bldr = ife_variables.bldrc / ife_variables.blmatf[0, 8]

        # Area
        acurt = np.pi * (
            (ife_variables.chrad + ife_variables.bldr) ** 2.0 - ife_variables.chrad**2.0
        )

        # Mass Flow
        mdot = 512.0 * vel * ife_variables.blmatf[0, 8] * acurt

        # Pump Power (MW)
        ife_variables.lipmw = (
            1e-6
            * mdot
            * constants.ACCELERATION_GRAVITY
            * (
                ife_variables.chdzl
                + ife_variables.chdzu
                + ife_variables.bldzu
                + ife_variables.bldzl
            )
            / ife_variables.etali
        )

        # Fall Time
        ife_variables.taufall = (
            2.0 * (ife_variables.chdzl + ife_variables.chdzu + ife_variables.bldzu) / vel
        )

        ife_variables.rrmax = 1.0 / ife_variables.taufall

        # TBR and Emult model was for spherical lithium
        # Remove reactor head
        phi = np.arctan(ife_variables.chrad / ife_variables.chdzu)
        sang = 1.0 - np.cos(phi)
        li_frac = 1.0 - 0.5 * sang

        # TBR
        fwbs_variables.tbr = (
            3.7418
            * (1.0 / (1.0 + np.exp(-2.6366 * ife_variables.bldrc)) - 0.5)
            * li_frac
        )

        # Energy Multiplication
        fwbs_variables.f_p_blkt_multiplication = (
            2.2414
            * (1.0 / (1.0 + np.exp(-3.0038 * ife_variables.bldrc)) - 0.5)
            * li_frac
        )

        # Radial build

        ife_variables.r1 = ife_variables.chrad
        ife_variables.r2 = ife_variables.r1 + ife_variables.fwdr
        ife_variables.r3 = ife_variables.r2 + ife_variables.v1dr
        ife_variables.r4 = ife_variables.r3 + ife_variables.bldr
        ife_variables.r5 = ife_variables.r4 + ife_variables.v2dr
        ife_variables.r6 = ife_variables.r5 + ife_variables.shdr
        ife_variables.r7 = ife_variables.r6 + ife_variables.v3dr

        # Vertical build (below midplane)

        ife_variables.zl1 = ife_variables.chdzl
        ife_variables.zl2 = ife_variables.zl1 + ife_variables.fwdzl
        ife_variables.zl3 = ife_variables.zl2 + ife_variables.v1dzl
        ife_variables.zl4 = ife_variables.zl3 + ife_variables.bldzl
        ife_variables.zl5 = ife_variables.zl4 + ife_variables.v2dzl
        ife_variables.zl6 = ife_variables.zl5 + ife_variables.shdzl
        ife_variables.zl7 = ife_variables.zl6 + ife_variables.v3dzl

        # Vertical build (above midplane)

        ife_variables.zu1 = ife_variables.chdzu
        ife_variables.zu2 = ife_variables.zu1 + ife_variables.fwdzu
        ife_variables.zu3 = ife_variables.zu2 + ife_variables.v1dzu
        ife_variables.zu4 = ife_variables.zu3 + ife_variables.bldzu
        ife_variables.zu5 = ife_variables.zu4 + ife_variables.v2dzu
        ife_variables.zu6 = ife_variables.zu5 + ife_variables.shdzu

        ife_variables.v3dzu = (
            (ife_variables.zu6 + ife_variables.zl6)
            + buildings_variables.trcl
            + buildings_variables.stcl
            + 5.1
            + 9.41e-6 * 1.0e5
        )

        ife_variables.zu7 = ife_variables.zu6 + ife_variables.v3dzu

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part

        # Chamber

        chvol = (
            np.pi
            * ife_variables.r1
            * ife_variables.r1
            * (ife_variables.zu1 + ife_variables.zl1)
        )

        # First wall

        ife_variables.fwvol[0] = (
            np.pi
            * (ife_variables.r2 * ife_variables.r2 - ife_variables.r1 * ife_variables.r1)
            * (ife_variables.zu1 + ife_variables.zl1)
        )
        ife_variables.fwvol[1] = (
            np.pi
            * ife_variables.r2
            * ife_variables.r2
            * (ife_variables.zu2 - ife_variables.zu1)
        )
        ife_variables.fwvol[2] = (
            np.pi
            * ife_variables.r2
            * ife_variables.r2
            * (ife_variables.zl2 - ife_variables.zl1)
        )

        # First void

        ife_variables.v1vol[0] = (
            np.pi
            * (ife_variables.r3 * ife_variables.r3 - ife_variables.r2 * ife_variables.r2)
            * (ife_variables.zu2 + ife_variables.zl2)
        )
        ife_variables.v1vol[1] = (
            np.pi
            * ife_variables.r3
            * ife_variables.r3
            * (ife_variables.zu3 - ife_variables.zu2)
        )
        ife_variables.v1vol[2] = (
            np.pi
            * ife_variables.r3
            * ife_variables.r3
            * (ife_variables.zl3 - ife_variables.zl2)
        )

        # Blanket
        # Radial Blanket - between void 2 and chamber
        ife_variables.blvol[0] = (
            np.pi
            * (ife_variables.r4 * ife_variables.r4 - ife_variables.r3 * ife_variables.r3)
            * (ife_variables.zu3 + ife_variables.zl3)
        )
        # Upper Blanket - Pool radially between shield and
        # chamber of input height.
        ife_variables.blvol[1] = (
            np.pi
            * (ife_variables.r5 * ife_variables.r5 - ife_variables.r3 * ife_variables.r3)
            * ife_variables.bldzu
        )
        # Lower Blanket - Pool filling base of device
        ife_variables.blvol[2] = (
            np.pi
            * ife_variables.r5
            * ife_variables.r5
            * (ife_variables.zl4 - ife_variables.zl3)
        )

        # Second void

        ife_variables.v2vol[0] = (
            np.pi
            * (ife_variables.r5 * ife_variables.r5 - ife_variables.r4 * ife_variables.r4)
            * (ife_variables.chdzl + ife_variables.chdzu)
        )
        ife_variables.v2vol[1] = 0.0
        ife_variables.v2vol[2] = 0.0

        # Shield
        ife_variables.shvol[0] = (
            np.pi
            * (ife_variables.r6 * ife_variables.r6 - ife_variables.r5 * ife_variables.r5)
            * (ife_variables.zu5 + ife_variables.zl5)
        )
        # Top Section is in three parts to account for the dip at
        # the centre.  The first is the horizontal top, the second is the
        # horizontal
        ife_variables.shvol[1] = np.pi * (
            (
                (
                    ife_variables.r6 * ife_variables.r6
                    - (ife_variables.chrad - ife_variables.shdr)
                    * (ife_variables.chrad - ife_variables.shdr)
                )
                * ife_variables.shdzu
            )
            + (
                (
                    ife_variables.r1 * ife_variables.r1
                    - ife_variables.flirad * ife_variables.flirad
                )
                * ife_variables.shdzu
            )
            + (
                (
                    ife_variables.r1 * ife_variables.r1
                    - (ife_variables.r1 - ife_variables.shdzu)
                    * (ife_variables.r1 - ife_variables.shdzu)
                )
                * (ife_variables.bldzu - ife_variables.shdzu)
            )
        )
        ife_variables.shvol[2] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zl6 - ife_variables.zl5)
        )

        # Third void

        ife_variables.v3vol[0] = (
            np.pi
            * (ife_variables.r7 * ife_variables.r7 - ife_variables.r6 * ife_variables.r6)
            * (ife_variables.zu6 + ife_variables.zl6)
        )
        ife_variables.v3vol[1] = (
            np.pi
            * ife_variables.r7
            * ife_variables.r7
            * (ife_variables.zu7 - ife_variables.zu6)
            + np.pi
            * (
                (ife_variables.r1 - ife_variables.shdzu)
                * (ife_variables.r1 - ife_variables.shdzu)
                - ife_variables.flirad * ife_variables.flirad
            )
            * ife_variables.bldzu
        )
        ife_variables.v3vol[2] = (
            np.pi
            * ife_variables.r7
            * ife_variables.r7
            * (ife_variables.zl7 - ife_variables.zl6)
        )

        # Material volumes

        for i in range(ife_variables.MAXMAT + 1):
            ife_variables.chmatv[i] = max(0.0, chvol * ife_variables.chmatf[i])
            for j in range(3):
                ife_variables.fwmatv[j, i] = max(
                    0.0, ife_variables.fwvol[j] * ife_variables.fwmatf[j, i]
                )
                ife_variables.v1matv[j, i] = max(
                    0.0, ife_variables.v1vol[j] * ife_variables.v1matf[j, i]
                )
                ife_variables.blmatv[j, i] = max(
                    0.0, ife_variables.blvol[j] * ife_variables.blmatf[j, i]
                )
                ife_variables.v2matv[j, i] = max(
                    0.0, ife_variables.v2vol[j] * ife_variables.v2matf[j, i]
                )
                ife_variables.shmatv[j, i] = max(
                    0.0, ife_variables.shvol[j] * ife_variables.shmatf[j, i]
                )
                ife_variables.v3matv[j, i] = max(
                    0.0, ife_variables.v3vol[j] * ife_variables.v3matf[j, i]
                )

        # First wall area
        # The chamber is surrounded by liquid on three sides
        # with only the top being solid.  This is considered part
        # of the shield. There is a target injector tube at the
        # centre of this area.
        first_wall_variables.a_fw_total = np.pi * (
            ife_variables.r1 * ife_variables.r1
            - ife_variables.flirad * ife_variables.flirad
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

        ife_variables.r1 = ife_variables.chrad
        ife_variables.r2 = ife_variables.r1 + ife_variables.fwdr
        ife_variables.r3 = ife_variables.r2 + ife_variables.v1dr
        ife_variables.r4 = ife_variables.r3 + ife_variables.bldr
        ife_variables.r5 = ife_variables.r4 + ife_variables.v2dr
        ife_variables.r6 = ife_variables.r5 + ife_variables.shdr
        ife_variables.r7 = ife_variables.r6 + ife_variables.v3dr

        # Vertical build (below midplane)

        ife_variables.zl1 = ife_variables.chdzl
        ife_variables.zl2 = ife_variables.zl1 + ife_variables.fwdzl
        ife_variables.zl3 = ife_variables.zl2 + ife_variables.v1dzl
        ife_variables.zl4 = ife_variables.zl3 + ife_variables.bldzl
        ife_variables.zl5 = ife_variables.zl4 + ife_variables.v2dzl
        ife_variables.zl6 = ife_variables.zl5 + ife_variables.shdzl
        ife_variables.zl7 = ife_variables.zl6 + ife_variables.v3dzl

        # Vertical build (above midplane)

        ife_variables.zu1 = ife_variables.chdzu
        ife_variables.zu2 = ife_variables.zu1 + ife_variables.fwdzu
        ife_variables.zu3 = ife_variables.zu2 + ife_variables.v1dzu
        ife_variables.zu4 = ife_variables.zu3 + ife_variables.bldzu
        ife_variables.zu5 = ife_variables.zu4 + ife_variables.v2dzu
        ife_variables.zu6 = ife_variables.zu5 + ife_variables.shdzu
        ife_variables.zu7 = ife_variables.zu6 + ife_variables.v3dzu

        # Component volumes
        # The following notation applies below:
        # J=1 : side part
        # J=2 : top part
        # J=3 : bottom part

        # Chamber

        chvol = (
            np.pi
            * ife_variables.r1
            * ife_variables.r1
            * (ife_variables.zu1 + ife_variables.zl1)
        )

        # First wall

        ife_variables.fwvol[0] = (
            np.pi
            * (ife_variables.r2 * ife_variables.r2 - ife_variables.r1 * ife_variables.r1)
            * (ife_variables.zu1 + ife_variables.zl1)
        )
        ife_variables.fwvol[1] = (
            np.pi
            * ife_variables.r2
            * ife_variables.r2
            * (ife_variables.zu2 - ife_variables.zu1)
        )
        ife_variables.fwvol[2] = (
            np.pi
            * ife_variables.r2
            * ife_variables.r2
            * (ife_variables.zl2 - ife_variables.zl1)
        )

        # First void

        ife_variables.v1vol[0] = (
            np.pi
            * (ife_variables.r3 * ife_variables.r3 - ife_variables.r2 * ife_variables.r2)
            * (ife_variables.zu2 + ife_variables.zl2)
        )
        ife_variables.v1vol[1] = (
            np.pi
            * ife_variables.r3
            * ife_variables.r3
            * (ife_variables.zu3 - ife_variables.zu2)
        )
        ife_variables.v1vol[2] = (
            np.pi
            * ife_variables.r3
            * ife_variables.r3
            * (ife_variables.zl3 - ife_variables.zl2)
        )

        # Blanket

        ife_variables.blvol[0] = (
            np.pi
            * (ife_variables.r4 * ife_variables.r4 - ife_variables.r3 * ife_variables.r3)
            * (ife_variables.zu3 + ife_variables.zl3)
        )
        ife_variables.blvol[1] = (
            np.pi
            * ife_variables.r4
            * ife_variables.r4
            * (ife_variables.zu4 - ife_variables.zu3)
        )
        ife_variables.blvol[2] = (
            np.pi
            * ife_variables.r4
            * ife_variables.r4
            * (ife_variables.zl4 - ife_variables.zl3)
        )

        # Second void

        ife_variables.v2vol[0] = (
            np.pi
            * (ife_variables.r5 * ife_variables.r5 - ife_variables.r4 * ife_variables.r4)
            * (ife_variables.zu4 + ife_variables.zl4)
        )
        ife_variables.v2vol[1] = (
            np.pi
            * ife_variables.r5
            * ife_variables.r5
            * (ife_variables.zu5 - ife_variables.zu4)
        )
        ife_variables.v2vol[2] = (
            np.pi
            * ife_variables.r5
            * ife_variables.r5
            * (ife_variables.zl5 - ife_variables.zl4)
        )

        # Shield

        ife_variables.shvol[0] = (
            np.pi
            * (ife_variables.r6 * ife_variables.r6 - ife_variables.r5 * ife_variables.r5)
            * (ife_variables.zu5 + ife_variables.zl5)
        )
        ife_variables.shvol[1] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zu6 - ife_variables.zu5)
        )
        ife_variables.shvol[2] = (
            np.pi
            * ife_variables.r6
            * ife_variables.r6
            * (ife_variables.zl6 - ife_variables.zl5)
        )

        # Third void

        ife_variables.v3vol[0] = (
            np.pi
            * (ife_variables.r7 * ife_variables.r7 - ife_variables.r6 * ife_variables.r6)
            * (ife_variables.zu6 + ife_variables.zl6)
        )
        ife_variables.v3vol[1] = (
            np.pi
            * ife_variables.r7
            * ife_variables.r7
            * (ife_variables.zu7 - ife_variables.zu6)
        )
        ife_variables.v3vol[2] = (
            np.pi
            * ife_variables.r7
            * ife_variables.r7
            * (ife_variables.zl7 - ife_variables.zl6)
        )

        # Material volumes

        for i in range(ife_variables.MAXMAT + 1):
            ife_variables.chmatv[i] = max(0.0, chvol * ife_variables.chmatf[i])
            for j in range(3):
                ife_variables.fwmatv[j, i] = max(
                    0.0, ife_variables.fwvol[j] * ife_variables.fwmatf[j, i]
                )
                ife_variables.v1matv[j, i] = max(
                    0.0, ife_variables.v1vol[j] * ife_variables.v1matf[j, i]
                )
                ife_variables.blmatv[j, i] = max(
                    0.0, ife_variables.blvol[j] * ife_variables.blmatf[j, i]
                )
                ife_variables.v2matv[j, i] = max(
                    0.0, ife_variables.v2vol[j] * ife_variables.v2matf[j, i]
                )
                ife_variables.shmatv[j, i] = max(
                    0.0, ife_variables.shvol[j] * ife_variables.shmatf[j, i]
                )
                ife_variables.v3matv[j, i] = max(
                    0.0, ife_variables.v3vol[j] * ife_variables.v3matf[j, i]
                )

        # First wall area

        first_wall_variables.a_fw_total = (
            2.0
            * np.pi
            * ife_variables.r1
            * ((ife_variables.zu1 + ife_variables.zl1) + ife_variables.r1)
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
        match ife_variables.ifedrv:
            case -1:
                # Target gain and driver efficiency dependencies on
                # driver energy are input
                ife_variables.gain, ife_variables.etadrv = self.driver(
                    ife_variables.edrive, ife_variables.gainve, ife_variables.etave
                )
            case 0:  # Target gain and driver efficiency are input
                ife_variables.gain = ife_variables.tgain
                ife_variables.etadrv = ife_variables.drveff
            case 1:  # Laser driver based on SOMBRERO design
                ife_variables.gain, ife_variables.etadrv = self.lasdrv(
                    ife_variables.edrive
                )
            case 2:  # Heavy-ion beam driver based on OSIRIS design
                ife_variables.gain, ife_variables.etadrv = self.iondrv(
                    ife_variables.edrive
                )
            case 3:
                ife_variables.etadrv = ife_variables.drveff
            case _:
                raise ProcessValueError(
                    f"ifedrv={ife_variables.ifedrv} is an invalid option"
                )

        if ife_variables.ifedrv != 3:
            # Repetition rate (Hz)
            ife_variables.reprat = ife_variables.pdrive / ife_variables.edrive
            # Fusion power (MW)
            physics_variables.p_fusion_total_mw = (
                1.0e-6 * ife_variables.pdrive * ife_variables.gain
            )
        else:
            # Driver Power
            ife_variables.reprat = ife_variables.rrin
            ife_variables.pdrive = ife_variables.reprat * ife_variables.edrive
            # Gain
            physics_variables.p_fusion_total_mw = ife_variables.pfusife
            ife_variables.gain = physics_variables.p_fusion_total_mw / (
                1.0e-6 * ife_variables.pdrive
            )

        # Wall load (assume total fusion power applies)

        if ife_variables.ifetyp == 1:
            # OSIRIS-type build: First wall subtends a solid angle of 2 pi * SANG

            phi = 0.5 * np.pi + np.arctan(ife_variables.zl1 / ife_variables.r1)
            sang = 1.0 - np.cos(phi)
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_fusion_total_mw
                * 0.5
                * sang
                / first_wall_variables.a_fw_total
            )

        elif ife_variables.ifetyp == 4:
            # 2019 build only has first wall at the top which has a tube at
            # its centre.  This calculates solid angle and removes tube.

            phi = np.arctan(ife_variables.r1 / ife_variables.zu1)
            sang = 1.0 - np.cos(phi)
            phi = np.arctan(ife_variables.flirad / ife_variables.zu1)
            sang = sang - (1.0 - np.cos(phi))
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_fusion_total_mw
                * 0.5
                * sang
                / first_wall_variables.a_fw_total
            )

        else:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_fusion_total_mw / first_wall_variables.a_fw_total
            )

        if not output:
            return

        process_output.oheadr(self.outfile, "Physics / Driver Issues")

        match ife_variables.ifedrv:
            case -1 | 0:
                process_output.ocmmnt(self.outfile, "Driver type : generic")
            case 1:
                process_output.ocmmnt(self.outfile, "Driver type : laser")
            case 2:
                process_output.ocmmnt(self.outfile, "Driver type : heavy ion beam")

        process_output.oblnkl(self.outfile)

        process_output.ovarre(
            self.outfile, "Driver energy (J)", "(edrive)", ife_variables.edrive
        )
        process_output.ovarre(
            self.outfile, "Driver efficiency", "(etadrv)", ife_variables.etadrv
        )
        process_output.ovarre(
            self.outfile,
            "Driver power reaching target (W)",
            "(pdrive)",
            ife_variables.pdrive,
        )
        process_output.ovarre(
            self.outfile,
            "Driver repetition rate (Hz)",
            "(reprat)",
            ife_variables.reprat,
        )
        process_output.ovarre(self.outfile, "Target gain", "(gain)", ife_variables.gain)
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

        structure_variables.aintmass = 0.0
        structure_variables.clgsmass = 0.0
        structure_variables.coldmass = 0.0
        structure_variables.fncmass = 0.0
        structure_variables.gsmass = 0.0

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
        ife_variables.tfacmw = ife_variables.ptargf * (ife_variables.reprat / 6.0) ** 0.7

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
            fwbs_variables.den_steel,
            2300.0,
            2020.0,
            2010.0,
            2400.0,
            1.517,
            55.0,
            512.0,
        ]

        # Material masses
        for i in range(ife_variables.MAXMAT + 1):
            den = matden[i]
            ife_variables.chmatm[i] = ife_variables.chmatv[i] * den
            for j in range(3):
                ife_variables.fwmatm[j, i] = ife_variables.fwmatv[j, i] * den
                ife_variables.v1matm[j, i] = ife_variables.v1matv[j, i] * den
                ife_variables.blmatm[j, i] = ife_variables.blmatv[j, i] * den
                ife_variables.v2matm[j, i] = ife_variables.v2matv[j, i] * den
                ife_variables.shmatm[j, i] = ife_variables.shmatv[j, i] * den
                ife_variables.v3matm[j, i] = ife_variables.v3matv[j, i] * den

        # Total masses of components (excluding coolant)
        fwbs_variables.m_fw_total = 0.0
        fwbs_variables.m_blkt_total = 0.0
        fwbs_variables.whtshld = 0.0
        for i in range(5):
            for j in range(3):
                fwbs_variables.m_fw_total = (
                    fwbs_variables.m_fw_total + ife_variables.fwmatm[j, i]
                )
                fwbs_variables.m_blkt_total = (
                    fwbs_variables.m_blkt_total + ife_variables.blmatm[j, i]
                )
                fwbs_variables.whtshld = (
                    fwbs_variables.whtshld + ife_variables.shmatm[j, i]
                )

        # Other masses
        fwbs_variables.m_blkt_beryllium = 0.0
        fwbs_variables.m_blkt_vanadium = 0.0
        fwbs_variables.m_blkt_steel_total = 0.0
        fwbs_variables.m_blkt_li2o = 0.0
        fwbs_variables.m_blkt_lithium = 0.0

        for j in range(3):
            fwbs_variables.m_blkt_steel_total = (
                fwbs_variables.m_blkt_steel_total + ife_variables.blmatm[j, 1]
            )
            fwbs_variables.m_blkt_li2o = (
                fwbs_variables.m_blkt_li2o + ife_variables.blmatm[j, 4]
            )
            fwbs_variables.m_blkt_lithium = (
                fwbs_variables.m_blkt_lithium + ife_variables.blmatm[j, 8]
            )

        # Total mass of FLiBe
        ife_variables.mflibe = ife_variables.chmatm[3]
        for j in range(3):
            ife_variables.mflibe = (
                ife_variables.mflibe
                + ife_variables.fwmatm[j, 3]
                + ife_variables.v1matm[j, 3]
                + ife_variables.blmatm[j, 3]
                + ife_variables.v2matm[j, 3]
                + ife_variables.shmatm[j, 3]
                + ife_variables.v3matm[j, 3]
            )

        # A fraction FBREED of the total breeder inventory is outside the
        # core region, i.e. is in the rest of the heat transport system
        if (ife_variables.fbreed < 0.0) or (ife_variables.fbreed > 0.999):
            raise ProcessValueError("Illegal fbreed value", fbreed=ife_variables.fbreed)

        #  Following assumes that use of FLiBe and Li2O are
        # mutually exclusive
        ife_variables.mflibe = ife_variables.mflibe / (1.0 - ife_variables.fbreed)
        fwbs_variables.m_blkt_li2o = fwbs_variables.m_blkt_li2o / (
            1.0 - ife_variables.fbreed
        )
        fwbs_variables.m_blkt_lithium = fwbs_variables.m_blkt_lithium / (
            1.0 - ife_variables.fbreed
        )

        # Blanket and first wall lifetimes (HYLIFE-II: = plant life)
        if (ife_variables.ifetyp == 3) or (ife_variables.ifetyp == 4):
            life = cost_variables.life_plant
        else:
            life = min(
                cost_variables.life_plant,
                cost_variables.abktflnc
                / (
                    physics_variables.pflux_fw_neutron_mw
                    * cost_variables.f_t_plant_available
                ),
            )

        fwbs_variables.life_blkt_fpy = life
        fwbs_variables.life_fw_fpy = life

        if not output:
            return

        process_output.oheadr(self.outfile, "First Wall, Blanket, Shield")
        process_output.ovarre(
            self.outfile,
            "First wall area (m2)",
            "(a_fw_total)",
            first_wall_variables.a_fw_total,
        )
        process_output.ovarre(
            self.outfile,
            "First wall mass (kg)",
            "(m_fw_total)",
            fwbs_variables.m_fw_total,
        )
        process_output.ovarre(
            self.outfile,
            "Blanket mass (kg)",
            "(m_blkt_total)",
            fwbs_variables.m_blkt_total,
        )
        process_output.ovarre(
            self.outfile,
            "Blanket lithium mass (kg)",
            "(m_blkt_lithium)",
            fwbs_variables.m_blkt_lithium,
        )
        process_output.ovarre(
            self.outfile,
            "Total mass of FLiBe (kg)",
            "(mflibe)",
            ife_variables.mflibe,
        )
        process_output.ovarre(
            self.outfile, "Shield mass (kg)", "(whtshld)", fwbs_variables.whtshld
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

        pdrvmw = 1.0e-6 * ife_variables.pdrive

        # Primary nuclear heating (MW)
        # Total thermal power removed from fusion core

        heat_transport_variables.priheat = (
            fwbs_variables.f_p_blkt_multiplication * physics_variables.p_fusion_total_mw
        )

        # Useful (high-grade) thermal power (MW)

        heat_transport_variables.p_plant_primary_heat_mw = (
            heat_transport_variables.priheat * (1.0 - fwbs_variables.fhole)
        )

        # Assume 0.24 of thermal power is intercepted by the first wall
        # (Bourque et al)
        # HYLIFE-II case: Assume FLiBe flows intercept all fusion power
        # and provide the energy multiplication as though it were a
        # conventional blanket

        if (ife_variables.ifetyp != 3) and (ife_variables.ifetyp != 4):
            heat_transport_variables.p_fw_div_heat_deposited_mw = (
                0.24 * heat_transport_variables.p_plant_primary_heat_mw
            )
            fwbs_variables.p_blkt_nuclear_heat_total_mw = (
                heat_transport_variables.p_plant_primary_heat_mw
                - heat_transport_variables.p_fw_div_heat_deposited_mw
            )
        else:
            heat_transport_variables.p_fw_div_heat_deposited_mw = 0.0
            fwbs_variables.p_blkt_nuclear_heat_total_mw = (
                heat_transport_variables.p_plant_primary_heat_mw
            )

        fwbs_variables.p_shld_nuclear_heat_mw = 0.0

        # Lost fusion power (MW)

        fwbs_variables.pnucloss = (
            heat_transport_variables.priheat
            - heat_transport_variables.p_plant_primary_heat_mw
        )  # = priheat*fhole

        # Number of primary heat exchangers

        heat_transport_variables.n_primary_heat_exchangers = np.ceil(
            heat_transport_variables.p_plant_primary_heat_mw / 1000.0
        )

        # Secondary heat (some of it... rest calculated in IFEPW2)

        # Wall plug driver power (MW)

        heat_transport_variables.p_hcd_electric_total_mw = pdrvmw / ife_variables.etadrv

        # Waste driver power (MW)

        heat_transport_variables.p_hcd_electric_loss_mw = (
            heat_transport_variables.p_hcd_electric_total_mw - pdrvmw
        )

        # Cryogenic power (MW)
        # Cryogenic temperature is assumed to be 4.5K

        heat_transport_variables.p_cryo_plant_electric_mw = ife_variables.pifecr
        heat_transport_variables.helpow = (
            1.0e6
            * heat_transport_variables.p_cryo_plant_electric_mw
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
        heat_transport_variables.fachtmw = (
            heat_transport_variables.p_plant_electric_base_total_mw
        )

        # Total secondary heat
        heat_transport_variables.p_plant_secondary_heat_mw = (
            heat_transport_variables.p_hcd_electric_loss_mw
            + fwbs_variables.pnucloss
            + heat_transport_variables.fachtmw
            + heat_transport_variables.vachtmw
            + heat_transport_variables.p_tritium_plant_electric_mw
            + ife_variables.tdspmw
            + ife_variables.tfacmw
            + heat_transport_variables.p_cryo_plant_electric_mw
            + ife_variables.htpmw_ife
        )

        # Calculate powers relevant to a power-producing plant
        if cost_variables.ireactor == 1:
            # Gross electric power
            heat_transport_variables.p_plant_electric_gross_mw = (
                heat_transport_variables.p_plant_primary_heat_mw
                * heat_transport_variables.eta_turbine
            )

            # Balance of plant recirculating power fraction
            heat_transport_variables.fgrosbop = min(
                0.5,
                (
                    ife_variables.fauxbop
                    / (heat_transport_variables.p_plant_electric_gross_mw / 1000.0)
                    ** 0.6
                ),
            )

            # Total recirculating power
            heat_transport_variables.p_plant_electric_recirc_mw = (
                heat_transport_variables.fgrosbop
                * heat_transport_variables.p_plant_electric_gross_mw
            ) + heat_transport_variables.pacpmw

            # Net electric power
            heat_transport_variables.p_plant_electric_net_mw = (
                heat_transport_variables.p_plant_electric_gross_mw
                - heat_transport_variables.p_plant_electric_recirc_mw
            )

            if not output:
                return

            process_output.oheadr(self.outfile, "Power / Heat Transport")
            process_output.ovarre(
                self.outfile,
                "Fusion power escaping via holes (MW)",
                "(pnucloss)",
                fwbs_variables.pnucloss,
            )
            process_output.ovarre(
                self.outfile,
                "Power multiplication factor",
                "(f_p_blkt_multiplication)",
                fwbs_variables.f_p_blkt_multiplication,
            )
            if ife_variables.ifetyp == 4:
                process_output.ovarre(
                    self.outfile, "Tritium Breeding Ratio", "(tbr)", fwbs_variables.tbr
                )
                process_output.ovarre(
                    self.outfile,
                    "Lithium Fall Time (s)",
                    "(taufall)",
                    ife_variables.taufall,
                )

            process_output.ovarre(
                self.outfile,
                "Driver wall plug power (MW)",
                "(p_hcd_electric_total_mw)",
                heat_transport_variables.p_hcd_electric_total_mw,
            )
            process_output.ovarre(
                self.outfile,
                "First wall nuclear heating (MW)",
                "(p_fw_div_heat_deposited_mw)",
                heat_transport_variables.p_fw_div_heat_deposited_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Blanket nuclear heating (MW)",
                "(p_blkt_nuclear_heat_total_mw)",
                fwbs_variables.p_blkt_nuclear_heat_total_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Primary heat (MW)",
                "(p_plant_primary_heat_mw)",
                heat_transport_variables.p_plant_primary_heat_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Secondary heat (MW)",
                "(p_plant_secondary_heat_mw)",
                heat_transport_variables.p_plant_secondary_heat_mw,
            )
            process_output.oblnkl(self.outfile)
            process_output.ovarre(
                self.outfile,
                "Heat removal from driver power (MW)",
                "(p_hcd_electric_loss_mw)",
                heat_transport_variables.p_hcd_electric_loss_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from cryogenic plant (MW)",
                "(p_cryo_plant_electric_mw)",
                heat_transport_variables.p_cryo_plant_electric_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from vacuum pumps (MW)",
                "(vachtmw)",
                heat_transport_variables.vachtmw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from target factory (MW)",
                "(tfacmw)",
                ife_variables.tfacmw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from delivery system (MW)",
                "(tdspmw)",
                ife_variables.tdspmw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from tritium plant (MW)",
                "(p_tritium_plant_electric_mw)",
                heat_transport_variables.p_tritium_plant_electric_mw,
            )
            process_output.ovarre(
                self.outfile,
                "Heat removal from facilities (MW)",
                "(fachtmw)",
                heat_transport_variables.fachtmw,
            )
            process_output.ovarin(
                self.outfile,
                "Number of primary heat exchangers",
                "(n_primary_heat_exchangers)",
                heat_transport_variables.n_primary_heat_exchangers,
            )

            if cost_variables.ireactor == 1:
                process_output.osubhd(self.outfile, "Reactor powers :")
                process_output.ovarre(
                    self.outfile,
                    "Gross electric power (MW)",
                    "(p_plant_electric_gross_mw)",
                    heat_transport_variables.p_plant_electric_gross_mw,
                )
                process_output.ovarre(
                    self.outfile,
                    "Net electric power (MW)",
                    "(p_plant_electric_net_mw)",
                    heat_transport_variables.p_plant_electric_net_mw,
                )
                process_output.ovarre(
                    self.outfile,
                    "Balance of plant aux. power fraction",
                    "(fgrosbop)",
                    heat_transport_variables.fgrosbop,
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

        basemw = heat_transport_variables.p_plant_electric_base * 1e-6

        # Power needed per floor area, MW/m2

        pmwpm2 = heat_transport_variables.pflux_plant_floor_electric * 1e-6

        # Total pulsed power system load, MW

        heat_transport_variables.pacpmw = (
            heat_transport_variables.p_cryo_plant_electric_mw
            + heat_transport_variables.vachtmw
            + ife_variables.tdspmw
            + ife_variables.tfacmw
            + (ife_variables.htpmw_ife * ife_variables.reprat / 6.0)
            + heat_transport_variables.p_tritium_plant_electric_mw
            + heat_transport_variables.p_hcd_electric_total_mw
            + basemw
            + (buildings_variables.a_plant_floor_effective * pmwpm2)
            + ife_variables.lipmw
        )

        # Total baseline power to facility loads, MW

        heat_transport_variables.p_plant_electric_base_total_mw = basemw + (
            buildings_variables.a_plant_floor_effective * pmwpm2
        )

        # Estimate of the total low voltage power, MW

        heat_transport_variables.tlvpmw = (
            heat_transport_variables.p_plant_electric_base_total_mw
            + heat_transport_variables.p_tritium_plant_electric_mw
            + (ife_variables.htpmw_ife * ife_variables.reprat / 6.0)
            + heat_transport_variables.vachtmw
            + 0.5 * heat_transport_variables.p_cryo_plant_electric_mw
            + ife_variables.tfacmw
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
            buildings_variables.a_plant_floor_effective,
        )
        process_output.ovarre(
            self.outfile, "Power/floor area (MW/m2)", "(pmwpm2)", pmwpm2
        )
        process_output.ovarre(
            self.outfile,
            "Driver power supplies (MW)",
            "(p_hcd_electric_total_mw)",
            heat_transport_variables.p_hcd_electric_total_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Target delivery system (MW)",
            "(tdspmw)",
            ife_variables.tdspmw,
        )
        process_output.ovarre(
            self.outfile, "Target factory (MW)", "(tfacmw)", ife_variables.tfacmw
        )
        process_output.ovarre(
            self.outfile,
            "Tritium processing plant (MW)",
            "(p_tritium_plant_electric_mw)",
            heat_transport_variables.p_tritium_plant_electric_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Vacuum pump motors (MW)",
            "(vachtmw)",
            heat_transport_variables.vachtmw,
        )
        process_output.ovarre(
            self.outfile,
            "Cryogenic comp motors (MW)",
            "(p_cryo_plant_electric_mw)",
            heat_transport_variables.p_cryo_plant_electric_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Heat transport system pump motors (MW)",
            "(htpmw_ife*reprat/6)",
            ife_variables.htpmw_ife * ife_variables.reprat / 6.0,
        )
        if ife_variables.ifetyp == 4:
            process_output.ovarre(
                self.outfile,
                "Lithium Pump Power (MW)",
                "(lipmw)",
                ife_variables.lipmw,
            )
        process_output.oblnkl(self.outfile)
        process_output.ovarre(
            self.outfile,
            "Total pulsed power (MW)",
            "(pacpmw)",
            heat_transport_variables.pacpmw,
        )
        process_output.ovarre(
            self.outfile,
            "Total base power reqd at all times (MW)",
            "(p_plant_electric_base_total_mw)",
            heat_transport_variables.p_plant_electric_base_total_mw,
        )
        process_output.ovarre(
            self.outfile,
            "Total low voltage power (MW)",
            "(tlvpmw)",
            heat_transport_variables.tlvpmw,
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
        hrbi = ife_variables.zl7 + ife_variables.zu7

        # Distance from centre of device to wall
        buildings_variables.wrbi = ife_variables.r7

        # Internal volume (square floor)
        vrci = (2.0 * buildings_variables.wrbi) ** 2 * hrbi

        # External dimensions
        # RBWT = wall thickness
        # RBRT = roof thickness
        # FNDT = foundation thickness
        rbw = 2.0 * (ife_variables.r7 + buildings_variables.rbwt)
        rbl = rbw
        rbh = hrbi + buildings_variables.rbrt + buildings_variables.fndt

        # External volume
        rbv = rbw * rbl * rbh

        # Maintenance building
        # ====================

        # The reactor maintenance building includes the hot cells, the
        # decontamination chamber, the transfer corridors, and the waste
        # treatment building.  The dimensions of these areas are scaled
        # from a reference (tokamak) design based on the shield sector size.

        # Shield height

        shh = ife_variables.zl6 + ife_variables.zu6

        # Transport corridor size

        tcw = ife_variables.r6 + 4.0 * buildings_variables.trcl
        tcl = 5.0 * tcw + 2.0 * buildings_variables.hcwt

        # Decontamination cell size

        dcw = 2.0 * tcw + 1.0

        # Hot cell size

        hcw = ife_variables.r6 + 3.0 * buildings_variables.hccl + 2.0
        hcl = 3.0 * ife_variables.r6 + 4.0 * buildings_variables.hccl + tcw

        # Radioactive waste treatment
        # rww = dcw
        # rwl = hcl - dcl - buildings_variables.hcwt

        # Maintenance building dimensions

        rmbw = hcw + dcw + 3.0 * buildings_variables.hcwt
        rmbl = hcl + 2.0 * buildings_variables.hcwt

        # Height

        if buildings_variables.wgt2 > 1.0:
            wgts = buildings_variables.wgt2
        else:
            wgts = fwbs_variables.whtshld

        cran = 9.41e-6 * wgts + 5.1
        rmbh = (
            10.0
            + (ife_variables.zl6 + ife_variables.zu6)
            + buildings_variables.trcl
            + cran
            + 5.1
            + buildings_variables.stcl
            + buildings_variables.fndt
        )
        tch = shh + buildings_variables.stcl + buildings_variables.fndt

        # Volume

        fac2 = 2.8
        rmbv = fac2 * rmbw * rmbl * rmbh + tcw * tcl * tch

        # Warm shop and hot cell gallery

        wsa = (rmbw + 7.0) * 20.0 + rmbl * 7.0
        fac3 = 1.9
        wsv = fac3 * wsa * rmbh

        # Cryogenic building volume

        cryv = 55.0 * np.sqrt(heat_transport_variables.helpow)

        # Electrical building volume
        # (set equal to power injection (i.e. driver) building volume)

        elev = buildings_variables.pibv

        # Calculate effective floor area for ac power module

        buildings_variables.a_plant_floor_effective = (
            rbv
            + rmbv
            + wsv
            + buildings_variables.triv
            + elev
            + buildings_variables.conv
            + cryv
            + buildings_variables.admv
            + buildings_variables.shov
        ) / 6.0

        # Convert local into global variables

        buildings_variables.admvol = buildings_variables.admv
        buildings_variables.convol = buildings_variables.conv
        buildings_variables.elevol = elev
        buildings_variables.rbvol = rbv
        buildings_variables.rmbvol = rmbv
        buildings_variables.shovol = buildings_variables.shov
        buildings_variables.volrci = vrci
        buildings_variables.wsvol = wsv

        # Total volume of nuclear buildings

        buildings_variables.volnucb = vrci + rmbv + wsv + buildings_variables.triv + cryv

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
            buildings_variables.wrbi,
        )
        process_output.ovarre(
            self.outfile,
            "Effective floor area (m2)",
            "(a_plant_floor_effective)",
            buildings_variables.a_plant_floor_effective,
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
            buildings_variables.triv,
        )
        process_output.ovarre(
            self.outfile, "Electrical building volume (m3)", "(elev)", elev
        )
        process_output.ovarre(
            self.outfile,
            "Control building volume (m3)",
            "(conv)",
            buildings_variables.conv,
        )
        process_output.ovarre(
            self.outfile, "Cryogenics building volume (m3)", "(cryv)", cryv
        )
        process_output.ovarre(
            self.outfile,
            "Administration building volume (m3)",
            "(admv)",
            buildings_variables.admv,
        )
        process_output.ovarre(
            self.outfile, "Shops volume (m3)", "(shov)", buildings_variables.shov
        )
        process_output.ovarre(
            self.outfile,
            "Total volume of nuclear buildings (m3)",
            "(volnucb)",
            buildings_variables.volnucb,
        )

    def ifevac(self):
        """Routine to calculate parameters of the vacuum system for an
        Inertial Fusion Energy power plant

        This routine calculates the parameters of the vacuum system for an
        Inertial Fusion Energy power plant.
        <P>The calculated values are hard-wired; they are based loosely
        on those for a tokamak of 6m major radius. F/MI/PJK/LOGBOOK12, p.87
        """
        vacuum_variables.dlscal = 2.0
        vacuum_variables.n_vv_vacuum_ducts = 16
        vacuum_variables.m_vv_vacuum_duct_shield = 0.0
        vacuum_variables.dia_vv_vacuum_ducts = 0.3
        vacuum_variables.n_vac_pumps_high = 32


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
