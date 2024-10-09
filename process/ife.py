from process.fortran import ife_module as ife, ife_variables, constants, process_output

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
    author: P J Knight, CCFE, Culham Science Centre
    N/A
    This module contains routines for calculating the
    parameters of an Inertial Fusion Energy power plant.
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

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        # write to output file
        if output:
            # Costs
            self.costs.run()
            self.costs.output()

            # Plant availability
            self.availability.avail(output=True)

            # IFE physics
            ife.ifephy(self.outfile, 1)

            # Device build
            self.ifebld(output)

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
        self.ifebld()

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
        self.costs.run()

    def ifebld(self, output: bool = False):
        """Routine to create the build of an inertial fusion energy device
        and to calculate the material volumes for the device core
        author: P J Knight, CCFE, Culham Science Centre

        This routine constructs the build of an inertial fusion energy device
        and calculates the material volumes for the device core.

        :param output: boolean to control writing of output to outfile/mfile
        :type output: bool
        """
        match ife_variables.ifetyp:
            case 1:
                ife.osibld()
            case 2:
                ife.sombld()
            case 3:
                ife.hylbld()
            case 4:
                ife.bld2019()
            case _:
                ife.genbld()

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
        for title, name, thickness, radius in radial_build_data:
            process_output.obuild(self.outfile, title, thickness, radius)

        for title, name, thickness, radius in radial_build_data:
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

            for title, name, thickness, radius in vertical_build_data[:11]:
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
            for title, name, thickness, radius in vertical_build_data:
                process_output.obuild(self.outfile, title, thickness, radius)
        else:
            process_output.oheadr(self.outfile, "Vertical build")
            process_output.write(
                self.outfile, "\t" * 20 + "Thickness (m)" + "\t" * 3 + "Height (m)"
            )
            for title, name, thickness, radius in vertical_build_data:
                process_output.obuild(self.outfile, title, thickness, radius)
            for title, name, thickness, radius in vertical_build_data:
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
