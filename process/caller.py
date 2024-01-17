from process import fortran as ft


class Caller:
    """Calls physics and engineering models."""

    def __init__(self, models, x):
        """Initialise all physics and engineering models.

        To ensure that, at the start of a run, all physics/engineering
        variables are fully initialised with consistent values, the models are
        called with the initial optimisation paramters, x.

        :param models: physics and engineering model objects
        :type models: process.main.Models
        :param x: optimisation parameters
        :type x: np.ndarray
        """
        self.models = models
        self.call_models(x)

    def call_models(self, xc):
        """Call the physics and engineering models.

        This method is the principal caller of all the physics and
        engineering models. Some are Fortran subroutines within modules, others
        will be methods on Python model objects.
        :param xc: Array of optimisation parameters
        :type xc: np.array
        """
        # Number of active iteration variables
        nvars = xc.shape[0]

        # Increment the call counter
        ft.numerics.ncalls = ft.numerics.ncalls + 1

        # Convert variables
        ft.define_iteration_variables.convxc(xc, nvars)

        # Perform the various function calls
        # Stellarator caller
        if ft.stellarator_variables.istell != 0:
            self.models.stellarator.run(output=False)
            # TODO Is this return safe?
            return

        # Inertial Fusion Energy calls
        if ft.ife_variables.ife != 0:
            self.models.ife.run(output=False)
            return

        # Tokamak calls
        # Plasma geometry model
        self.models.plasma_geom.geomty()

        # Machine Build Model
        # Radial build
        if ft.build_variables.tf_in_cs == 1:
            self.models.build.tf_in_cs_bore_calc()
        self.models.build.radialb(output=False)

        # Vertical build
        self.models.build.vbuild(output=False)

        self.models.physics.physics()

        # startup model (not used)
        # call startup(ft.constants.nout,0)  !  commented-out for speed reasons

        # Toroidal field coil model

        self.models.tfcoil.run()

        # Toroidal field coil superconductor model
        if ft.tfcoil_variables.i_tf_sup == 1:
            self.models.sctfcoil.run(output=False)

        # Poloidal field and central solenoid model
        self.models.pfcoil.run()

        # Pulsed reactor model
        self.models.pulse.run(output=False)

        # Blanket model
        """Blanket switch values
        No.  |  model
        ---- | ------
        1    |  CCFE HCPB model
        2    |  KIT HCPB model
        3    |  CCFE HCPB model with Tritium Breeding Ratio calculation
        4    |  KIT HCLL model
        5    |  DCLL model
        """
        if ft.fwbs_variables.iblanket == 1:
            # CCFE HCPB model
            self.models.ccfe_hcpb.run(output=False)
        # iblanket = 2, KIT HCPB removed
        elif ft.fwbs_variables.iblanket == 3:
            # CCFE HCPB model with Tritium Breeding Ratio calculation
            self.models.ccfe_hcpb.run(output=False)
            ft.fwbs_variables.tbr = self.models.ccfe_hcpb.tbr_shimwell(
                ft.fwbs_variables.breeder_f,
                ft.fwbs_variables.li6enrich,
                ft.fwbs_variables.iblanket_thickness,
                output=False,
            )
        # iblanket = 4, KIT HCLL removed
        elif ft.fwbs_variables.iblanket == 5:
            # DCLL model
            self.models.dcll.run(output=False)

        self.models.divertor.run(output=False)

        # Structure Model
        self.models.structure.run(output=False)

        # Tight aspect ratio machine model
        if ft.physics_variables.itart == 1 and ft.tfcoil_variables.i_tf_sup != 1:
            self.models.tfcoil.cntrpst()

        # Toroidal field coil power model
        self.models.power.tfpwr(output=False)

        # Poloidal field coil power model
        self.models.power.pfpwr(output=False)

        # Plant heat transport part 1
        self.models.power.power1()

        # Vacuum model
        self.models.vacuum.run(output=False)

        # Buildings model
        self.models.buildings.run(output=False)

        # Plant AC power requirements
        self.models.power.acpow(output=False)

        # Plant heat transport pt 2 & 3
        self.models.power.power2(output=False)

        self.models.power.power3(output=False)

        # Availability model
        self.models.availability.run(output=False)

        # Water usage in secondary cooling system
        self.models.water_use.run(output=False)

        # Costs model
        """Cost switch values
        No.  |  model
        ---- | ------
        0    |  1990 costs model
        1    |  2015 Kovari model
        2    |  2019 STEP model
        """
        if ft.cost_variables.cost_model == 0:
            self.models.costs.run(output=False)
        elif ft.cost_variables.cost_model == 1:
            self.models.costs_2015.run(output=False)
        elif ft.cost_variables.cost_model == 2:
            self.models.costs_step.run()

        # FISPACT and LOCA model (not used)- removed
