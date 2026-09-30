"""Runtime-only veqpy equilibrium state (not serialised to IN.DAT or MFILE)."""


class VeqpyRuntime:
    """Holds the latest veqpy equilibrium object and axis iteration results."""

    __slots__ = ("equilibrium", "ne_axis_m3", "te_axis_kev")

    def __init__(self) -> None:
        self.equilibrium = None
        self.ne_axis_m3 = 0.0
        self.te_axis_kev = 0.0

    def clear(self) -> None:
        self.equilibrium = None
        self.ne_axis_m3 = 0.0
        self.te_axis_kev = 0.0
