import logging

from process.core import constants
from process.core.model import Model

logger = logging.getLogger(__name__)


class PlasmaFuelling(Model):
    """Class to hold plasma fuelling calculations for plasma processing."""

    def __init__(self):

        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def output(self):
        """This model doesn't output to the output file, but it does generate contour
        plots of plasma fuel flow rates vs recycling and fuelling efficiency.
        """

    @staticmethod
    def calculate_fuel_burnup_fraction(
        fusrat_total: float, molflow_plasma_fuelling_vv_injected: float
    ) -> float:
        """Calculate the fuel burnup fraction

        Parameters
        ----------
        fusrat_total : float
            Total fusion rate (particles/s).
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate into vacuum vessel (particles/s).

        Returns
        -------
        Fuel burnup fraction (dimensionless).

        Notes
        -----
        The fusion rate is multiplied by two to convert from nucleus pairs to particles,
        as the fuelling rate is in particles/s.

        """
        return 2 * fusrat_total / molflow_plasma_fuelling_vv_injected

    @staticmethod
    def calculate_tritium_burnup_fraction(
        fusrat_dt_total: float,
        molflow_plasma_fuelling_vv_injected: float,
        f_molflow_plasma_fuelling_tritium: float,
    ) -> float:
        """Calculate the tritium burnup fraction

        Parameters
        ----------
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate into vacuum vessel (particles/s).
        f_molflow_plasma_fuelling_tritium : float
            Fraction of tritium in the plasma fuelling.

        Returns
        -------
        Tritium burnup fraction (dimensionless).

        Notes
        -----
        The fusion rate is multiplied by two to convert from nucleus pairs to particles,
        as the fuelling rate is in particles/s.

        """
        return fusrat_dt_total / (
            molflow_plasma_fuelling_vv_injected * f_molflow_plasma_fuelling_tritium
        )

    @staticmethod
    def calculate_deuterium_burnup_fraction(
        fusrat_dt_total: float,
        molflow_plasma_fuelling_vv_injected: float,
        f_molflow_plasma_fuelling_deuterium: float,
        fusrat_plasma_dd_total: float,
        fusrat_plasma_dhe3: float,
    ) -> float:
        """Calculate the deuterium burnup fraction

        Parameters
        ----------
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate into vacuum vessel (particles/s).
        f_molflow_plasma_fuelling_deuterium : float
            Fraction of deuterium in the plasma fuelling.
        fusrat_plasma_dd_total : float
            Total deuterium consumption rate from DD fusion (particles/s).
        fusrat_plasma_dhe3 : float
            Deuterium consumption rate from D-He3 fusion (particles/s).

        Returns
        -------
        Deuterium burnup fraction (dimensionless).

        Notes
        -----
        The fusion rate is multiplied by two to convert from nucleus pairs to particles,
        as the fuelling rate is in particles/s.

        """
        return (fusrat_dt_total + 2 * fusrat_plasma_dd_total + fusrat_plasma_dhe3) / (
            molflow_plasma_fuelling_vv_injected * f_molflow_plasma_fuelling_deuterium
        )

    @staticmethod
    def calculate_plasma_tritium_flow_rate(
        f_molflow_plasma_fuelling_tritium: float,
        eta_plasma_fuelling: float,
        molflow_plasma_fuelling_vv_injected: float,
        fusrat_dt_total: float,
        fusrat_plasma_dd_triton: float,
        t_energy_confinement: float,
        f_plasma_particles_lcfs_recycled: float,
        nd_plasma_fuel_ions_vol_avg: float,
        vol_plasma: float,
        f_plasma_fuel_tritium: float,
    ) -> float:
        """Calculate the tritium flow rate in the plasma exhaust.

        Parameters
        ----------
        f_molflow_plasma_fuelling_tritium : float
            Fraction of tritium in the plasma fuelling.
        eta_plasma_fuelling : float
            Fuelling rate efficiency.
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate (particles/s).
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        fusrat_plasma_dd_triton : float
            Tritium production rate from DD fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_plasma_particles_lcfs_recycled : float
            Fraction of plasma particles recycled at the LCFS.
        nd_plasma_fuel_ions_vol_avg : float
            Volume-averaged density of fuel ions in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).
        f_plasma_fuel_tritium : float
            Fraction of tritium in the plasma fuel.

        Returns
        -------
        float
            Tritium flow rate in the plasma exhaust (particles/s).

        """
        return (
            (
                f_molflow_plasma_fuelling_tritium
                * eta_plasma_fuelling
                * molflow_plasma_fuelling_vv_injected
            )
            - fusrat_dt_total
            + fusrat_plasma_dd_triton
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_tritium)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_deuterium_flow_rate(
        f_molflow_plasma_fuelling_deuterium: float,
        eta_plasma_fuelling: float,
        molflow_plasma_fuelling_vv_injected: float,
        fusrat_dt_total: float,
        fusrat_plasma_dhe3: float,
        fusrat_plasma_dd_total: float,
        t_energy_confinement: float,
        f_plasma_particles_lcfs_recycled: float,
        nd_plasma_fuel_ions_vol_avg: float,
        vol_plasma: float,
        f_plasma_fuel_deuterium: float,
    ) -> float:
        """Calculate the deuterium flow rate in the plasma exhaust.

        Parameters
        ----------
        f_molflow_plasma_fuelling_deuterium : float
            Fraction of deuterium in the plasma fuelling.
        eta_plasma_fuelling : float
            Fuelling rate efficiency.
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate (particles/s).
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        fusrat_plasma_dhe3 : float
            Deuterium consumption rate from D-He3 fusion (particles/s).
        fusrat_plasma_dd_total : float
            Total deuterium consumption rate from DD fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_plasma_particles_lcfs_recycled : float
            Fraction of plasma particles recycled at the LCFS.
        nd_plasma_fuel_ions_vol_avg : float
            Volume-averaged density of fuel ions in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).
        f_plasma_fuel_deuterium : float
            Fraction of deuterium in the plasma fuel.

        Returns
        -------
        float
            Deuterium flow rate in the plasma exhaust (particles/s).


        """
        return (
            (
                f_molflow_plasma_fuelling_deuterium
                * eta_plasma_fuelling
                * molflow_plasma_fuelling_vv_injected
            )
            - fusrat_dt_total
            - 2 * fusrat_plasma_dd_total
            - fusrat_plasma_dhe3
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_deuterium)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_helium3_flow_rate(
        f_molflow_plasma_fuelling_helium3: float,
        eta_plasma_fuelling: float,
        molflow_plasma_fuelling_vv_injected: float,
        fusrat_plasma_dhe3: float,
        t_energy_confinement: float,
        f_plasma_particles_lcfs_recycled: float,
        nd_plasma_fuel_ions_vol_avg: float,
        vol_plasma: float,
        f_plasma_fuel_helium3: float,
    ) -> float:
        """Calculate the helium-3 flow rate in the plasma exhaust.

        Parameters
        ----------
        f_molflow_plasma_fuelling_helium3 : float
            Fraction of helium-3 in the plasma fuelling.
        eta_plasma_fuelling : float
            Fuelling rate efficiency.
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate (particles/s).
        fusrat_plasma_dhe3 : float
            Deuterium consumption rate from D-He3 fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_plasma_particles_lcfs_recycled : float
            Fraction of plasma particles recycled at the LCFS.
        nd_plasma_fuel_ions_vol_avg : float
            Volume-averaged density of fuel ions in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).
        f_plasma_fuel_helium3 : float
            Fraction of helium-3 in the plasma fuel.

        Returns
        -------
        float
            Helium-3 flow rate in the plasma exhaust (particles/s).

        """
        return (
            (
                f_molflow_plasma_fuelling_helium3
                * eta_plasma_fuelling
                * molflow_plasma_fuelling_vv_injected
            )
            + fusrat_plasma_dhe3
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_helium3)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_alphas_flow_rate(
        fusrat_dt_total: float,
        fusrat_plasma_dhe3: float,
        t_energy_confinement: float,
        f_t_alpha_energy_confinement: float,
        nd_plasma_alphas_vol_avg: float,
        vol_plasma: float,
    ) -> float:
        """Calculate the alpha particle flow rate in the plasma exhaust.

        Parameters
        ----------
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        fusrat_plasma_dhe3 : float
            Deuterium consumption rate from D-He3 fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_t_alpha_energy_confinement : float
            Ratio of alpha particle confinement time to energy confinement time (dimensionless).
        nd_plasma_alphas_vol_avg : float
            Volume-averaged density of alpha particles in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).

        Returns
        -------
        float
            Alpha particle flow rate in the plasma exhaust (particles/s).

        """
        # Alpha particle balance

        return (
            fusrat_dt_total
            + fusrat_plasma_dhe3
            - (nd_plasma_alphas_vol_avg * vol_plasma)
            / (t_energy_confinement * f_t_alpha_energy_confinement)
        )
