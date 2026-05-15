"""Water use model for calculating water usage during plant operation."""

import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.model import Model

SECDAY = 86400e0


class WaterUse(Model):
    """Model to calculate water use during plant operation, based on the amount of heat
    that needs to be rejected through the cooling system, and the cooling mechanism
    used.

    The water use is calculated for two different cooling mechanisms: cooling towers,
    and cooling through water bodies (pond, lake, river). The water use for cooling
    through water bodies is calculated for both recirculating and once-through systems.

    The water use for cooling towers is calculated as a single value, as the water use
    for cooling towers is not expected to differ significantly between recirculating
    and once-through systems.
    """

    def __init__(self):
        self.outfile = constants.NOUT

    def output(self):
        """Routine to call the water usage calculation routines and write output to
        file.
        """
        self.run(output=True)

    def run(self, output: bool = False):
        """Routine to call the water usage calculation routines.
        This routine calls the different water usage routines.

        Parameters
        ----------
        output :
            indicate whether output should be written to the output file, or not
        """
        rejected_heat = self.data.heat_transport.p_plant_primary_heat_mw * (
            1 - self.data.heat_transport.eta_turbine
        )

        wastethermeng = rejected_heat * SECDAY

        if output:
            po.oheadr(
                self.outfile, "Water usage during plant operation (secondary cooling)"
            )
            po.ocmmnt(
                self.outfile,
                "Estimated amount of water used through different cooling system "
                "options:",
            )
            po.ocmmnt(self.outfile, "1. Cooling towers")
            po.ocmmnt(
                self.outfile,
                "2. Water bodies (pond, lake, river): recirculating or once-through",
            )

        # call subroutines for cooling mechanisms:

        # cooling towers
        self.cooling_towers(wastethermeng, output=output)

        # water-body cooling
        self.cooling_water_body(wastethermeng, output=output)

    def cooling_towers(self, wastetherm: float, output: bool):
        """Water used in cooling towers

        Parameters
        ----------
        wastetherm:
            thermal energy (MJ) to be cooled by this system
        output:

        """
        self.data.water_use.evapratio = 1.0e0 - (
            (
                -0.000279e0 * self.data.water_use.airtemp**3
                + 0.00109e0 * self.data.water_use.airtemp**2
                - 0.345e0 * self.data.water_use.airtemp
                + 26.7e0
            )
            / 100.0e0
        )
        # Diehl et al. USGS Report 2013-5188, http://dx.doi.org/10.3133/sir20135188

        self.data.water_use.volheat = (
            self.data.water_use.waterdens * self.data.water_use.latentheat
        )

        self.data.water_use.energypervol = (
            self.data.water_use.volheat / self.data.water_use.evapratio
        )

        self.data.water_use.volperenergy = (
            1.0e0 / self.data.water_use.energypervol * 1000000.0e0
        )

        self.data.water_use.evapvol = wastetherm * self.data.water_use.volperenergy

        # find water withdrawn from external source
        self.data.water_use.waterusetower = 1.4e0 * self.data.water_use.evapvol
        # Estimated as a ratio to evaporated water (averaged across observed dataset)
        #  as per Diehl et al. USGS Report 2014-5184, http://dx.doi.org/10.3133/sir20145184

        # end break

        #  Output section
        if output:
            po.ovarre(
                self.outfile,
                "Volume used in cooling tower (m3/day)",
                "(waterusetower)",
                self.data.water_use.waterusetower,
                "OP ",
            )

    def cooling_water_body(self, wastetherm: float, output: bool):
        """Water evaporated in cooling through water bodies
        Based on spreadsheet from Diehl et al. USGS Report 2013-5188, which includes
        cooling coefficients found through fits across a dataset containing a wide range
        of temperatures, windspeeds, and heat loading:
        http://pubs.usgs.gov/sir/2013/5188/appendix/sir2013-5188_appendix4_fews_version_3.104.xlsx


        Parameters
        ----------
        wastetherm:
            thermal energy (MJ) to be cooled by this system
        output:

        """
        evapsum = 0.0e0

        for icool in range(1, 4):
            if icool == 1:
                # small pond as a cooling body
                # heat loading, MW/acre, based on estimations from US power plants
                heatload = 0.35e0
                # coefficients as per Brady et al. 1969:
                # wind function coefficients
                a = 2.47e0
                b = 0e0
                c = 0.12e0
                # fitted coefficients of heat loading
                d = 3061.331e0
                e = -48.810e0
                f = -78.559e0
                g = -291.820e0
                h = 0.267e0
                i = -0.610e0
                j = 33.497e0

            elif icool == 2:
                # large lake or reservoir as a cooling body
                # heat loading, MW/acre, based on estimations from US power plants
                heatload = 0.10e0
                # coefficients as per Webster et al. 1995:
                # wind function coefficients
                a = 1.04e0
                b = 1.05e0
                c = 0.0e0
                # fitted coefficients of heat loading
                d = 3876.843e0
                e = -49.071e0
                f = -295.246e0
                g = -327.935e0
                h = 0.260e0
                i = 10.528e0
                j = 40.188e0

            elif icool == 3:
                # stream or river as a cooling body
                # heat loading, MW/acre, based on estimations from US power plants
                heatload = 0.20e0
                # coefficients as per Gulliver et al. 1986:
                # wind function coefficients
                a = 2.96e0
                b = 0.64e0
                c = 0.0e0
                # fitted coefficients of heat loading
                d = 2565.009e0
                e = -43.636e0
                f = -93.834e0
                g = -203.767e0
                h = 0.257e0
                i = 2.408e0
                j = 20.596e0

            # Unfortunately, the source spreadsheet was from the US, so the fits for
            # water body heating due to heat loading and the cooling wind functions
            # are in non-metric units, hence the conversions required here.
            # Limitations: maximum wind speed of ~5 m/s; initial
            # self.data.water_use.watertemp < 25 degC

            # convert self.data.water_use.windspeed to mph
            self.data.water_use.windspeedmph = self.data.water_use.windspeed * 2.237e0

            # convert heat loading into cal/(cm2.sec)
            heatloadimp = heatload * 1000000.0e0 * 0.239e0 / 40469000.0e0

            # estimate how heat loading will raise temperature, for this water body
            heatratio = (
                d
                + (e * self.data.water_use.watertemp)
                + (f * self.data.water_use.windspeedmph)
                + (g * heatload)
                + (h * self.data.water_use.watertemp**2)
                + (i * self.data.water_use.windspeedmph**2)
                + (j * heatload**2)
            )

            # estimate resultant heated water temperature
            self.data.water_use.watertempheated = self.data.water_use.watertemp + (
                heatloadimp * heatratio
            )

            # find wind function, m/(day.kPa), applicable to this water body:
            windfunction = (
                a
                + (b * self.data.water_use.windspeed)
                + (c * self.data.water_use.windspeed**2)
            ) / 1000.0e0

            # difference in saturation vapour pressure (Clausius-Clapeyron approximation)
            satvapdelta = (
                0.611e0
                * np.exp(
                    (17.27e0 * self.data.water_use.watertempheated)
                    / (237.3e0 + self.data.water_use.watertempheated)
                )
            ) - (
                0.611e0
                * np.exp(
                    (17.27e0 * self.data.water_use.watertemp)
                    / (237.3e0 + self.data.water_use.watertemp)
                )
            )

            # find 'forced evaporation' driven by heat inserted into system
            deltae = (
                self.data.water_use.waterdens
                * self.data.water_use.latentheat
                * windfunction
                * satvapdelta
            )

            # convert heat loading to J/(m2.day)
            heatloadmet = heatload * 1000000.0e0 / 4046.85642e0 * SECDAY

            # find evaporation ratio: ratio of the heat used to evaporate water
            #   to the total heat discharged through the tower
            self.data.water_use.evapratio = deltae / heatloadmet
            # Diehl et al. USGS Report 2013-5188, http://dx.doi.org/10.3133/sir20135188

            self.data.water_use.volheat = (
                self.data.water_use.waterdens * self.data.water_use.latentheat
            )

            self.data.water_use.energypervol = (
                self.data.water_use.volheat / self.data.water_use.evapratio
            )

            self.data.water_use.volperenergy = (
                1.0e0 / self.data.water_use.energypervol * 1000000.0e0
            )

            self.data.water_use.evapvol = wastetherm * self.data.water_use.volperenergy

            # using this method the estimates for pond, lake and river evaporation
            # produce similar results, the average will be taken and used in the next
            # stage of calculation
            evapsum += self.data.water_use.evapvol

        evapsum /= icool

        # water volume withdrawn from external source depends on recirculation or
        # 'once-through' system choice. Estimated as a ratio to evaporated water
        # (averaged across observed dataset) as per Diehl et al. USGS Report 2014-5184,
        # http://dx.doi.org/10.3133/sir20145184

        # recirculating water system:
        self.data.water_use.wateruserecirc = 1.0e0 * evapsum

        # once-through water system:
        self.data.water_use.wateruseonethru = 98.0e0 * evapsum

        # end break

        #  Output section
        if output:
            po.ovarre(
                self.outfile,
                "Volume used in recirculating water system (m3/day)",
                "(wateruserecirc)",
                self.data.water_use.wateruserecirc,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Volume used in once-through water system (m3/day)",
                "(wateruseonethru)",
                self.data.water_use.wateruseonethru,
                "OP ",
            )
