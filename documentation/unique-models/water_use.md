# Power Plant Water Use

The routines documented here can be used to estimate the amount of water used in the secondary cooling system of a  fusion power plant for electricity generation. 

The water *used* by a power plant is defined as that withdrawn from the environment, and includes the water evaporated [^1].

The water evaporated can be estimated from the waste heat remaining after electricity generation: equal to the high-grade heat useful for electricity production (`p_plant_primary_heat_mw`) minus the gross electricity produced (`p_plant_electric_gross_mw`).  

### Cooling Towers

In a cooling tower some heat is absorbed to evaporate the water (latent heat), and some is used to heat the water and the surrounding air (sensible heat).  Ref [^2] gives the evaporation ratio, <em>ER</em>, the ratio of the latent heat absorbed to the total heat discharged through the tower as a function of the ambient air temperature, $T_a$:

$$
ER = 1 - 0.01 (-0.000279 \, {T_a}^3 + 0.00109 \, {T_a}^2 - 0.345 \, {T_a} + 26.7)
$$

For reasons that are not clear the total water usage is given as 1.4 times the evaporated water.  It seems unlikely that leakage could account for this difference.

## Cooling using water bodies

Cooling using lakes or rivers is implemented through either a recirculating system (applicable for cooling ponds), or a once-through system (applicable for rivers) [^3].

The calculation for water evaporated from a cooling water body is based upon heat budget models developed for USGS Report 2013–5188 [^2], in which the evaporation is found from water temperature, wind speed and experimentally-derived wind functions. Refer to that report's spreadsheet model for details: [spreadsheet](https://pubs.usgs.gov/sir/2013/5188/appendix/sir2013-5188_appendix4_fews_version_3.104_edit_20141106.xlsx)

In this estimation, the amount of evaporation from ponds, lakes and rivers is similar. The current routines therefore use an average value across these water bodies to calculate water use.

- For a recirculating water system the total water used is defined as the amount evaporated. 
- For a once-through system the total water used is estimated as 98.3 x the water evaporated.

## Coastal plants

Most nuclear power stations in the UK today are on the coast and use seawater for the condenser.  The water flow rate (usage) will be huge, but very will be evaporated.  The overall water usage is probably negligible.

## Caveats and Limitations

- For simplicity, the properties of water in these calculations are those at 21$^{\circ}$C.  
- Any water used for the rejection of low grade heat is neglected.

(Updated 12/4/23)

[^1]: J Macknick and R Newmark and G Heath and K C Hallett, Operational water consumption and withdrawal factors for electricity generating technologies: a review of existing literature, 2012, Environ. Res. Lett. 7 045802, https://doi.org/10.1088/1748-9326/7/4/045802

[^2]: Diehl, T.H., Harris, M.A., Murphy, J.C., Hutson, S.S., and Ladd, D.E., 2013, Methods for estimating water consumption for thermoelectric power plants in the United States: U.S. Geological Survey Scientific Investigations Report 2013–5188, 78 p., http://dx.doi.org/10.3133/sir20135188

[^3]: Diehl, T.H., and Harris, M.A., 2014, Withdrawal and consumption of water by thermoelectric power plants in the United States, 2010: U.S. Geological Survey Scientific Investigations Report 2014–5184, 28 p., http://dx.doi.org/10.3133/sir20145184
