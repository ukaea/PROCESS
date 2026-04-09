
# Armour, First Wall and Breeding Blanket

The surface facing the plasma is a thin layer of a material highly resistant to 
melting and erosion, such as tungsten, referred to "armour". It is cooled by 
conduction to the first wall underneath.

The first wall sits behind the armour, and is dedicated to removing the heat 
landing on the armour. It does not breed tritium. Due to the hostile environment 
the first wall and armour have only a short lifetime and therefore need to be 
replaced regularly. It is cooled either by gaseous helium or by pressurised 
liquid water, depending on the selection of blanket type using the switch 
`blkttype`.

## Wall Load Calculation

Switch `i_pflux_fw_neutron` determines whether the neutron wall load (power per unit area) 
should be calculated using the plasma surface area (`i_pflux_fw_neutron = 1`) or the first 
wall area (`i_pflux_fw_neutron = 2`) as the denominator. In the former case, input 
parameter `ffwal` (default value 0.92) can be used to scale the neutron power 
reaching the first wall.

The breeding blanket performs a number of tasks. An incoming neutron from a
deuterium-tritium (D-T) fusion reaction in the plasma loses energy in the
blanket. This energy is removed by the blanket coolant and used to produce
electricity. The neutron may also react with a lithium nucleus present in the
blanket to produce ("breed") a tritium nucleus which can be re-used as
fuel. The competing requirements of heating and tritium synthesis mean that a
neutron multiplier must be present, to ensure balance between tritium
destruction and creation. The blanket therefore contains beryllium to fulfil
this purpose. As with the first wall, the blanket has a relatively short
lifetime because of the high neutron fluence.

### Blanket Model Options

The models used for the thermoydraulics of the first wall, the profile of 
deposition of the neutron energy, tritium breeding, and conversion of heat to 
electricity have been revised extensively.

`i_blanket_type` -- This switch selects between different types of blanket.

- `== 1` -- CCFE HCPB (helium-cooled pebble bed) model. The energy 
    deposition in the armour and first wall, blanket and shield are calculated 
    using parametric fits to an MCNP neutron and photon transport model of a 
    sector of a tokamak. The blanket contains lithium orthosilicate 
    Li$_4$SiO$_4$, titanium beryllide TiBe$_{12}$, helium and Eurofer steel. 

`i_thermal_electric_conversion` -- This switch controls how the coolant pumping power in the 
first wall and blanket is determined, and also how the calculation of the plant's 
thermal to electric conversion efficiency (the secondary cycle thermal 
efficiency) proceeds.

## Thermo-hydraulic model for first wall and blanket

!!! Note "Note" 
    This is called for i_p_coolant_pumping = 2 and 3

Summary of key variables and switches:  

|                          |       First Wall        | Breeding Blanket Primary | Liquid Breeder/Coolant               |
| :----------------------: | :---------------------: | ------------------------ | ------------------------------------ |
|     Coolant Channels     |      :-----------:      | ------------------------ | --------------------------           |
|        length (m)        |   `len_fw_channel`   | ---                      | ---                                  |
|        width (m)         | `radius_fw_channel` (radius, cicular) | `radius_fw_channel`                    | `a_bz_liq`, `b_bz_liq` (rectangular) |
|    wall thickness (m)    |        `dr_fw_wall`        | dr_fw_wall                  | `th_wall_secondary`                  |
|        dx_fw_module (m)         |         `dx_fw_module`         | ---                      | ---                                  |
|    roughness epsilon     |       `roughness_fw_channel`       | ---                      | ---                                  |
|     peak FW temp (K)     |         `temp_fw_peak`         | ---                      | ---                                  |
|     maximum temp (K)     |       `temp_fw_max`       | ---                      | ---                                  |
|        FCI switch        |           ---           | ---                      | `i_blkt_liquid_breeder_channel_type`                               |
|         Coolant          |      :-----------:      | ------------------------ | --------------------------           |
|  primary coolant switch  |       `i_fw_coolant_type`       | `i_blkt_coolant_type`                 | ---                                  |
| secondary coolant switch |           ---           | ---                      | `i_blkt_liquid_breeder_type`                           |
|      inlet temp (K)      |        `temp_fw_coolant_in`        | `temp_blkt_coolant_in`             | `inlet_temp_liq`                     |
|     outlet temp (K)      |       `temp_fw_coolant_out`        | `temp_blkt_coolant_out`            | `outlet_temp_liq`                    |
|      pressure (Pa)       |      `pres_fw_coolant`       | `pres_blkt_coolant`             | `blpressure_liq`                     |

The default thermo-hydraulic model assumes that a solid breeder is in use, with both the first wall and the breeding blanket using helium as a coolant.
This can be changed using the switches detailed in the following subsection. 


--------------

## First wall

<img
    title="First wall"
    src="../../images/first_wall.png"
    alt="First wall"
    width="90%"
    >

Figure 1: *First wall concept with coolant channels*

The first wall is assumed to be thermally separate from the blanket (Figure 1).  No separation has been made between the structural part of the first wall and the armour.  A simple heuristic model has been used to estimate the peak temperature, as follows.

----------------

### Calculate FW temperature | `fw_temp()`

This function is used to calculate the first wall heating, it assumes the same coolant and geometry parameters for the inboard and outboard first wall though with different possible heat loads.

1. The coolant properties at the inlet and outlet of the first wall are determined using:
    - `i_fw_coolant_type` (coolant type),
    - `temp_fw_coolant_in` & `temp_fw_coolant_out` (inlet and outlet temperatures),
    - `pres_fw_coolant` (fixed coolant pressure).

2. Calculate the average coolant density and heat capacity by averaging the inlet and outlet values.

3. Compute the heat load per unit length $[\text{W/m}]$ of a first wall segment with its own cooling pipe:

    $$
    \text{Load} = \left(\mathtt{pden\_fw\_nuclear} \Delta r_{FW} + \mathtt{pflux\_fw\_rad}\right) \times \mathtt{dx\_fw\_module}
    $$

4. Determine the mean mass flow rate:

    $$
    \dot{m} = \frac{L_{\text{FW}} \times \text{Load}}{c_{\text{average}} \left(T_{\text{out}} - T_{\text{in}}\right)}
    $$

5. Calculate the mass flux in a single channel:

    $$
    \text{Mass flux} = \frac{\dot{m}}{A_{\text{channel}}}
    $$

6. Compute the coolant velocity:

    $$
    \mathtt{vel\_fw\_coolant\_average} = \frac{\text{Mass flux}}{\rho_{\text{outlet}}}
    $$

7. Estimate the mean temperature between the outlet coolant and the peak FW structure temperature:

    $$
    \mathtt{temp_k} = \frac{T_{\text{outlet}} + T_{\text{FW,peak}}}{2}
    $$

8. Calculate the FW thermal conductivity at $\mathtt{temp_k}$ using the [`fw_thermal_conductivity()`](#fw-thermal-conductivity--fw_thermal_conductivity) function.

9. Determine the heat transfer coefficient using the [`heat_transfer()`](#fw-heat-transfer--heat_transfer) function.

10. Compute the worst-case load.

    $$
    \text{Oneload} = \mathtt{f\_fw\_peak} \times \frac{\mathtt{pden\_fw\_nuclear} \times \mathtt{dx\_fw\_module} \times \frac{\Delta r_{FW}}{4}}{\mathtt{pflux\_fw\_rad} \times \times \mathtt{dx\_fw\_module}}
    $$

11. Set the effective heat transfer area equal to the pipe diameter

12. Calculate the temperature drop in the first wall material

    $$
    \Delta T_{\text{FW}} = \frac{\text{Oneload} \times \texttt{dr_fw}}{k 2r_{\text{channel}}}
    $$

13. Maximum distance traveled by surface heat load = $\texttt{diagonal}$

    $$
    \texttt{diagonal}=\sqrt{(\texttt{radius_fw_channel}+\texttt{fw} \_ \texttt{wall})^2 + \left(\frac{\texttt{dx_fw_module}}{2}-\texttt{radius_fw_channel}\right)^2 }
    $$

14. Typical distance travelled by surface heat load:

    $$
    \texttt{mean} \_ \texttt{distance}=\frac{\texttt{fw} \_ \texttt{wall}+\texttt{diagonal}}{2}
    $$

15. 

    $$ 
    \texttt{mean_width} = \frac{\texttt{dx_fw_module} + \pi \times \texttt{radius_fw_channel}}{2}
    $$

------------------

Minimum distance traveled by surface heat load = $\texttt{fw} \_ \texttt{wall}$



The energy travels over a cross-section which is initially $= \texttt{dx_fw_module}$
It spreads out, arriving at the coolant pipe over an area of half the circumference.
We use the mean of these values:



The temperature difference between the plasma-facing surface and the coolant is then:

$$ 
\texttt{deltat} \_ \texttt{solid} = \frac{\texttt{onedload} \times \texttt{mean} \_ \texttt{distance}}{\texttt{tkfw} \times \texttt{mean} \_ \texttt{width}}
$$

where $\texttt{tkfw}$ is the thermal conductivity of the first wall material and $\texttt{onedload}$ is the heat load per unit length.

-------------

### FW heat transfer | `heat_transfer()`

1. **Calculate the Reynolds number:**

    $$
    \mathrm{Re} = \frac{\rho v \left(2r_{\text{channel}}\right)}{\mu}
    $$

    where $\rho$ is the coolant density and $\mu$ is the coolant viscosity.

2. **Calculate the Prandtl number:**

    $$
    \mathrm{Pr} = \frac{c_{\text{p}}\mu}{k}
    $$

    were $c_{\text{p}}$ is the coolant heat capacity and $k$ is the coolant thermal conductivity.

3. **Calculate the Darcy friction factor using the [`darcy_friction_haaland()`](#fw-coolant-friction--darcy_friction_haaland) method:**

    $$
    f = \texttt{darcy_friction_haaland()}
    $$

4. **Calculate the Nusselt number using the [Gnielinski correlation](https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation):**

    $$
    \mathrm{Nu_D}  = \frac{\left(f/8\right)\left(\mathrm{Re}-1000\right)\mathrm{Pr}}{1+12.7\left(f/8\right)^{0.5}\left(\mathrm{Pr}^{2/3}-1\right)}
    $$

    The relation is valid for:

    $$
    0.5 \le \mathrm{Pr} \le 2000 \\
    3000 \le \mathrm{Re} \le 5 \times 10^6
    $$

5. **Calculate the heat transfer coefficient with the Nusselt number:**

    $$
    h = \frac{\mathrm{Nu_D}k}{2r_{\text{channel}}}
    $$


--------------

### FW coolant friction | `darcy_friction_haaland()`

 The pressure drop is based on the Darcy fraction factor, using the [Haaland equation](https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation), an approximation to the implicit Colebrook–White equation. 

$$
\frac{1}{\sqrt{f}} = -1.8 \log{\left[ \left(\frac{\epsilon / D}{3.7}\right)^{1.11} \frac{6.9}{\text{Re}} \right]}
$$

------------

### FW thermal conductivity | `fw_thermal_conductivity()`

The thermal conductivity of the first wall is assumed to be that of Eurofer97 using the relation below[^1] [^2]:
 
$$
K_{\text{Eurofer97}} = 5.4308 + 0.13565T - 0.00023862T^2 + 1.3393 \times 10^{-7} T^3
$$

!!! warning Thermal conductivity validity bounds

    The sources for the stated thermal conductivity relation above state that the relation is only valid up to 800K [^1] [^2].

-------------


### Model Switches

There are three blanket model options, chosen by the user to match their selected blanket design using the switch 'i_blkt_dual_coolant' (default=0):
    0.   Solid breeder - nuclear heating in the blanket is extracted by the primary coolant.
    1.   Liquid metal breeder, single-coolant 
        - nuclear heating in the blanket is extracted by the primary coolant.
        - liquid metal is circulated for tritium extraction, specified by number of circulations/day.
    2.   Liquid metal breeder, dual-coolant
        - nuclear heating in the liquid breeder/coolant is extracted by the liquid breeder/coolant.
        - nuclear heating in the blanket structure is extracted by the primary coolant

The default assuption for all blanket models is that the first wall and breeding blanket have the same coolant (flow = FW inlet -> FW outlet -> BB inlet-> BB outlet). 
It is possible to choose a different coolant for the FW and breeding blanket, in which case the mechanical pumping powers for the FW and BB are calculated separately. 
The model has three mechanical pumping power options, chosen by the user to match their selected blanket design using the switch 'i_fw_blkt_shared_coolant' (default=0): 
    0.   Same coolant for FW and BB ('i_fw_coolant_type`=`i_blkt_coolant_type`)
    1.   Different coolant for FW and BB ('i_fw_coolant_type`/=`i_blkt_coolant_type`) 

!!! Note "Note" 
    For the dual-coolant blanket the 'i_fw_blkt_shared_coolant' switch is relevant for the blanket structure coolant and not the liquid metal breeder/coolant choice.  

The user can select the number poloidal and toroidal modules for the IB and OB BB. The 'i_blkt_module_segmentation' switch can be set to 1 for a single-module-segment blanket (default=0):
    0. Multi-module segment 
    1. Single-module-segment

|   Variable   | Units | Itvar. | Default | Description                                              |
| :----------: | :---: | ------ | ------- | -------------------------------------------------------- |
| `n_blkt_inboard_modules_poloidal` |  ---  |        | 7       | Number of inboard blanket modules in poloidal direction  |
| `n_blkt_outboard_modules_poloidal` |  ---  |        | 8       | Number of outboard blanket modules in poloidal direction |
| `n_blkt_inboard_modules_toroidal` |  ---  |        | 32      | Number of inboard blanket modules in toroidal direction  |
| `n_blkt_outboard_modules_toroidal` |  ---  |        | 48      | Number of outboard blanket modules in toroidal direction |

#### Liquid Breeder or Dual Coolant

There are two material options for the liquid breeder/coolant, chosen by the user to match their selected blanket design using the switch 'i_blkt_liquid_breeder_type' (default=0):
    0.  Lead-Lithium 
    1.  Lithium (needs testing)    
Both options use the mid-temperature of the metal to find the following properties: density, specific heat, thermal conductivity, dynamic viscosity and electrical conductivity. 
The Hartmann number is also calculated (using the magnetic feild strength in the centre of the inboard or outboard blanket module). 

|       Variable        | Units | Scanvar. | Usage         | Default | Description                                           |
| :-------------------: | :---: | -------- | ------------- | ------- | ----------------------------------------------------- |
|   `blpressure_liq`    |  Pa   | 70       | idualcool=1,2 | 1.7D6   | liquid metal breeder/coolant pressure                 |
|   `inlet_temp_liq`    |   K   | 68       | idualcool=1,2 | 570     | Inlet temperatute of liquid metal breeder/coolant     |
|   `outlet_temp_liq`   |   K   | 69       | idualcool=1,2 | 720     | Outlet temperatute of liquid metal breeder/coolant    |
|    `n_liq_recirc`     |  ---  | 71       | idualcool=1   | 10      | Number of liquid metal breeder recirculations per day |
| `f_nuc_pow_bz_struct` |  ---  | 73       | i_blanket_type=5    | 0.34    | FW nuclear power as fraction of total                 |
|  `f_nuc_pow_bz_liq`   |  ---  | 74       | i_blanket_type=5    | 0.66    | Fraction of BZ power cooled by primary coolant        |

#### Flow Channel Inserts for Liquid Metal Breeder

There are three model options, chosen by the user to match their selected blanket design using the switch 'i_blkt_liquid_breeder_channel_type' (default=0):
    0.   No FCIs used. Conductivity of Eurofer steel is assumed for MHD pressure drop calculations in the liquid metal breeder.
    1.   FCIs used, assumed to be perfectly electrically insulating.  
    2.   FCIs used, with conductivity chosen by the user (`bz_channel_conduct_liq`).

|         Variable         |   Units   | Itvar. | Usage       | Default | Description                                                         |
| :----------------------: | :-------: | ------ | ----------- | ------- | ------------------------------------------------------------------- |
| `bz_channel_conduct_liq` | A V-1 m-1 | 72     | ifci = 0, 2 | 8.33D5  | Liquid metal coolant/breeder thin conductor or FCI wall conductance |


[^1]: A. A. Tavassoli et al., “Materials design data for reduced activation martensitic steel type EUROFER,” Journal of Nuclear Materials, vol. 329–333, pp. 257–262, Aug. 2004, doi: https://doi.org/10.1016/j.jnucmat.2004.04.020.

[^2]: Tavassoli, F. "Fusion Demo Interim Structural Design Criteria (DISDC)/Appendix A Material Design Limit Data/A3. S18E Eurofer Steel." CEA, EFDA_TASK_TW4-TTMS-005-D01 (2004).