# Style Guide

`PROCESS` uses [`ruff`](https://github.com/astral-sh/ruff) to enforce code formatting and detect common code mistakes/smells. Unless there is a justifiable reason, the `ruff` linter should not be ignored or overridden. 

You can check if your code meets the `ruff` standards by running:

1. `ruff check` and `ruff format` (the later may modify your source code, this is completely safe)
2. `pre-commit run --all-files`

!!! Info "ruff automatic fixes"
    Many problems that `ruff` can detect it can also safely fix. Simply run `ruff check --fix` and it will try and automatically fix the mistakes, rewriting your source code. 


--------------------

## Naming conventions

!!! quote "The Zen of Python"

    *“Explicit is better than implicit.”*

    *"Readability counts."*

### Functions

Use a lowercase word or words. Separate words by underscores(`_`) to improve readability.

---------------------

### Switches

Switches should start with the `i_` prefix in their name, be of integer type and should be indexed from 0.

---------------------

### Constants

Use an uppercase single letter, word, or words. Separate words with underscores to improve readability.

Refrain from declaring or typing known numerical constants directly in the code. Instead call the value from `constants.f90`
If the constants doesn't exist then add it with a source link and uncertainty value.

---------------------

### Variables

!!! note

    Variable names are slowly being converted to their new verbose form as development progresses. Therefore you may come across variables that do not match the style mentioned below.

Use a lowercase single letter, word, or words. Separate words with underscores to improve readability.
If converting between units it may be required to have some capital letters at the end of the variable name to differentiate between different orders of magnitude, `m` and `M` etc.

The agreed upon style is to name the variables by the following scheme:

$\mathtt{var = \ <data type>\_<system>\_<description>\_<units>}$

**Having the units at the end of the name is only necessary when the variable is not in SI/standard units.**

For example, a variable can look like:

$$
\mathtt{var = \ a_\_tf\_wp}
$$

Which represents:

$$
\mathtt{var = \ \text{Area}:\text{Toroidal Field system}:\text{Winding Pack}}
$$

So for above; `a` is the data type representing area, `tf` is the system representing the TF coils and `wp` is the secondary description or system. This variable thus represents the cross-sectional area of the TF coil winding pack.

It may also be useful to use several data type prefixes to greater greater clarity.

For example, a variable can look like:

$$
\mathtt{var = \ f\_a_\_tf\_wp}
$$

In system designation it means this:

$$
\mathtt{var = \ \text{Fraction}:\text{Area}:\text{Toroidal Field system}:\text{Winding Pack}}
$$

This means the variable represents the fraction of the TF coil area taken up by the winding pack.

!!! note "Naming conventions with limit variables"

    For naming variables which represent either upper or lower limits the words `_max` and `_min` should be used in the variable name. Though if a variable represents the highest of a measured value then the variable name should use the word `_peak`.


--------------

#### System designations

Below are a few shorthand designations for different systems that should be used in variable names

- Toroidal field coils: `_tf_`
- Poloidal field coils: `_pf_`
- Vacuum Vessel: `_vv_`
- First wall: `_fw_`
- Divertor: `_div_`
- Blanket: `_blkt_`
- Shield: `_shld_`
- Central Solenoid: `_cs_`
- Heating & Current Drive: `_hcd_`
  - Electron cyclotron current drive: `_eccd_`
  - Ion cyclotron current drive: `_iccd_`
  - Electron Bernstein Wave: `_ebw_`
  - Neutral Beam: `_nb_`
- Centre post: `_cp_` Should only be used for ST's
- Vacuum pumping systems: `_vac_`
- The fusion power plant as a whole: `_plant_`

If the variables are physics variables and do not belong to a system then:

`var = <physics variable>_<description>`

The data types of different variables can be seen below:

---------------------

#### Data types

---------------------

##### Radii and thicknesses

- Radial positions should start with the `r_` prefix
- Radial thicknesses should start with the `dr_` prefix

- Vertical positions should start with the `z_` prefix
- Vertical thicknesses should start with the `dz_` prefix

---------------------

##### Integer countable items

- Integer countable items should start with the `n_` prefix

Example, the total number of TF coils: `n_tf_coils`

---------------------

##### Number densities

- Number density items should start with the `nd_` prefix

---------------------

##### Areas

- Areas should start with the `a_` prefix

Example, the area of the TF winding pack: `a_tf_wp`

---------------------

##### Volumes

- Volumes should start with the `vol_` prefix

---------------------

##### Lengths

- Lengths should start with the `len_` prefix

---------------------

##### Radii

- Radii should start with the `radius_` prefix

---------------------

##### Diameters

- Diameters should start with the `dia_` prefix

---------------------

##### Circumferences

- Circumferences should start with the `circ_` prefix

---------------------

##### Velocities

- Velocities should start with the `vel_` prefix

---------------------

##### Mass

- Masses should start with the `m_` prefix

---------------------

##### Molecule flow rates

- Molecule flow rates should start with the `molflow_` prefix

This should be used for units of $\text{molecules} \cdot \text{s}^{-1}$

This unit is normally used for fuelling or vacuum pumping models

---------------------

##### Volume flow rates

- Volume flow rates should start with the `volflow_` prefix

This should be used for units of $\text{m}^{3} \cdot \text{s}^{-1}$

This unit is normally used for fuelling or vacuum pumping models

---------------------

##### Outgassing rates

- Outgassing rates should start with the `outgrat_` prefix

This should be used for units of $\text{Pa} \cdot \text{m}\text{s}^{-1}$

This is normally used for vacuum systems and outgassing from its surfaces.

---------------------

##### Mass flow rates

- Mass flow rates should start with the `mflow_` prefix

This should be used for units of $\text{kg} \cdot \text{s}^{-1}$

---------------------

##### Mass fluxes rates

- Mass fluxes should start with the `mflux_` prefix

This should be used for units of $\text{kg} \cdot \text{m}^{-2}\text{s}^{-1}$

---------------------

##### Fluid Reynolds numbers

- Reynolds numbers should start with the `reynolds_` prefix

---------------------

##### Darcy Friction Factors

- Darcy Friction Factors for fluid flow numbers should start with the `darcy_frict_` prefix

---------------------

##### Pressures

- Pressures should start with the `pres_` prefix
- Pressure changes or drops should start with the `dpres_` prefix

---------------------

##### Densities

- Densities should start with the `den_` prefix

---------------------

##### Voltages

- Voltages should start with the `v_` prefix

---------------------

##### Resistances

- Resistances should start with the `res_` prefix

---------------------

##### Resistivity

- Resistivity variables should start with the `rho_` prefix

---------------------

##### Residual resistance ratios

- Residual resistance ratios should start with the `rrr_` prefix

---------------------

##### Currents

- Currents should start with the `c_` prefix

---------------------

##### Inductances

- Inductances should start with the `ind_` prefix

---------------------

##### Electric charges

- Electric charges should start with the `charge_` prefix

--------------------

###### Charge numbers

- Electric charge numbers should start with the `n_charge_` prefix.

------------------

##### Current densities

- Current densities should start with the `j_` prefix

---------------------

##### Powers

- Powers should start with the `p_` prefix

---------------------

##### Power densities

- Power densities should start with the `pden_` prefix

---------------------

##### Power fluxes

- Power fluxes should start with the `pflux_` prefix

---------------------

##### Energies

- Energies should start with the `e_` prefix

---------------------

##### Particle fluence

- Fluences should start with the `flu_` prefix

---------------------

##### Temperatures

- Temperatures should start with the `temp_` prefix

---------------------

##### Temperature Changes

- Relative rises or falls in temperatures should start with the `dtemp_` prefix

---------------------

##### Times

- Times should start with the `t_` prefix
- Time intervals should start with the `dt_` prefix

---------------------

##### Magnetic field strengths

- Magnetic field strengths should start with the `b_` prefix

---------------------

##### Magnetic flux

- Magnetic fluxes can start with the `web_` prefix representing Webers.

- Since magnetic flux units are more commonly used in inductive current drive it may be more appropriate
    to use the `vs_` prefix instead representing a $\text{Vs}$.

---------------------

##### Frequencies

- Frequencies should start with the `freq_`

---------------------

##### Angles

- Angles should start with the `deg_` or `rad_` depending on the units used

---------------------

##### Solid Angles

- Solid angles should start with the `ster_` prefix. Short for steradians.

---------------------

##### Lifetimes

- Lifetimes of components should start with the `life_` prefix.

The default units for lifetimes is in years.

The unit declaration `_fpy` can be used to specify that it is the full-power year lifetime.

---------------------

##### Viscosities

- Viscosities should start with the `visc_` prefix.

---------------------

##### Stress

- Stresses should start with the `stress_` prefix followed by the type of stress, for example `stress_shear_`. Normal stresses are assumed in the `stress_` only case.


---------------------

##### Strain

- Strains due to normal stresses should start with the `strain_` prefix.

- Strains due to shear stresses should start with the `strain_shear_` prefix.

---------------------

##### Young's Modulus

- The Young's modulus of materials should start with the `youngmod_` prefix.

---------------------

##### Poisson's ratios

- Poisons ratios for material under stress should start with the `f_poisson_`

--------------------

##### Forces

- Forces should start with the `forc_` prefix.

---------------------

##### Torques

- Torques should start with the `torq_` prefix

---------------------

##### Current drive efficiencies

Absolute current drive efficiencies ($\eta_{\text{CD}}$) representing Amps driven per Watt of injected power start with the `eta_cd` prefix.

$$
\eta_{\text{CD}} = \frac{I_{\text{driven}}}{P_{\text{injected}}}
$$

Normalized current drive efficiecnies using major radius and volume averaged electron temperature start with the `eta_cd_norm` prefix

$$
\eta_{\text{CD,norm}} = R_0 n_{\text{e,20}} \eta_{\text{CD}}
$$

$\eta_{\text{CD,norm}}$ has the units of $\frac{1\times 10^{20} \text{A}}{\text{W} \text{m}^2}$

The above is concurrent with that of general efficiencies given [below](#efficiencies).

--------------


##### Fusion rates

- Fusion rates should start with the `fusrat_` prefix.

This should be used for units of reactions per second.

---------------------

##### Fusion rate densities

- Fusion rate densities should start with the `fusden_` prefix.

This should be used for units of reactions per cubic metre, per second ($\text{reactions} \ \text{m}^{-3}\text{s}^{-1}$).

---------------------

##### Variables representing fractions

If a variable is intended to demonstrate a fraction of a value or distribution etc. Then it should start with the `f_` prefix.

###### Efficiencies

Similar to this is variables representing efficiencies.

If a variable is intended to represent an engineering efficiency then it should start with the `eta_` prefix to represent $\eta$

---------------------


### Variable Length

Try to keep names to a sensible length while also keeping the name explicit and descriptive.

---------------------


### Physical Type

The physical type of the variable should form the first part of the variable name, e.g. for plasma resistance the variable should be named:

```fortran
res_plasma = 1.0
```

Another example would be pulse length

```fortran
time_pulse_length = 7200.0
```

------------------

### Units

Inside PROCESS all variables should be in SI units unless otherwise stated. For example:

```fortran
! Fusion power [W]
p_fusion_total = 1000.0d6
```

If a variable is not in SI units then its units should be put at the end of of the variable name.
Example:

```fortran
! Fusion power [MW]
p_fusion_total_mw = 1000.0d0
```

---------------------

### Coordinates and dimensions

Coordinates should be defined as

```fortran
r_plasma_centre = 9.0d0

z_plasma_centre = 0.0d0

theta_ = 
```

For dimensions

```fortran
dr_cs = 

dz_cs = 

dtheta_description =
```

---------------------

### Loop order

Loop variables that use I, j etc. should use

```fortran
ii
    jj
        kk
            mm
```

---------------------

### Examples

| Variable name | Description | Units |
| ------------- | ----------- | :---: |
| `c_plasma`    | Plasma current | A |
| `c_plasma_MA` | Plasma current | MA |
| `b_t_onaxis`  | Toroidal field on-axis | T |
| `b_t_max`     | Max toroidal field | T |
| `nd_electron_vol` | Volume average electron density | m-3 |
| `temp_electron_vol_eV` | Volume average electron temperature | eV |
| `m_steel` | Mass of steel | kg |
| `m_steel_tonne` | Mass of steel | tonne |
| `e_neutron_eV` | Energy of neutron | eV |
| `e_neutron_MeV` | Energy of neutron | MeV |
| `v_tf_dump` | TF dump voltage | V |
| `t_plant_life` | Plant lifetime | s |
| `t_plant_life_yrs` | Plant lifetime | years |
| `dr_tf_inboard_leg` | TF coil inboard leg radial thickness | m |
| `dr_blanket_inboard` | Inboard blanket thickness | m |
| `vel_coolant` | TF centrepost coolant velocity | m/s |
| `vol_plasma` | Plasma volume | m3 |
| `a_plasma` | Plasma area | m2 |
| `rad_div_target` | Divertor target angle | radians |
| `deg_div_target` | Divertor target angle | deg |
| `s_shear_tf` | TF shear stress  | Pa |
| `` |  |  |

Please see issue [#940](https://github.com/ukaea/PROCESS/issues/940) to discuss new conventions.

## Type-Hints

It is greatly encouraged and recommended to include type hints for all inputs and outputs in Python. Please follow the guidelines set out in [PEP-484](https://peps.python.org/pep-0484/).

## Docstrings

The docstring style is that of the [Sphinx type](https://www.sphinx-doc.org/en/master/index.html). Though there are some additions for `Notes` and `References` in order to give mathematical reasoning and sources to some functions.

### Functions

If writing in new Python functions please use the docstring template below.

```python
def function_name(param1, param2):
    """
    Brief description of what the function does.

    Detailed description of the function. This can include information about the algorithm,
    any important notes, and other relevant details.

    :param type param1: Description of the first parameter.
    :param type param2: Description of the second parameter.
    :returns: Description of the return value.
    :rtype: return_type
    :raises ExceptionType: Description of the exception raised (if any).

    :notes:
        - Additional notes about the function.
        - Any important considerations or caveats.

    :references:
        - Reference 1: Description of the reference.
        - Reference 2: Description of the reference.
    """
```

### Classes

If writing in new Python classes please use the docstring template below.

```python
class ExampleClass:
    """
    Brief description of the class.

    Detailed description of the class. This can include information about the purpose
    of the class, how it should be used, and any other relevant details.

    Attributes:
        attribute1 (type): Description of attribute1.
        attribute2 (type): Description of attribute2.
        attribute3 (type): Description of attribute3.

    Methods:
        method1(param1, param2): Description of method1.
        method2(param1, param2): Description of method2.
    """

    def __init__(self, attribute1, attribute2, attribute3):
        """
        Initializes the ExampleClass with the given attributes.

        :param type attribute1: Description of attribute1.
        :param type attribute2: Description of attribute2.
        :param type attribute3: Description of attribute3.
        """
        self.attribute1 = attribute1
        self.attribute2 = attribute2
        self.attribute3 = attribute3
```

## Comments

- **Comments that contradict the code are worse than no comments**

- Comments should be complete sentences. The first word should be capitalized, unless it is an identifier that begins with a lower case letter.

- Use inline comments sparingly.

- Comments above apply to code below.