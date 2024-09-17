# Style Guide

## Line Length

For optimal readability, a limit of 100 characters for maximum line length has been set. This is 
below the maximum line length of 132 characters for Fortran (to prevent compilation errors) and 
prevents long lines that run on past the edge of the screen wasting programmers time with scrolling.

## Double declarations

PROCESS uses the Fortran 2008+ intrinsic precision module as shown in the example below. The
use statement will need to be at the module level. See the 
[fortran wiki](http://fortranwiki.org/fortran/show/Real+precision) for more information.

```fortran
use, intrinsic :: iso_fortran_env, only: dp=>real64

real(dp) :: b
!! Variable description

```

## Naming conventions


### Case

All variables should be lower case.

### Length

Try to keep variable names to a sensible length. Abbreviations of some parts of the name are suitable e.g. div for divertor. Use underscores to separate words.

### Physical Type

The physical type of the variable should form the first part of the variable name, e.g. for plasma resistance the variable should be named:

```fortran
res_plasma = 1.0
```

Another example would be pulse length

```fortran
time_pulse_length = 7200.0
```

### Units

Inside PROCESS all variables should be in SI units unless otherwise stated. For example:

```fortran
! Fusion power [W]
p_fusion = 1000.0d6

! Fusion power [MW]
p_fusion_mw = 1000.0d0
```

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

### Loop order

Loop variables that use I, j etc. should use

```fortran
ii
    jj
        kk
            mm
```

### Examples

| Variable name | Description | Units |
| ------------- | ----------- | :---: |
| `i_plasma`    | Plasma current | A |
| `i_plasma_ma` | Plasma current | MA |
| `b_t_onaxis`  | Toroidal field on-axis | T |
| `b_t_max`     | Max toroidal field | T |
| `n_electron_vol` | Volume average electron density | m-3 |
| `t_electron_vol_ev` | Volume avgerage electron temperature | eV |
| `m_steel` | Mass of steel | kg |
| `m_steel_tonne` | Mass of steel | tonne |
| `e_neutron_ev` | Energy of neutron | eV |
| `e_neutron_mev` | Energy of neutron | MeV |
| `v_tf_dump` | TF dump voltage | V |
| `time_plant_life` | Plant lifetime | s |
| `time_plant_life_yrs` | Plant lifetime | years |
| `dr_tf_inboard_leg` | TF coil inboard leg radial thickness | m |
| `dr_blanket_inboard` | Inboard blanket thickness | m |
| `velocity_coolant` | TF centrepost coolant velocity | m/s |
| `vol_plasma` | Plasma volume | m3 |
| `a_plasma` | Plasma area | m2 |
| `angle_div_target` | Divertor target angle | radians |
| `angle_div_target_deg` | Divertor target angle | deg |
| `sig_tf_r` | TF radial stress  | Pa |
| `` |  |  |

Please see issue 940 <a href = "https://github.com/ukaea/PROCESS/issues/940"> to discuss new conventions.

# Code Documentation Using FORD

PROCESS uses FORD (FORtran Documentation) to automatically generate documentation from comments 
in the FORTRAN code. FORD parses FORTRAN source to understand the structure of the project as well 
as picking up "docmarked" comments in the source to create the documentation.

Regular Fortran comments are prefixed with a "!"; these are ignored by FORD and don't go into 
the documentation. FORD comments are prefixed by a "!!", called a docmark; these are picked up 
by FORD and go into the documentation.

The "!!" docmark goes after the statement it documents. For example, to document variables:

```fortran
real(kind(1.0D0)) :: alphan = 0.25D0
!! Density profile index

real(kind(1.0D0)) :: alphap = 0.0D0
!! Pressure profile index

real(kind(1.0D0)) :: alpharate = 0.0D0
!! Alpha particle production rate (particles/m3/sec)
```

...and to document modules:
```fortran
module global_variables
  !! Module containing miscellaneous global variables
  !! This module contains miscellaneous global variables not
  !! well-suited to any of the other 'variables' modules.
```

This documentation will appear in the 
[FORD docs](http://process.gitpages.ccfe.ac.uk/process/ford_site/index.html) section in the 
left-hand navigation bar. Within this site, the "Variables" section in the top navigation bar 
provides variable descriptions in the same manner as the original "vardes" page. 

To document a statement before it occurs in the source, use "!>". However, it is encouraged to 
use "!!" for consistency. The rationale behind this and further information is included on the 
[FORD wiki](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation).

The FORD project on github can be found [here](https://github.com/Fortran-FOSS-Programmers/ford).

## Example of FORD documentation for a subroutine (constraint equation)

```fortran

subroutine constraint_eqn_001(args)
  !! author: J Morris
  !! category: equality constraint
  !!
  !! Relationship between beta, temperature (keV) and density
  !!
  !! \begin{equation} 
  !! c_i = 1 - \frac{1}{\beta}\left( \beta_{ft} + \beta_{NBI} + 2 \times 10^3 \mu_0 e
  !! \left( \frac{n_e T_e + n_i T_i}{B_{tot}^2} \right) \right)
  !! \end{equation}
  !!
  !! - \( \beta \) -- total plasma beta
  !! - \( \beta_{ft} \) -- fast alpha beta component
  !! - \( \beta_{NBI} \) -- neutral beam beta component
  !! - \( n_e \) -- electron density [m\(^3\)]
  !! - \( n_i \) -- total ion density [m\(^3\)]
  !! - \( T_e \) -- density weighted average electron temperature [keV]
  !! - \( T_i \) -- density weighted average ion temperature [keV]
  !! - \( B_{tot} \) -- total toroidal + poloidal field [T]

  use physics_variables, only: betaft, betanb, dene, ten, dnitot, tin, btot, beta
  use constants, only: echarge,rmu0

  implicit none

  type(constraint_args_type), intent(out) :: args
  !! constraint derived type

    args%cc = 1.0D0 - (betaft + betanb + &
      2.0D3*rmu0*echarge * (dene*ten + dnitot*tin)/btot**2 )/beta
    args%con = beta * (1.0D0 - args%cc)
    args%err = beta * args%cc
    args%symbol = '='
    args%units  = ''

end subroutine constraint_eqn_001

```

Creates:

<img
    src="../../images/ford_example_1.png"
    alt="alt text"
    width="700px"
    >
