# Plasma Geometry

The plasma geometric major radius $R_0$ (`rmajor`) and aspect ratio $A$ (`aspect`) 
define the size of the plasma torus. The plasma minor radius $a$ (`rminor`) is 
calculated from these values. The inverse aspect ratio is given by, $\epsilon$ (`eps`) = $1/A$ .

The shape of the plasma cross-section is given by the
elongation of the last closed flux surface (LCFS) $\kappa$ (`kappa`) and the triangularity of the LCFS 
$\delta$ (`triang`), which can be scaled automatically with the aspect ratio if 
required using switch `ishape`:

## Elongation & Triangularity

- `ishape = 0` -- `kappa` and `triang` **must** be input.  The elongation and triangularity of the 95% 
    flux surface are calculated as follows, based on the 1989 ITER guidelines [^1]:
  $$
   \kappa_{95} = \kappa / 1.12
  $$
  $$
   \delta_{95} = \delta / 1.5
  $$
-------------------------------------------------------------------
- `ishape = 1` -- `kappa` and `triang` **must not** be input.  They are calculated by the following equations, 
  which estimate the largest elongation and triangularity achievable for 
  low aspect ratio machines based on the STAR code[^2]:
  
  $$
  \kappa = 2.05 \, \left(1 + 0.44 \, \epsilon^{2.1}\right)
  $$

  $$
  \delta = 0.53 \, \left(1 + 0.77 \, \epsilon^3\right)
  $$

  The values for the plasma shaping parameters at the 95% flux surface are calculated using a fit 
  to a family of equilibria calculated using the FIESTA code, equivalent to that used in `ishape = 8`.

  $$
   \kappa_{95} = \frac{(\kappa - 0.39467)}{0.90698}
  $$

  $$
   \delta_{95} = \frac{(\delta - 0.048306)}{1.3799}
  $$
---------------------------------------------------------------------
- `ishape = 2` -- the Zohm ITER scaling [^3] is used to calculate the elongation, where input variable `fkzohm` $= F_{kz}$ may be used to adjust the scaling, while the input 
  value of the triangularity is used unchanged
  
  $$
  \kappa = F_{kz} \, \times \, \mathrm{minimum} \left( 2.0, \, \, 1.5 + \frac{0.5}{A-1} \right)
  $$

  The elongation and triangularity of the 95% flux surface are calculated as follows, based on the 1989 ITER guidelines [^1]:
  
  $$
   \kappa_{95} = \kappa / 1.12
  $$
  $$
   \delta_{95} = \delta / 1.5
  $$

---------------------------------------------------------------------
- `ishape = 3` -- the Zohm ITER scaling[^3] is used to calculate the elongation (as for `ishape = 2` 
  above), but the triangularity at the 95% flux surface is input via variable `triang95`, and the 
  LCFS triangularity `triang` is calculated from it, rather than the other way round.
---------------------------------------------------------------------  
- `ishape = 4` -- the 95% flux surface values `kappa95` and `triang95` are both used as inputs, 
  and the LCFS values are calculated from them by inverting the equations given above 
  for ``ishape = 0``.
---------------------------------------------------------------------
- `ishape = 5` -- the 95% flux surface values `kappa95` and `triang95` are both used as inputs and 
  the LCFS values are calculated from a fit to MAST data:
  
  $$
  \kappa = 0.913 \, \kappa_{95} + 0.38654
  $$

  $$
  \delta = 0.77394 \, \delta_{95} + 0.18515 
  $$
---------------------------------------------------------------------
- `ishape = 6` -- the input values for `kappa` and `triang` are used directly and the 95% flux 
  surface values are calculated using the MAST scaling from `ishape = 5`.
---------------------------------------------------------------------
- `ishape = 7` -- the 95% flux surface values `kappa95` and `triang95` are both used as inputs and 
  the LCFS values are calculated from a fit to FIESTA runs:
  
  $$
  \kappa = 0.91 \, \kappa_{95} + 0.39
  $$

  $$
  \delta = 1.38 \, \delta_{95} + 0.05 
  $$
---------------------------------------------------------------------
- `ishape = 8` -- the input values for `kappa` and `triang` are used directly and the 95% flux 
  surface values are calculated using the FIESTA fit from `ishape = 7`.
---------------------------------------------------------------------
- `ishape = 9` -- the input values for `triang` and `rli` are used, `kappa` and the 95% flux 
  surface values are calculated.

  $$
   \kappa = \left(\left(1.09+\frac{0.26}{l_i}\right)\left(\frac{1.5}{A}\right)^{0.4}\right)
  $$

  The elongation and triangularity of the 95% flux surface are calculated as follows, based on the 1989 ITER guidelines [^1]:
  
  $$
   \kappa_{95} = \kappa / 1.12
  $$
  $$
   \delta_{95} = \delta / 1.5
  $$
---------------------------------------------------------------------
- `ishape = 10` -- the input values for  `triang` are used directly to calculate 95% flux surface values. `kappa` is calculated to a fit from CREATE data for  a EU-DEMO type machine ($2.6\le A \le 3.6$). 

  $$
  \kappa_{95} = \frac{(-18.84 -(0.87 \times A)) - \sqrt{4.84A^2 -28.77 A+52.51+14.72 m_{s limit}})}{2a}
  $$

  Values rounded to 2 dp

  If $\kappa_{95}>1.77$ then:

  $$
  \kappa_{95} =  \kappa_{95}^{\frac{1.77}{\kappa_{95}}} + \frac{0.3(\kappa_{95}-1.77)}{\frac{1.77}{\kappa_{95}}}
  $$

  The elongation and the triangularity of the 95% flux surface is calculated as follows, based on the 1989 ITER guidelines [^1]:
  
  $$
   \kappa = 1.12\kappa_{95}
  $$
  $$
   \delta_{95} = \delta / 1.5
  $$
  
---------------------------------------------------------------------
- `ishape = 11` -- the elongation is calculated directly dependant on the aspect ratio for spherical tokamak aspect ratios.[^4]

    $$
    \kappa = 0.95 \left(1.9+\frac{1.9}{A^{1.4}}\right)
    $$

  The elongation and triangularity of the 95% flux surface are calculated as follows, based on the 1989 ITER guidelines [^1]:
    
  $$
  \kappa_{95} = \kappa / 1.12
  $$
  $$
  \delta_{95} = \delta / 1.5
  $$

---------------------------------------------------------------------
An explicit constraint relating to the plasma's vertical stability may be turned on if
required. In principle, the inner surface of the outboard shield could be used
as the location of a conducting shell to mitigate the vertical
displacement growth rate of plasmas with significant elongation [^5]. The 
maximum permissible distance $r_{\text{shell, max}}$ of this shell from the geometric 
centre of the plasma may be set using input parameter `cwrmax`, such that 
$r_{\text{shell, max}} =$ `cwrmax*rminor`. Constraint equation 
no. 23 should be turned on with iteration variable no.\ 104 (`fcwr`) to enforce 
this. 

The plasma surface area, cross-sectional area and volume are calculated using
formulations that approximate the LCFS as a revolution of two arcs which
intersect the plasma X-points and the plasma midplane outer and inner
radii. (This is a reasonable assumption for double-null diverted plasmas, but
will be inaccurate for single-null plasmas, `snull = 1`).
 
## Geometrical properties

## Volume

## Surface Area

## Perimeter

## Cross-section

[^1]: N.A. Uckan and ITER Physics Group, *ITER Physics Design Guidelines: 1989*, 
ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)
[^2]: J.D. Galambos, *'STAR Code : Spherical Tokamak Analysis and Reactor Code'*,
Unpublished internal Oak Ridge document.
[^3]: H. Zohm et al, *'On the Physics Guidelines for a Tokamak DEMO'*,
FTP/3-3, Proc. IAEA Fusion Energy Conference, October 2012, San Diego
[^4]: Menard, J.E. & Brown, T. & El-Guebaly, L. & Boyer, M. & Canik, J. & Colling, Bethany & Raman, Roger & Wang, Z. & Zhai, Yunbo & Buxton, Peter & Covele, B. & D’Angelo, C. & Davis, Andrew & Gerhardt, S. & Gryaznevich, M. & Harb, Moataz & Hender, T.C. & Kaye, S. & Kingham, David & Woolley, R.. (2016). *Fusion nuclear science facilities and pilot plants based on the spherical tokamak.* Nuclear Fusion. 56. 106023. 10.1088/0029-5515/56/10/106023. 
[^5]: H.S. Bosch and G.M. Hale, *Improved Formulas for Fusion Cross-sections* 
and Thermal Reactivities', Nuclear Fusion 32 (1992) 611


