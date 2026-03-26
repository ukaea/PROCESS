# Magnetic Fields | `PlasmaFields()`

In a tokamak there are two main magnetic fields that we are concerned about, the toroidal magnetic field ($B_{\text{T}}$) and the poloidal magnetic field ($B_{\text{P}}$). The former created from the electric current in the toroidal field coils and the latter from the toroidal plasma current.



## Toroidal Field

In `PROCESS` the toroidal magnetic field at the plasma centre  $(B_{\text{T}}(R_0))$ (`b_plasma_toroidal_on_axis`) is normally an iteration variable and is a key paramter is most plasma physics and engineering models.

The toroidal field decreases as $\propto \frac{1}{R}$ from the edge of the inboard toroidal field coil winding pack across the plasma.

------------------------------------

### Plasma Inboard Toroidal Field | `calculate_plasma_inboard_toroidal_field()`

The toroidal field at the plasma inboard is given as:

$$
\overbrace{B_{\text{T}}(R_0-a)}^{\texttt{b_plasma_inboard_toroidal}} = \frac{R_0 B_{\text{T}}(R_0)}{R_0 -a}
$$

------------------------------------

### Plasma Outboard Toroidal Field | `calculate_plasma_outboard_toroidal_field()`

The toroidal field at the plasma outboard is given as:

$$
\overbrace{B_{\text{T}}(R_0+a)}^{\texttt{b_plasma_outboard_toroidal}} = \frac{R_0 B_{\text{T}}(R_0)}{R_0 +a}
$$

------------------------------------

### Plasma Toroidal Field Profile | `calculate_toroidal_field_profile()`

The full toroidal profile across the plasma can be given as:

$$
\overbrace{B_{\text{T}}(\rho)}^{\texttt{b_plasma_toroidal_profile}} = \frac{R_0 B_{\text{T}}(R_0)}{\rho}
$$

------------------------------------

## Poloidal Field

The poloidal field in `PROCESS` is always an output as it is calculated directly from the plasma current. Currently only the average poloidal field on the plasma surface can be calculated. To know the poloidal field within the plasma the Grad-Sharanov equation must be solved, which is beyond the scope of most systems codes. This would provide a correct toroidal current density profile that can be used to find the poloidal field at any point inside the plasma.

------------------------------------

### Plasma Surface Averaged Field | `calculate_surface_averaged_poloidal_field()`

As the total toroidal plasma current is calculated, the Biot-Savart law can be used to find the poloidal field at the plasma surface by using the plasmas calculated poloidal perimeter:

$$
\overbrace{\langle B_{\text{p}} (a) \rangle}^{\texttt{b_plasma_surface_poloidal_average}} = \frac{I_{\text{p}}}{L_{\text{plasma,perimeter}}}
$$

As most plasmas are non ciruclar, the poloidal field thus varies with poloidal angle so only the average value can be inferred from this method. If the plasma was a perfect torus then this would be the poloidal field at any point of the plasma surface

------------------------------------