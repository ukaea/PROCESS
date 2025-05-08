# Plasma Inductance

## Setting the normalised internal inductance

The value of the normalised internal inductance $l_i$ can be set with the `i_ind_plasma_internal_norm` switch.

----------

### User input


The user can specify the value of $l_i$ directly by stating `i_ind_plasma_internal_norm = 0` in the input file.

```python
IN.DAT

i_ind_plasma_internal_norm = 0
ind_plasma_internal_norm = 1.0


```

----------

### Wesson relation

This can be activated by stating `i_ind_plasma_internal_normx = 1` in the input file.

`ind_plasma_internal_norm` is set to `ind_plasma_internal_norm_wesson` using:  

$$
\texttt{ind_plasma_internal_norm_wesson} = \ln{\left(1.65+0.89\alpha_{\text{J}}\right)}
$$

This relation is based off of data taken from DIII-D shots[^1].

This is only recommended for high aspect ratio tokamaks[^2].

**It is recommended to use this switch with [`i_alphaj = 1`](../plasma_current/plasma_current.md#setting-the-current-profile-index) and [`i_beta_norm_max = 1`](../plasma_beta/plasma_beta.md#wesson-relation) as they are self-consistent with each other.**


---------


#### Menard Inductance Relation

This can be activated by stating `ind_plasma_internal_norm = 2` in the input file.

`ind_plasma_internal_norm` is set to `ind_plasma_internal_norm_menard` using[^3]:


$$
\texttt{ind_plasma_internal_norm_menard} = 3.4 - \kappa
$$

This relation is based off of data from NSTX for $l_i$ in the range of 0.4-0.85. This model should be used for $\kappa \ge 2.5$


**This is only recommended for spherical tokamaks**

**It is recommended to use this switch with [`i_beta_norm_max = 3`](../plasma_beta/plasma_beta.md#menard-beta-relation) as they are self-consistent with each other.**

[^1]: T. T. S et al., “Profile Optimization and High Beta Discharges and Stability of High Elongation Plasmas in the DIII-D Tokamak,” Osti.gov, Oct. 1990. https://www.osti.gov/biblio/6194284 (accessed Dec. 19, 2024).

[^2]: Tokamaks 4th Edition, Wesson, page 116

[^3]: J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,” Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016, doi: https://doi.org/10.1088/0029-5515/56/10/106023.

