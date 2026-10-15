import numpy as np
import scipy.optimize
import scipy.special as sp

def from_gw_density_fraction_to_avg_density(current, rminor, gw_density_fraction):
    """由gw密度分数得到平均密度。"""
    gw_density = current / np.pi / rminor**2 * 1e14
    avg_density = gw_density * gw_density_fraction
    return avg_density


def from_line_avg_density_to_axis_density(line_avg_density, alphan):
    """由线平均密度和抛物剖面指数反解轴密度。"""
    return line_avg_density / (
        (sp.gamma(0.5) / 2.0)
        * sp.gamma(alphan + 1.0)
    ) * sp.gamma(alphan + 1.5)


def from_line_avg_and_axis_density_to_alphan(line_avg_density, axis_density):
    """由线平均密度和轴密度反解抛物剖面指数。"""
    target = np.log(axis_density / line_avg_density) + sp.gammaln(0.5) - np.log(2.0)

    def residual(alphan):
        return sp.gammaln(alphan + 1.5) - sp.gammaln(alphan + 1.0) - target

    lower = -1.0 + 1.0e-9
    upper = 1.0
    while residual(upper) < 0.0:
        upper *= 2.0
    return scipy.optimize.brentq(residual, lower, upper)

def from_beta_toroidal_to_alphat(beta_toroidal, ndensity_vol_avg, alphan, temp_on_axis, f_temp, f_ndensity, b_toroidal_on_axis):
    ndensity_on_axis = ndensity_vol_avg * (1.0 + alphan)
    pres_on_axis = ndensity_on_axis * temp_on_axis * 1e3 * 1.602e-19 * (1.0 + f_ndensity * f_temp)
    pres_vol_avg = beta_toroidal * b_toroidal_on_axis**2 / (2.0 * 4.0 * np.pi * 1.0e-7)
    return pres_on_axis / pres_vol_avg - (1.0 + alphan)
    


current = 35e6
rminor = 5.0/1.5
line_avg_density = from_gw_density_fraction_to_avg_density(current, rminor, 2.06)
on_axis_density = 3.1e20
alphan = from_line_avg_and_axis_density_to_alphan(line_avg_density, on_axis_density)
print("alphan = ", alphan)
nd_plasma_electrons_vol_avg = on_axis_density / (1.0 + alphan)
f_nd_fuel_electron = (1.0 - 0.02 * 2.0) / 1.4
f_ndensity = 0.02 + f_nd_fuel_electron
f_temp = 4.0
temp_on_axis = 35.0
beta_toroidal = 0.22
b_toroidal_on_axis = 5.0
alphat = from_beta_toroidal_to_alphat(beta_toroidal, nd_plasma_electrons_vol_avg, alphan, temp_on_axis, f_temp, f_ndensity, b_toroidal_on_axis)
print("ndensity_vol_avg = ", on_axis_density / (1.0 + alphan))
print("alphat = ", alphat)
print("temp_vol_avg = ", temp_on_axis / (1.0 + alphat))
