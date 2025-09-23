import numpy as np
size_scaling = 2
current_scaling = 2
bt=3 * current_scaling
c_tf_total=1.3e6*50 * current_scaling
rmajor=5.2 * size_scaling
rminor=0.92 * size_scaling
tdmptf=1
dr_vv_inboard = 0.014
dr_vv_outboard = 0.014
rad_vv_out = rmajor + rminor + 0.5
rad_vv_in = rmajor - rminor - 0.5

# Stellarator version is working on the W7-X scaling, so we should use actual vv r_major
# plasma r_major is just an approximation, but exact calculations require 3D geometry
# Maybe it can be added to the stella_config file in the future
rad_vv = rmajor

# Actual VV force density
# Based on reference values from W-7X:
# Bref = 3;
# Iref = 1.3*50;
# aref = 0.92;
# \[Tau]ref = 1.;
# Rref = 5.2;
# dref = 14*10^-3;

# NOTE: original implementation used taucq which used a EUROfusion
# constant in the calculation. This was the minimum allowed quench time.
# Replacing with the actual quench time.
# MN/m^3
f_vv_actual = (
    2.54e6
    * (3e0 * 1.3e0 * 50e0 * 0.92e0**2e0)
    / (1e0 * 5.2e0 * 0.014e0)
    / (
        bt
        * c_tf_total
        * rminor**2
        / (
            (dr_vv_inboard + dr_vv_outboard)
            / 2
            * tdmptf
            * rad_vv
        )
    )

)
print(f_vv_actual)

f_vv_actual = (
    2.54
    *   (3e0 / bt
        * 1.3e6 * 50e0 / c_tf_total
        * 0.92e0**2e0 / rminor**2
        ) **(-1)
    *   (
            1e0 / tdmptf
            * 5.2e0 / rad_vv
            * 0.014e0 / ((dr_vv_inboard + dr_vv_outboard) / 2)
        )
) 

print(f_vv_actual)

# This is not correct - it gives pressure on the vv wall, not stress
# N/m^2
# is the vv width the correct length to multiply by to turn the
# force density into a stress?
# sctfcoil_module.vv_stress_quench = (
#     f_vv_actual
#     * 1e6
#     * ((dr_vv_inboard + dr_vv_outboard) / 2)
# )

# This approach merge stress model from tokamaks with induced force calculated from W7-X scaling
a_vv = (rad_vv_out + rad_vv_in) / (rad_vv_out - rad_vv_in)
zeta = 1/np.pi + ((a_vv - 1) * np.log((a_vv + 1) / (a_vv - 1)) / (2 * a_vv))
zeta1 = 1/np.pi
zeta2 = ((a_vv - 1) * np.log((a_vv + 1) / (a_vv - 1)) / (2 * a_vv))

print('dump time: ', tdmptf)
print('toroidal stress: ', zeta1 * f_vv_actual  * rad_vv_in)
print('z stress: ', zeta2 * f_vv_actual  * rad_vv_in)
print('total stress: ', zeta * f_vv_actual  * rad_vv_in)

vv_stress_quench =  zeta * f_vv_actual * 1e6 * rad_vv_in