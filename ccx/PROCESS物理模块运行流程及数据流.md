# 物理模块运行流程

本文按 `Physics.run`（`process/models/physics/physics.py`）的调用顺序整理。几何量（`rmajor`、`rminor`、`kappa`、`vol_plasma`、`a_plasma_surface` 等）已由前面的 `plasma_geom` 算完。

## 1. `plasma_composition`

`i_fusion_reactions`；

`nd_plasma_electrons_vol_avg`；

燃料：

    `f_plasma_fuel_deuterium + f_plasma_fuel_tritium + f_plasma_fuel_helium3 = 1` 或 `f_plasma_fuel_proton + f_plasma_fuel_boron11 = 1`

氦灰&杂质：

    `f_nd_alpha_thermal_electron`；
    `f_nd_impurity_electrons`;

中性束：

    `f_nd_beam_electron`(仅`i_plasma_ignited = 0`时需要);`f_beam_tritium`

质子：  

    `f_nd_protium_electrons`(仅氘氚氦时需要)；
    `i_nd_plasma_protons`(仅`i_fusion_reactions = "p-b11"`时需要);
    `f_nd_protons_electrons_input`(仅`i_nd_plasma_protons = 1`时需要)

## 2. `clear_veqpy_equilibrium_state`
## 3. `calculate_plasma_masses`
## 4. 库仑对数
## 5. `calculate_plasma_current` and `qstar` and 极向场

`q95`

`b_plasma_toroidal_on_axis`

| switch | 模型，需要用户输入 |
| --- | --- |
| `i_plasma_current = 0` | 用户输入, `plasma_current_user_input` |
| `i_plasma_current = 1` | Peng analytic fit |
| `i_plasma_current = 2` | Peng divertor scaling |
| `i_plasma_current = 3` | Simple ITER scaling（柱近似） |
| `i_plasma_current = 4` | IPDG89 scaling |
| `i_plasma_current = 5` | Todd empirical scaling I |
| `i_plasma_current = 6` | Todd empirical scaling II |
| `i_plasma_current = 7` | Connor-Hastie model |
| `i_plasma_current = 8` | Sauter scaling |
| `i_plasma_current = 9` | FIESTA ST scaling |

## 7. `alphaj`（电流剖面指数）

`q0`

| switch         | 说明，需要用户输入               |
| -------------- | -------------------- |
| `i_alphaj = 0` | `alphaj`     |
| `i_alphaj = 1` | 采用 Wesson 结果 |

## 8. `inductance.run`

| switch         | 说明，需要用户输入               |
| -------------- | -------------------- |
| `i_ind_plasma_internal_norm = 0` | user input, `ind_plasma_internal_norm`   |
| `i_ind_plasma_internal_norm = 1` | wesson scaling |
| `i_ind_plasma_internal_norm = 2` | menard scaling |

## 9. `plasma_profile.run`

| switch | 需要用户输入 |
| --- | --- |
| `i_plasma_pedestal = 0` | 抛物剖面。`alphan`、`alphat`；`nd_plasma_electrons_vol_avg`、`temp_plasma_electron_vol_avg_kev` |
| `i_plasma_pedestal = 1` | 台基剖面。除抛物剖面的输入外，另需 `tbeta`、`radius_plasma_pedestal_density_norm`、`radius_plasma_pedestal_temp_norm`、`temp_plasma_pedestal_electron_kev`、`temp_plasma_separatrix_electron_kev` |
| `i_nd_plasma_pedestal_separatrix = 0`(仅台基剖面) | `nd_plasma_pedestal_electron`、`nd_plasma_separatrix_electron` |
| `i_nd_plasma_pedestal_separatrix = 1`(仅台基剖面) | `f_nd_plasma_pedestal_greenwald`、`f_nd_plasma_separatrix_greenwald` |
| `f_temp_plasma_ion_electron > 0` | `f_temp_plasma_ion_electron` |
| `f_temp_plasma_ion_electron ≤ 0` | `temp_plasma_ion_vol_avg_kev` |
| — | `n_plasma_profile_elements`：剖面网格点数，可选 |

## 10. 平衡（`i_equilibrium_solve`）

| switch                    | 需要用户输入                                               |
| ------------------------- | ---------------------------------------------------- |
| `i_equilibrium_solve = 0` | |
| `i_equilibrium_solve = 1` | |

## 11. 磁场
## 12. `beta.run` and 储能

`beta_total_vol_avg`

`c_beta`

| switch | input |
| --- | --- |
| `i_beta_norm_max = 0` | user input, `beta_norm_max` |
| `i_beta_norm_max = 1` | wesson scaling |
| `i_beta_norm_max = 2` | original scaling |
| `i_beta_norm_max = 3` | menard scaling |
| `i_beta_norm_max = 4` | tholerus scaling |
| `i_beta_norm_max = 5` | stambaugh scaling |


## 13. 脉冲时间

`t_plant_pulse_burn`、 `t_plant_pulse_fusion_ramp` 、 `t_plant_pulse_dwell` 

| switch | 需要用户输入 |
| --- | --- |
| `i_pulsed_plant ≠ 1` 且 `i_t_current_ramp_up = 0` | |
| `i_pulsed_plant = 1` 且 `pulsetimings = 0` | `t_plant_pulse_coil_precharge` |
| `i_pulsed_plant = 1` 且 `pulsetimings ≠ 0` | `t_plant_pulse_plasma_current_ramp_up` 为输入或迭代变量，`t_plant_pulse_coil_precharge` 为爬升和预充磁较大值 |

## 14. `dia_current.run`

| switch                      | 需要用户输入            |
| --------------------------- | ----------------- |
| `i_diamagnetic_current = 0` | 不计算抗磁电流           |
| `i_diamagnetic_current = 1` | TART 定标           |
| `i_diamagnetic_current = 2` | SCENE 定标          |
| `i_equilibrium_solve = 1`   | 抗磁电流份额改由平衡积分，覆盖定标 |




## 15. Pfirsch–Schlüter 电流


| switch                           | 需要用户输入                    |
| -------------------------------- | ------------------------- |
| `i_pfirsch_schluter_current = 0` | 份额保持原值（默认 0）              |
| `i_pfirsch_schluter_current = 1` | 用 SCENE 公式，输入是上一步的总 \beta |




## 16. `plasma_bootstrap_current.run`

自举 + 抗磁 + Pfirsch–Schlüter 合成内部电流份额。超过 `f_c_plasma_non_inductive` 时被截断。感应份额 = 1 - 非感应份额；辅助驱动份额 = 非感应 − 内部。

`i_bootstrap_current`
`f_c_plasma_bootstrap`(`i_bootstrap_current = 0`时)
`f_c_plasma_bootstrap_max`：自举份额上限(非用户输入自举电流时)；`f_c_plasma_non_inductive`：非感应份额；`cboot`：部分定标的系数 

## 17. `current_drive`


| switch                   | 需要用户输入                                     |
| ------------------------ | ------------------------------------------ |
| `i_hcd_calculations = 0` | 不算加热与电流驱动                                  |
| `i_hcd_calculations = 1` | `i_hcd_primary` `i_hcd_secondary` `feffcd` |


`p_hcd_secondary_extra_heat_mw`(`i_hcd_secondary != 0` 时)
`p_hcd_secondary_injected_mw`
`p_hcd_promary_extra_heat_mw`

以 `i_hcd_primary` 为例，hcd模型额外需要的输入


| switch | 模型                               | 需要用户输入                                       |
| ------ | -------------------------------- | -------------------------------------------- |
| 0      | 不计算                              | -                                            |
| 1      | Fenstermacher Lower Hybrid       | -                                            |
| 2      | IPDG89 Ion Cyclotron             | -                                            |
| 3      | Fenstermacher Electron Cyclotron | -                                            |
| 4      | EHST Lower Hybrid                | -                                            |
| 5      | ITER Neutral Beam                | `f_radius_beam_tangency_rmajor` `e_beam_kev` |
| 6      | Culham Lower Hybrid              | -                                            |
| 7      | Culham Electron Cyclotron        | -                                            |
| 8      | Culham Neutral Beam              | `f_radius_beam_tangency_rmajor` `e_beam_kev` |
| 9      | none                             | -                                            |
| 10     | User Input Electron Cyclotron    | `eta_cd_norm_ecrh`                           |
| 11     | none                             | -                                            |
| 12     | User Input Electron Bernstein    | `n_ecrh_harmonic` `xi_ebw`                   |
| 13     | Freethy Electron Cyclotron       | `n_ecrh_harmonic` `i_ecrh_wave_mode`         |


根据模型选择，需要：
`eta_lowhyb_injector_wall_plug` 或
`eta_icrh_injector_wall_plug` 或
`eta_ecrh_injector_wall_plug` 或
`eta_ebw_injector_wall_plug` 或
`f_p_beam_orbit_loss` `eta_beam_injector_wall_plug`

## 18. 聚变反应

| switch | 需要用户输入 |
| --- | --- |
| `i_pb_reactivity_model = 1` | Nevins 模型 |
| `i_pb_reactivity_model = 2` | Tentori–Belloni 模型 |

## 19. 束靶聚变

`beamfus0`、`betbm0` 

没有束或已点火时，总聚变速率等于热等离子体的速率。

## 20. `set_fusion_powers` 与快 α 比压

`f_p_alpha_plasma_deposited`

| switch                         | 需要用户输入                                        |
| ------------------------------ | --------------------------------------------- |
| `i_beta_fast_alpha = 0`        | Uckan 快 α 比压                                  |
| `i_beta_fast_alpha = 1`        | Ward 修正                                       |
| `i_fusion_reactions = "p-b11"` | 上面的结果被 `fast_alpha_beta_pb` 覆盖                |

## 21. 离子–电子能量交换 `rether`
## 22. 辐射

`f_sync_reflect`

`radius_plasma_core_norm`

`f_p_plasma_core_rad_reduction`
| switch                                | 需要用户输入                                                           |
| ------------------------------------- | ---------------------------------------------------------------- |
| `i_calculate_radiation`（仅 p-b11 路径传入） | 选择辐射计算方式                                                         |

## 23. `plasma_ohmic_heating`

`plasma_res_factor`

## 24. `plasma_transition.run`（L–H 阈值）

`i_l_h_threshold` 决定写入 `p_l_h_threshold_mw` 的那一个，供约束 15 使用。

氢硼且阈值开关落在 1,2,3,4,5,15,16,17,18 时，初始化会警告该定标没有 Z_eff 或离子质量修正。

## 25. 分界面功率

P_{\mathrm{sep}} = 沉积的 α 功率 + 非 α 带电粒子功率 + 欧姆功率 − 总辐射。未点火时再加上注入功率。中子不进入。


| switch                 | 需要用户输入                   |
| ---------------------- | ------------------------ |
| `i_plasma_ignited = 1` | 注入功率不进入 P_{\mathrm{sep}} |


算出的 P_{\mathrm{sep}} 会再做一个正值平滑，避免密度上限里出现负功率。

## 26. `scrape_off_layer.run` and 偏滤器 and 电阻扩散时间

`f_len_sol_power_decay_inboard_outboard`

`f_p_div_outboard_separatrix`

| switch                               | 需要用户输入                                                |
| ------------------------------------ | ---------------------------------------------------------- |
| `i_len_sol_outboard_power_decay = 0` | user input, `len_sol_outboard_power_decay`                 |
| `= 1`                                | Eich 2013                                                  |
| `= 2 / 3`                            | MAST 2014 的两条                                           |
| `= 4 / 5`                            | Eich 2011 JET / JET+ASDEX                                  |
| `n_divertors = 2`                    | `f_p_div_lower_separatrix`：下偏滤器分到的功率份额           |

## 27. `density_limit.run`

需要用户输入：`i_density_limit`， `prn1`(边界密度与体平均密度之比，来自偏滤器数据)

八个上限都算。只有被选中的那个写入 `nd_plasma_electrons_max`，给约束方程 5。Greenwald 那一项始终另存，用来算 Greenwald 份额。


| switch                | 说明                                          |
| --------------------- | ------------------------------------------- |
| `i_density_limit = 1` | 旧 ASDEX。用 P_{\mathrm{sep}}、磁场、q_{95}、`prn1` |
| `= 2 / 3`             | Borrass ITER I / II                         |
| `= 4`                 | JET 边界辐射。用 Z_{\mathrm{eff}}、注入功率、q_*        |
| `= 5`                 | JET 简化                                      |
| `= 6`                 | Hugill–Murakami                             |
| `= 7`                 | Greenwald。约束 5 比较线平均密度                      |
| `= 8`（默认）             | 新 ASDEX。用注入功率、电流、q_{95}、`prn1`              |




## 28. `calculate_confinement_time` 与输运损失

需要用户输入：`hfact` `i_confinement_time` `i_plasma_ignited` `f_p_alpha_plasma_deposited` `i_rad_loss`


| switch                   | 需要用户输入                                       |
| ------------------------ | -------------------------------------------- |
| `i_confinement_time = 0` | `tauee_in`                                   |
| `= 34`（默认）               | IPB98(y,2)。质量用燃料质量                  |
| `= 41`                   | Petty 2008，几乎不随 \beta 变，无质量项                 |
| `= 46 / 47 / 48`         | Menard NSTX / 与 Petty 的混合 / Buxton gyro-Bohm |
| `= 53`                   | Kurskiev 球形环定标，无质量项                          |
| `i_rad_loss = 0`         | 损失功率减去全部辐射                                   |
| `i_rad_loss = 1`（默认）     | 只减去芯部辐射                                      |
| `i_rad_loss = 2`         | 不减辐射                                         |
| `i_plasma_ignited = 1`   | 注入功率不进入损失功率                                  |
| —                        | `hfact`：乘在选出的 \tau_E 上                       |


恒星器定标（21–23、37–38）不能在托卡马克上用。

## 29. 伏秒、燃料燃耗份额、加料率与三乘积


| switch | 需要用户输入                                           |
| ------ | ------------------------------------------------ |
| —      | `ejima_coeff`、`csawth`：伏秒。燃烧时间和爬升时间来自第 13 步      |
| —      | `burnup_in`：≤ 10^{-9} 时由 α 灰和燃料密度算燃耗，否则直接用该输入    |
| —      | `tauratio`：氦与加料粒子约束时间之比，`burnup_in`：≤ 10^{-9}时需要 |

## 30. 偏滤器功率分配与辐射份额


| switch | 需要用户输入                                                        |
| ------ | ------------------------------------------------------------- |
| —      | `f_p_div_lower_separatrix` 先换成中平面分离面偏移，再按 Brunner 拆成上下、内外四路功率 |
| —      | `rad_fraction_sol`：刮削层辐射份额，和芯部辐射份额合成总辐射份额                     |




## 31. Reinke 分界面温度


| switch        | 需要用户输入       |
| ------------- | ---------------- |
| 约束 `icc` 含 78 | 用 P_{\mathrm{sep}}/P_{LH}、Greenwald 份额、`q95`、`kappa` 重写 `temp_plasma_separatrix_electron_kev`。另需 `reinke.lhat`、`reinke.impvardiv` 以及约束：`reinke.fzmin` `reinke.fzactual` |
| 未打开 78        | 跳过。分界面温度保持剖面输入     |


# 注意事项

1. `beta_thermal_vol_avg` 是根据总 beta 和 束 beta 、 fast alpha bata 计算得到的，并非真正的热比压

2. HCD 功率 `p_hcd_injected_total_mw` 是通过辅助加热电流和电流驱动效率得到的，不是从能量平衡方程得到