"""Api functions for PROCESS summary plots."""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Literal

import matplotlib.backends.backend_pdf as bpdf
import matplotlib.pyplot as plt

from process.core.io.mfile import MFile
from process.core.io.plot.summary.common import color_key
from process.core.io.plot.summary.constants import RADIAL_BUILD, vertical_lower
from process.core.io.plot.summary.geometry.build import (
    plot_geometry_info,
    plot_radial_build,
)
from process.core.io.plot.summary.geometry.misc import (
    plot_blkt_pipe_bends,
    plot_blkt_structure,
)
from process.core.io.plot.summary.geometry.poloidal import (
    plot_first_wall_poloidal_cross_section,
    plot_first_wall_top_down_cross_section,
    plot_full_machine_poloidal_cross_section,
    poloidal_cross_section,
)
from process.core.io.plot.summary.geometry.toroidal import toroidal_cross_section
from process.core.io.plot.summary.magnets.cables import (
    plot_cable_in_conduit_cable,
    plot_hts_tape_geometry,
)
from process.core.io.plot.summary.magnets.cs import (
    plot_cs_coil_structure,
    plot_cs_turn_structure,
    plot_magnetics_info,
    plot_pf_cs_plasma_mutual_inductance,
    plot_physics_info,
)
from process.core.io.plot.summary.magnets.pf import plot_pf_dimensions
from process.core.io.plot.summary.magnets.tf import (
    plot_corc_cable_geometry,
    plot_quench_time_evolution,
    plot_resistive_tf_info,
    plot_resistive_tf_wp,
    plot_superconducting_tf_wp,
    plot_tf_cable_in_conduit_turn,
    plot_tf_coil_structure,
    plot_tf_corc_cable_summary_box,
    plot_tf_croco_turn,
    plot_tf_stress,
)
from process.core.io.plot.summary.plasma.confinement import (
    plot_brunner_divertor_power_split_comparison_stackplot,
    plot_confinement_time_comparison,
    plot_sol_power_decay_length_comparison,
)
from process.core.io.plot.summary.plasma.current_drive import plot_bootstrap_comparison
from process.core.io.plot.summary.plasma.overview import (
    plot_detailed_plasma_parameters,
    plot_main_plasma_information,
)
from process.core.io.plot.summary.plasma.physics import (
    plot_magnetic_fields_in_plasma,
    plot_max_normalised_beta_comparison,
    plot_plasma_coloumb_logarithms,
    plot_plasma_current_comparison,
    plot_plasma_outboard_toroidal_ripple_map,
)
from process.core.io.plot.summary.power_flow import (
    plot_main_power_flow,
    plot_power_info,
)
from process.core.io.plot.summary.profiles.atomic import (
    plot_collision_frequency_profile,
    plot_collision_time_profile,
    plot_debye_length_profile,
    plot_electron_frequency_profile,
    plot_ion_charge_profile,
    plot_ion_frequency_profile,
    plot_ion_slowing_down_time_profile,
    plot_mean_free_path_profile,
    plot_resistivity_profile,
    plot_velocity_profile,
)
from process.core.io.plot.summary.profiles.misc import (
    plot_line_brem_loss_function_profile,
    plot_line_brem_power_density_profile,
    plot_line_brem_power_profile,
)
from process.core.io.plot.summary.profiles.plasma import (
    plot_beta_profiles,
    plot_cumulative_plasma_thermal_energy_profiles,
    plot_fusion_rate_contours,
    plot_fusion_rate_profiles,
    plot_jprofile,
    plot_n_profiles,
    plot_plasma_effective_charge_profile,
    plot_plasma_poloidal_pressure_contours,
    plot_plasma_pressure_profiles,
    plot_plasma_thermal_energy_profiles,
    plot_qprofile,
    plot_t_profiles,
)
from process.core.io.plot.summary.profiles.radiation import (
    plot_cs_radial_stress_contour_profile,
    plot_cs_radial_stress_profile,
    plot_larmor_radius_profile,
    plot_plasma_pressure_gradient_profiles,
    plot_rad_density_contour,
)
from process.core.io.plot.summary.profiles.stress import (
    plot_cs_hoop_stress_contour_profile,
    plot_cs_hoop_stress_profile,
    plot_cs_stress_time_profile,
    plot_cs_tresca_2d_contour,
    plot_cs_vertical_stress_profile,
    plot_cs_von_mises_2d_contour,
    plot_vertical_stress_contour_profile,
)
from process.core.io.plot.summary.reporting.constraints import (
    plot_equality_constraint_equations,
    plot_inequality_constraint_equations,
)
from process.core.io.plot.summary.reporting.layouts import plot_upper_vertical_build
from process.core.io.plot.summary.reporting.misc import (
    RadialBuild,
    plot_density_limit_comparison,
    plot_ebw_ecrh_coupling_graph,
    plot_fw_90_deg_pipe_bend,
    plot_h_threshold_comparison,
    plot_iteration_variables,
    plot_lower_vertical_build,
)
from process.core.io.plot.summary.reporting.panels import (
    plot_cover_page,
    plot_header,
    plot_separatrix_power_split,
)
from process.core.io.plot.summary.time_profiles import (
    plot_current_profiles_over_time,
    plot_system_power_profiles_over_time,
)
from process.models.physics.plasma_geometry import PlasmaShapeModelType
from process.models.tfcoil.base import TFConductorModel
from process.models.tfcoil.superconducting import SuperconductingTFTurnType


def main_plot(
    m_file: MFile,
    scan: int,
    imp: str = "../data/lz_non_corona_14_elements/",
    demo_ranges: bool = False,
    colour_scheme: Literal[1, 2] = 1,
) -> list[plt.Figure]:
    """Function to create radial and vertical build plot on given figure.

    Parameters
    ----------
    m_file :
        MFILE
    scan :
        scan to read from MFILE
    imp :
        path to impurity data
    demo_ranges: bool :
         (Default value = False)
    colour_scheme:

    """
    # Checking the impurity data folder
    # Get path to impurity data dir
    # TODO use Path objects throughout module, not strings

    with resources.path(
        "process.data.lz_non_corona_14_elements", "Ar_lz_tau.dat"
    ) as imp_path:
        imp = str(imp_path.parent) + "/"

    i_shape = int(m_file.get("i_plasma_shape", scan=scan))
    # Setup params for text plots
    plt.rcParams.update({"font.size": 8})

    pages = {}

    def _add_page(name: str | None = None):
        """Add a page to the dictionary of pages. If no name is provided, then assign the
        lowest unused number.

        Raises
        ------
        KeyError
            If a page number has already been used
        """
        if name is None:
            prev_index = max((int(k) for k in pages if k.isnumeric()), default=0)
            name = str(prev_index + 1)
        if name in pages:
            raise KeyError(f"Name collision: {name} already in `pages`!")
        pages[name] = plt.figure(figsize=(12, 9), dpi=80)
        return pages[name]

    radial_build = create_thickness_builds(m_file, scan)

    plot_cover_page(
        _add_page("cover").add_subplot(111),
        m_file,
        scan,
        pages["cover"],
        radial_build,
        colour_scheme,
    )

    # Plot header info
    plot_header(_add_page("first").add_subplot(231), m_file, scan)

    # Geometry
    plot_geometry_info(pages["first"].add_subplot(232), m_file, scan)

    # Physics
    plot_physics_info(pages["first"].add_subplot(233), m_file, scan)

    # Magnetics
    plot_magnetics_info(pages["first"].add_subplot(234), m_file, scan)

    # power/flow economics
    plot_power_info(pages["first"].add_subplot(235), m_file, scan)

    # Current drive
    # plot_current_drive_info(pages["first"].add_subplot(236), m_file_data, scan)
    pages["first"].subplots_adjust(wspace=0.25, hspace=0.25)

    ax7 = _add_page().add_subplot(111)
    ax7.set_position([0.25, 0.1, 0.7, 0.8])  # Move plot slightly to the right
    plot_iteration_variables(ax7, m_file, scan)

    ax7_5 = _add_page().add_subplot(313)
    ax7_5.set_position([0.25, 0.1, 0.7, 0.8])
    plot_equality_constraint_equations(ax7_5, m_file, scan)
    ax7_6 = _add_page().add_subplot(111)
    ax7_6.set_position([0.3, 0.1, 0.65, 0.8])
    plot_inequality_constraint_equations(ax7_6, m_file, scan)

    # Plot main plasma information
    plot_main_plasma_information(
        _add_page("plasma_info").add_subplot(111, aspect="equal"),
        m_file,
        scan,
        colour_scheme,
        pages["plasma_info"],
    )

    # Plot density profiles
    plot_n_profiles(_add_page("profiles"), demo_ranges, m_file, scan)

    # Plot temperature profiles
    ax10 = pages["profiles"].add_subplot(232)
    ax10.set_position([0.375, 0.575, 0.25, 0.375])
    plot_t_profiles(ax10, demo_ranges, m_file, scan)

    # Plot impurity profiles
    ax11 = pages["profiles"].add_subplot(233)
    ax11.set_position([0.7, 0.45, 0.25, 0.5])

    plot_line_brem_power_density_profile(
        axis=ax11, mfile=m_file, scan=scan, impp=imp, demo_ranges=demo_ranges
    )

    # Plot current density profile
    ax12 = pages["profiles"].add_subplot(4, 3, 10)
    ax12.set_position([0.375, 0.105, 0.25, 0.15])
    plot_jprofile(ax12, m_file, scan)

    # Plot q profile
    ax13 = pages["profiles"].add_subplot(4, 3, 12)
    ax13.set_position([0.7, 0.105, 0.25, 0.15])
    plot_qprofile(ax13, demo_ranges, m_file, scan)

    ax_line_brem = _add_page("rad_contour").add_subplot(325)
    plot_line_brem_loss_function_profile(
        axis=ax_line_brem,
        mfile=m_file,
        scan=scan,
        impp=imp,
    )

    ax_zeff = pages["rad_contour"].add_subplot(321, sharex=ax_line_brem)
    plot_plasma_effective_charge_profile(ax_zeff, m_file, scan)
    ax_zeff.set_xlabel("")
    ax_zeff.tick_params(
        axis="x", which="both", bottom=True, top=False, labelbottom=False
    )

    ax_ion_charge = pages["rad_contour"].add_subplot(323, sharex=ax_line_brem)
    plot_ion_charge_profile(ax_ion_charge, m_file, scan)
    ax_ion_charge.set_xlabel("")
    ax_ion_charge.tick_params(
        axis="x", which="both", bottom=True, top=False, labelbottom=False
    )

    if i_shape == 1:
        plot_rad_density_contour(
            pages["rad_contour"].add_subplot(122, aspect="equal"), m_file, scan, imp
        )

    if i_shape != 1:
        msg = (
            "Radiation contour plots require a closed (Sauter) plasma"
            " boundary (i_plasma_shape == 1). Current i_plasma_shape ="
            f" {i_shape}. Contour plots are skipped; see the 1D radiation"
            " plots for available information."
        )
        # Add explanatory text to both figures reserved for contour outputs
        pages["rad_contour"].text(
            0.75, 0.5, msg, ha="center", va="center", wrap=True, fontsize=12
        )
    plot_line_brem_power_profile(
        _add_page("line_brem_power").add_subplot(121), m_file, scan, imp
    )
    plot_fusion_rate_profiles(
        _add_page("fusion_rate").add_subplot(122),
        pages["fusion_rate"],
        m_file,
        scan,
    )

    _add_page("rx_1_2"), _add_page("rx_3_4")
    if m_file.get("i_plasma_shape", scan=scan) == PlasmaShapeModelType.SAUTER:
        plot_fusion_rate_contours(pages["rx_1_2"], pages["rx_3_4"], m_file, scan)

    if i_shape != PlasmaShapeModelType.SAUTER:
        msg = (
            "Fusion-rate contour plots require a closed (Sauter) plasma"
            " boundary (i_plasma_shape == 1). Current i_plasma_shape ="
            f" {i_shape}. Contour plots are skipped; see the 1D fusion"
            " rate/profile plots for available information."
        )
        # Add explanatory text to both figures reserved for contour outputs
        pages["rx_1_2"].text(
            0.5, 0.5, msg, ha="center", va="center", wrap=True, fontsize=12
        )
        pages["rx_3_4"].text(
            0.5, 0.5, msg, ha="center", va="center", wrap=True, fontsize=12
        )

    plot_plasma_pressure_profiles(
        _add_page("pressure_profile").add_subplot(222), m_file, scan
    )
    plot_plasma_pressure_gradient_profiles(
        pages["pressure_profile"].add_subplot(224), m_file, scan
    )
    # Currently only works with Sauter geometry as plasma has a closed surface

    if i_shape == PlasmaShapeModelType.SAUTER:
        plot_plasma_poloidal_pressure_contours(
            pages["pressure_profile"].add_subplot(121, aspect="equal"),
            m_file,
            scan,
        )
    else:
        ax = pages["pressure_profile"].add_subplot(131, aspect="equal")
        msg = (
            "Plasma poloidal pressure contours require a closed (Sauter)"
            " plasma boundary (i_plasma_shape =="
            f" {PlasmaShapeModelType.SAUTER}). Current i_plasma_shape ="
            f" {i_shape}. Contour plots are skipped; see the 1D"
            " pressure/profile plots for available information."
        )
        ax.text(
            0.5,
            0.5,
            msg,
            ha="center",
            va="center",
            wrap=True,
            fontsize=10,
            transform=ax.transAxes,
        )
        ax.axis("off")

    plot_magnetic_fields_in_plasma(
        _add_page("beta").add_subplot(122, aspect="equal"), m_file, scan
    )
    plot_beta_profiles(pages["beta"].add_subplot(321), m_file, scan)

    ax_thermal_energy = pages["beta"].add_subplot(325)
    plot_plasma_thermal_energy_profiles(ax_thermal_energy, m_file, scan)
    ax_thermal_energy_cumulative = pages["beta"].add_subplot(
        323, sharex=ax_thermal_energy
    )
    ax_thermal_energy_cumulative.set_position([0.127, 0.35, 0.35, 0.2])
    plot_cumulative_plasma_thermal_energy_profiles(
        ax_thermal_energy_cumulative, m_file, scan
    )

    plot_ebw_ecrh_coupling_graph(_add_page().add_subplot(111), m_file, scan)

    plot_bootstrap_comparison(_add_page("current").add_subplot(221), m_file, scan)
    plot_plasma_current_comparison(pages["current"].add_subplot(224), m_file, scan)
    plot_h_threshold_comparison(
        _add_page("plasma_compare_1").add_subplot(224), m_file, scan
    )
    plot_density_limit_comparison(
        pages["plasma_compare_1"].add_subplot(221), m_file, scan
    )

    plot_max_normalised_beta_comparison(
        _add_page("plasma_compare_2").add_subplot(221), m_file, scan
    )
    plot_confinement_time_comparison(
        pages["plasma_compare_2"].add_subplot(224), m_file, scan
    )

    plot_sol_power_decay_length_comparison(
        _add_page("plasma_compare_3").add_subplot(221), m_file, scan
    )

    plot_brunner_divertor_power_split_comparison_stackplot(
        _add_page("plasma_exhaust").add_subplot(121), m_file, scan
    )

    plot_separatrix_power_split(
        pages["plasma_exhaust"].add_subplot(122), m_file, scan, colour_scheme
    )

    plot_debye_length_profile(
        _add_page("microscopic_quantities").add_subplot(232), m_file, scan
    )
    plot_velocity_profile(pages["microscopic_quantities"].add_subplot(233), m_file, scan)
    plot_plasma_coloumb_logarithms(
        pages["microscopic_quantities"].add_subplot(231), m_file, scan
    )
    plot_collision_time_profile(
        pages["microscopic_quantities"].add_subplot(234), m_file, scan
    )
    plot_collision_frequency_profile(
        pages["microscopic_quantities"].add_subplot(235), m_file, scan
    )
    plot_mean_free_path_profile(
        pages["microscopic_quantities"].add_subplot(236), m_file, scan
    )

    plot_ion_slowing_down_time_profile(
        _add_page("detailed_params").add_subplot(231), m_file, scan
    )

    plot_resistivity_profile(pages["detailed_params"].add_subplot(232), m_file, scan)

    plot_detailed_plasma_parameters(
        pages["detailed_params"].add_subplot(233),
        fig=pages["detailed_params"],
        mfile=m_file,
        scan=scan,
    )

    ax_electron_freq = _add_page("freq").add_subplot(211)
    plot_electron_frequency_profile(ax_electron_freq, m_file, scan)

    ax_ion_freq = pages["freq"].add_subplot(413, sharex=ax_electron_freq)
    plot_ion_frequency_profile(ax_ion_freq, m_file, scan)

    ax_larmor = pages["freq"].add_subplot(414, sharex=ax_electron_freq)
    plot_larmor_radius_profile(ax_larmor, m_file, scan)

    pages["freq"].subplots_adjust(hspace=0.5)

    # Plot poloidal cross-section
    poloidal_cross_section(
        _add_page("tokamak_cross_section").add_subplot(121, aspect="equal"),
        m_file,
        scan,
        demo_ranges,
        radial_build,
        colour_scheme,
    )

    # Plot toroidal cross-section
    toroidal_cross_section(
        pages["tokamak_cross_section"].add_subplot(122, aspect="equal"),
        m_file,
        scan,
        demo_ranges,
        colour_scheme,
    )

    # Plot color key
    ax17 = pages["tokamak_cross_section"].add_subplot(222)
    ax17.set_position([0.5, 0.5, 0.5, 0.5])
    color_key(ax17, m_file, scan, colour_scheme)

    plot_full_machine_poloidal_cross_section(
        _add_page().add_subplot(111, aspect="equal"),
        m_file,
        scan,
        radial_build,
        colour_scheme,
    )

    ax_full_toroidal = _add_page("full_machine_toroidal").add_subplot(
        111, aspect="equal"
    )
    toroidal_cross_section(
        ax_full_toroidal,
        m_file,
        scan,
        demo_ranges,
        colour_scheme,
    )
    ax_full_toroidal.set_ylim(
        -ax_full_toroidal.get_ylim()[1],
        ax_full_toroidal.get_ylim()[1],
    )
    ax_full_toroidal.set_xlim(
        -ax_full_toroidal.get_xlim()[1],
        ax_full_toroidal.get_xlim()[1],
    )

    ax18 = _add_page().add_subplot(211)
    ax18.set_position([0.1, 0.33, 0.8, 0.6])
    plot_radial_build(ax18, m_file, colour_scheme)

    # Make each axes smaller vertically to leave room for the legend
    ax185 = _add_page("vertical_build").add_subplot(211)
    ax185.set_position([0.1, 0.61, 0.8, 0.32])

    ax18b = pages["vertical_build"].add_subplot(212)
    ax18b.set_position([0.1, 0.13, 0.8, 0.32])
    plot_upper_vertical_build(ax185, m_file, colour_scheme)
    plot_lower_vertical_build(ax18b, m_file, colour_scheme)

    # Can only plot WP and turn structure if superconducting coil at the moment
    if m_file.get("i_tf_sup", scan=scan) == TFConductorModel.SUPERCONDUCTING:
        # TF coil with WP
        ax19 = _add_page("tf_wp").add_subplot(221, aspect="equal")
        ax19.set_position([
            0.025,
            0.45,
            0.5,
            0.5,
        ])  # Half height, a bit wider, top left
        plot_superconducting_tf_wp(ax19, m_file, scan, pages["tf_wp"])

        _add_page("cable")
        if (
            m_file.get("i_tf_turn_type", scan=scan)
            == SuperconductingTFTurnType.CROSS_CONDUCTOR
        ):
            ax20 = pages["cable"].add_subplot(325, aspect="equal")
            ax20.set_position([0.025, 0.5, 0.4, 0.4])
            plot_tf_croco_turn(ax20, pages["cable"], m_file, scan)
        elif (
            m_file.get("i_tf_turn_type", scan=scan)
            == SuperconductingTFTurnType.CABLE_IN_CONDUIT
        ):
            # TF coil turn structure
            ax20 = pages["cable"].add_subplot(325, aspect="equal")
            ax20.set_position([0.025, 0.5, 0.4, 0.4])
            plot_tf_cable_in_conduit_turn(ax20, pages["cable"], m_file, scan)

        if (
            m_file.get("i_tf_turn_type", scan=scan)
            == SuperconductingTFTurnType.CROSS_CONDUCTOR
        ):
            plot_205 = pages["cable"].add_subplot(223, aspect="equal")
            plot_205.set_position([0.075, 0.1, 0.3, 0.3])
            plot_corc_cable_geometry(
                plot_205,
                r_centre=0.0,
                z_centre=0.0,
                dia_croco_strand=m_file.get("dia_tf_turn_croco_cable", scan=scan),
                dx_croco_strand_copper=m_file.get(
                    "dx_tf_croco_strand_copper", scan=scan
                ),
                dr_hts_tape=m_file.get("dr_tf_hts_tape", scan=scan),
                dx_croco_strand_tape_stack=m_file.get(
                    "dx_tf_croco_strand_tape_stack", scan=scan
                ),
                n_croco_strand_hts_tapes=m_file.get(
                    "n_tf_croco_strand_hts_tapes", scan=scan
                ),
                dx_hts_tape_rebco=m_file.get("dx_tf_hts_tape_rebco", scan=scan),
                dx_hts_tape_copper=m_file.get("dx_tf_hts_tape_copper", scan=scan),
                dx_hts_tape_hastelloy=m_file.get("dx_tf_hts_tape_hastelloy", scan=scan),
                show_legend=True,
            )
            plot_tf_corc_cable_summary_box(plot_205, pages["cable"], m_file, scan)
            ax_hts_tape = pages["cable"].add_subplot(339)
            ax_hts_tape.set_position([0.75, 0.1, 0.2, 0.2])
            plot_hts_tape_geometry(
                axis=ax_hts_tape,
                r_left=0.0,
                z_bottom=0.0,
                dr_hts_tape=m_file.get("dr_tf_hts_tape", scan=scan),
                dx_hts_tape_rebco=m_file.get("dx_tf_hts_tape_rebco", scan=scan),
                dx_hts_tape_copper=m_file.get("dx_tf_hts_tape_copper", scan=scan),
                dx_hts_tape_hastelloy=m_file.get("dx_tf_hts_tape_hastelloy", scan=scan),
                show_legend=True,
            )
        elif (
            m_file.get("i_tf_turn_type", scan=scan)
            == SuperconductingTFTurnType.CABLE_IN_CONDUIT
        ):
            plot_205 = pages["cable"].add_subplot(223, aspect="equal")
            plot_205.set_position([0.075, 0.1, 0.3, 0.3])
            plot_cable_in_conduit_cable(plot_205, pages["cable"], m_file, scan)
            plot_quench_time_evolution(
                tau_discharge=m_file.get("t_tf_superconductor_quench", scan=scan),
                b_peak=m_file.get("b_tf_inboard_peak_with_ripple", scan=scan),
                f_a_cable_copper=m_file.get("f_a_tf_turn_cable_copper", scan=scan),
                f_a_cable_space_helium=m_file.get(
                    "f_a_tf_turn_cable_space_cooling", scan=scan
                ),
                temp_he_peak=m_file.get("tftmp", scan=scan),
                temp_quench_max=m_file.get("temp_tf_conductor_quench_max", scan=scan),
                cu_rrr=m_file.get("rrr_tf_cu", scan=scan),
                t_quench_detection=m_file.get("t_tf_quench_detection", scan=scan),
                fluence=m_file.get("flu_tf_neutron_fast_max", scan=scan),
                j_operating=m_file.get("j_tf_wp", scan=scan),
                a_tf_turn_cable_space=m_file.get(
                    "a_tf_turn_cable_space_no_void", scan=scan
                ),
                a_tf_turn=m_file.get("a_tf_turn", scan=scan),
                axes_1=_add_page("quench_time_evo").add_subplot(211),
                axes_2=pages["quench_time_evo"].add_subplot(212),
            )
    else:
        ax19 = _add_page("tf_wp").add_subplot(211, aspect="equal")
        ax19.set_position([0.06, 0.55, 0.675, 0.4])
        plot_resistive_tf_wp(ax19, m_file, scan, pages["tf_wp"])
        plot_resistive_tf_info(ax19, m_file, scan, pages["tf_wp"])
    plot_tf_coil_structure(
        _add_page().add_subplot(111, aspect="equal"),
        m_file,
        scan,
        colour_scheme,
    )

    plot_plasma_outboard_toroidal_ripple_map(_add_page(), m_file, scan)

    plot_tf_stress(_add_page().subplots(nrows=3, ncols=1, sharex=True).flatten(), m_file)

    plot_pf_dimensions(
        axis=_add_page("pf_dimensions").add_subplot(121, aspect="equal"),
        mfile=m_file,
        scan=scan,
        colour_scheme=colour_scheme,
    )

    plot_current_profiles_over_time(_add_page().add_subplot(111), m_file, scan)

    plot_pf_cs_plasma_mutual_inductance(_add_page().add_subplot(111), m_file, scan)

    plot_cs_coil_structure(
        _add_page("cs_structure").add_subplot(121, aspect="equal"),
        pages["cs_structure"],
        m_file,
        scan,
    )
    plot_cs_turn_structure(
        pages["cs_structure"].add_subplot(326, aspect="equal"),
        pages["cs_structure"],
        m_file,
        scan,
    )

    plot_cs_stress_time_profile(
        axis=_add_page("cs_stress").add_subplot(337), mfile=m_file, scan=scan
    )

    ax_332 = pages["cs_stress"].add_subplot(332)
    plot_cs_hoop_stress_profile(
        axis=ax_332,
        mfile=m_file,
        scan=scan,
        j_cs=m_file.get("j_cs_pulse_start", scan=scan),
        b_cs_inner=m_file.get("b_cs_peak_pulse_start", scan=scan),
    )

    ax_333 = pages["cs_stress"].add_subplot(333)
    plot_cs_radial_stress_profile(
        axis=ax_333,
        mfile=m_file,
        scan=scan,
        j_cs=m_file.get("j_cs_pulse_start", scan=scan),
        b_cs_inner=m_file.get("b_cs_peak_pulse_start", scan=scan),
    )

    ax_334 = pages["cs_stress"].add_subplot(334)
    ax_334_position = ax_334.get_position()
    cbar_ax_334 = pages["cs_stress"].add_axes([
        ax_334_position.x1 - 0.01,
        ax_334_position.y0,
        0.012,
        ax_334_position.height,
    ])

    ax_336 = pages["cs_stress"].add_subplot(336, sharex=ax_333, sharey=ax_334)
    ax_336_position = ax_336.get_position()
    cbar_ax_336 = pages["cs_stress"].add_axes([
        ax_336_position.x1 + 0.01,
        ax_336_position.y0,
        0.012,
        ax_336_position.height,
    ])

    plot_cs_radial_stress_contour_profile(
        axis=ax_336,
        mfile=m_file,
        scan=scan,
        j_cs=m_file.get("j_cs_pulse_start", scan=scan),
        b_cs_inner=m_file.get("b_cs_peak_pulse_start", scan=scan),
        colorbar_axis=cbar_ax_336,
    )

    ax_331 = pages["cs_stress"].add_subplot(331)
    plot_cs_vertical_stress_profile(
        axis=ax_331,
        mfile=m_file,
        scan=scan,
    )
    plot_vertical_stress_contour_profile(
        axis=ax_334,
        mfile=m_file,
        scan=scan,
        colorbar_axis=cbar_ax_334,
    )

    ax_335 = pages["cs_stress"].add_subplot(335, sharex=ax_332, sharey=ax_334)
    ax_335_position = ax_335.get_position()
    cbar_ax_335 = pages["cs_stress"].add_axes([
        ax_335_position.x1 + 0.01,
        ax_335_position.y0,
        0.012,
        ax_335_position.height,
    ])
    plot_cs_hoop_stress_contour_profile(
        axis=ax_335,
        mfile=m_file,
        scan=scan,
        j_cs=m_file.get("j_cs_pulse_start", scan=scan),
        b_cs_inner=m_file.get("b_cs_peak_pulse_start", scan=scan),
        colorbar_axis=cbar_ax_335,
    )

    pages["cs_stress"].subplots_adjust(wspace=0.45, hspace=0.45)

    # Keep y-axis labeling on the left contour only when sharing y across contour
    # subplots.
    for axis in (ax_335, ax_336):
        axis.set_ylabel("")
        axis.tick_params(axis="y", labelleft=False)

    ax_338 = pages["cs_stress"].add_subplot(338, sharex=ax_332, sharey=ax_335)
    ax_338_position = ax_338.get_position()
    cbar_ax_338 = pages["cs_stress"].add_axes([
        ax_338_position.x1 + 0.01,
        ax_338_position.y0,
        0.012,
        ax_338_position.height,
    ])
    plot_cs_tresca_2d_contour(
        axis=ax_338,
        mfile=m_file,
        scan=scan,
        colorbar_axis=cbar_ax_338,
    )

    ax_339 = pages["cs_stress"].add_subplot(339, sharex=ax_332, sharey=ax_338)
    ax_339_position = ax_339.get_position()
    cbar_ax_339 = pages["cs_stress"].add_axes([
        ax_339_position.x1 + 0.01,
        ax_339_position.y0,
        0.012,
        ax_339_position.height,
    ])

    plot_cs_von_mises_2d_contour(
        axis=ax_339,
        mfile=m_file,
        scan=scan,
        colorbar_axis=cbar_ax_339,
    )

    plot_first_wall_top_down_cross_section(
        _add_page("fw_td_cross_section").add_subplot(221, aspect="equal"),
        m_file,
        scan,
    )
    plot_first_wall_poloidal_cross_section(
        pages["fw_td_cross_section"].add_subplot(122), m_file, scan
    )
    plot_fw_90_deg_pipe_bend(pages["fw_td_cross_section"].add_subplot(337), m_file, scan)

    plot_blkt_pipe_bends(_add_page("blkt_pipe_bends"), m_file, scan)
    ax_blanket = pages["blkt_pipe_bends"].add_subplot(122, aspect="equal")
    plot_blkt_structure(
        ax_blanket,
        pages["blkt_pipe_bends"],
        m_file,
        scan,
        radial_build,
        colour_scheme,
    )

    plot_main_power_flow(
        _add_page("main_power_flow").add_subplot(111, aspect="equal"),
        m_file,
        scan,
        pages["main_power_flow"],
    )

    ax24 = _add_page("power_profile_over_time").add_subplot(111)
    # set_position([left, bottom, width, height]) -> height ~ 0.66 => ~2/3 of page height
    ax24.set_position([0.08, 0.35, 0.84, 0.57])
    plot_system_power_profiles_over_time(
        ax24, m_file, scan, pages["power_profile_over_time"]
    )
    return list(pages.values())


def create_thickness_builds(m_file, scan: int):
    """Create the dictionaries of radial and vertical build values and cumulative
    values
    """
    if int(m_file.get("i_single_null", scan=scan)) == 0:
        vertical_upper = [
            "z_plasma_xpoint_upper",
            "dz_fw_plasma_gap",
            "dz_divertor",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
        ]
    else:
        vertical_upper = [
            "z_plasma_xpoint_upper",
            "dz_fw_plasma_gap",
            "dz_fw_upper",
            "dz_blkt_upper",
            "dr_shld_blkt_gap",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
        ]

    radial = {}
    cumulative_radial = {}
    subtotal = 0
    for item in RADIAL_BUILD:
        if item in {"rminori", "rminoro"}:
            build = m_file.get("rminor", scan=scan)
        elif item in {"vvblgapi", "vvblgapo"}:
            build = m_file.get("dr_shld_blkt_gap", scan=scan)
        elif "dr_vv_inboard" in item:
            build = m_file.get("dr_vv_inboard", scan=scan)
        elif "dr_vv_outboard" in item:
            build = m_file.get("dr_vv_outboard", scan=scan)
        else:
            build = m_file.get(item, scan=scan)

        radial[item] = build
        subtotal += build
        cumulative_radial[item] = subtotal

    upper = {}
    cumulative_upper = {}
    subtotal = 0
    for item in vertical_upper:
        upper[item] = m_file.get(item, scan=scan)
        subtotal += upper[item]
        cumulative_upper[item] = subtotal

    lower = {}
    cumulative_lower = {}
    subtotal = 0
    for item in vertical_lower:
        lower[item] = m_file.get(item, scan=scan)
        subtotal -= lower[item]
        cumulative_lower[item] = subtotal

    return RadialBuild(
        upper,
        lower,
        radial,
        cumulative_upper,
        cumulative_lower,
        cumulative_radial,
    )


def plot_summary(
    mfile: Path,
    scan: int = -1,
    demo_ranges: bool = False,
    colour: Literal[1, 2] = 1,
    output_format: str = "pdf",
    show: bool = False,
):
    """Create the summary.pdf"""

    def add_page_footer(
        fig: plt.Figure, page_number: int, total_pages: int, run_label: str
    ):
        footer_text = f"{run_label}"
        fig.text(
            0.01,
            0.01,
            footer_text,
            fontsize=7,
            ha="left",
            va="bottom",
            color="dimgray",
        )
        fig.text(
            0.99,
            0.01,
            f"Page {page_number}/{total_pages}",
            fontsize=7,
            ha="right",
            va="bottom",
            color="dimgray",
        )

    # create main plot
    # Increase range when adding new page
    # run main_plot
    mfile_obj = MFile(mfile) if mfile else MFile("MFILE.DAT")
    run_label = (
        f"{mfile_obj.get('fileprefix', scan=-1)} | scan {scan or -1} |"
        f" {mfile_obj.get('date', scan=-1)} {mfile_obj.get('time', scan=-1)} |"
        f" {mfile_obj.get('tagno', scan=-1)} | Branch:"
        f" {mfile_obj.get('branch_name', scan=-1)}  "
    )
    pages_of_plots = main_plot(
        mfile_obj,
        scan=scan or -1,
        demo_ranges=demo_ranges,
        colour_scheme=colour,
    )

    if output_format == "pdf":
        with bpdf.PdfPages(mfile.with_name(mfile.name + "SUMMARY.pdf")) as pdf:
            for page_number, p in enumerate(pages_of_plots, start=1):
                add_page_footer(p, page_number, len(pages_of_plots), run_label)
                pdf.savefig(p)
    elif output_format == "png":
        folder = Path(mfile.with_name(mfile.stem + "_SUMMARY"))
        folder.mkdir(parents=True, exist_ok=True)
        for no, page in enumerate(pages_of_plots):
            add_page_footer(page, no + 1, len(pages_of_plots), run_label)
            page.savefig(Path(folder, f"page{no}.png"), format="png")

    # show fig if option used
    if show:
        plt.show(block=True)

    plt.close("all")


__all__ = ["create_thickness_builds", "main_plot", "plot_summary"]
