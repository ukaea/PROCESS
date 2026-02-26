from process.core import process_output as po
from process.data_structure import (
    build_variables,
    first_wall_variables,
    fwbs_variables,
    heat_transport_variables,
    physics_variables,
    stellarator_configuration,
)
from process.data_structure import (
    stellarator_variables as st,
)


def st_build(stellarator, f_output: bool):
    """Routine to determine the build of a stellarator machine

    This routine determines the build of the stellarator machine.
    The values calculated are based on the mean minor radius, etc.,
    as the actual radial and vertical build thicknesses vary with
    toroidal angle.

    Parameters
    ----------
    stellarator :

    f_output:


    """
    if fwbs_variables.blktmodel > 0:
        build_variables.dr_blkt_inboard = (
            build_variables.blbuith + build_variables.blbmith + build_variables.blbpith
        )
        build_variables.dr_blkt_outboard = (
            build_variables.blbuoth + build_variables.blbmoth + build_variables.blbpoth
        )
        build_variables.dz_shld_upper = 0.5e0 * (
            build_variables.dr_shld_inboard + build_variables.dr_shld_outboard
        )

    #  Top/bottom blanket thickness

    build_variables.dz_blkt_upper = 0.5e0 * (
        build_variables.dr_blkt_inboard + build_variables.dr_blkt_outboard
    )

    # First Wall
    build_variables.dr_fw_inboard = (
        2.0e0 * fwbs_variables.radius_fw_channel + 2.0e0 * fwbs_variables.dr_fw_wall
    )
    build_variables.dr_fw_outboard = build_variables.dr_fw_inboard

    build_variables.dr_bore = physics_variables.rmajor - (
        build_variables.dr_cs
        + build_variables.dr_cs_tf_gap
        + build_variables.dr_tf_inboard
        + build_variables.dr_shld_vv_gap_inboard
        + build_variables.dr_vv_inboard
        + build_variables.dr_shld_inboard
        + build_variables.dr_blkt_inboard
        + build_variables.dr_fw_inboard
        + build_variables.dr_fw_plasma_gap_inboard
        + physics_variables.rminor
    )

    #  Radial build to centre of plasma (should be equal to physics_variables.rmajor)
    build_variables.rbld = (
        build_variables.dr_bore
        + build_variables.dr_cs
        + build_variables.dr_cs_tf_gap
        + build_variables.dr_tf_inboard
        + build_variables.dr_shld_vv_gap_inboard
        + build_variables.dr_vv_inboard
        + build_variables.dr_shld_inboard
        + build_variables.dr_blkt_inboard
        + build_variables.dr_fw_inboard
        + build_variables.dr_fw_plasma_gap_inboard
        + physics_variables.rminor
    )

    # Bc stellarators cannot scale physics_variables.rminor reasonably well an additional constraint equation is required,
    # that ensures that there is enough space between coils and plasma.
    build_variables.required_radial_space = (
        build_variables.dr_tf_inboard / 2.0e0
        + build_variables.dr_shld_vv_gap_inboard
        + build_variables.dr_vv_inboard
        + build_variables.dr_shld_inboard
        + build_variables.dr_blkt_inboard
        + build_variables.dr_fw_inboard
        + build_variables.dr_fw_plasma_gap_inboard
    )

    # derivative_min_LCFS_coils_dist  for how strong the stellarator shape changes wrt to aspect ratio
    build_variables.available_radial_space = (
        st.r_coil_minor * st.f_coil_shape - physics_variables.rminor
    )
    +stellarator_configuration.stella_config_derivative_min_lcfs_coils_dist * (
        physics_variables.rminor
        - st.f_st_rmajor * stellarator_configuration.stella_config_rminor_ref
    )
    #  Radius to inner edge of inboard shield
    build_variables.r_shld_inboard_inner = (
        physics_variables.rmajor
        - physics_variables.rminor
        - build_variables.dr_fw_plasma_gap_inboard
        - build_variables.dr_fw_inboard
        - build_variables.dr_blkt_inboard
        - build_variables.dr_shld_inboard
    )

    #  Radius to outer edge of outboard shield
    build_variables.r_shld_outboard_outer = (
        physics_variables.rmajor
        + physics_variables.rminor
        + build_variables.dr_fw_plasma_gap_outboard
        + build_variables.dr_fw_outboard
        + build_variables.dr_blkt_outboard
        + build_variables.dr_shld_outboard
    )

    #  Thickness of outboard TF coil legs
    build_variables.dr_tf_outboard = build_variables.dr_tf_inboard

    #  Radius to centre of outboard TF coil legs

    build_variables.dr_shld_vv_gap_outboard = build_variables.gapomin
    build_variables.r_tf_outboard_mid = (
        build_variables.r_shld_outboard_outer
        + build_variables.dr_vv_outboard
        + build_variables.dr_shld_vv_gap_outboard
        + 0.5e0 * build_variables.dr_tf_outboard
    )

    #  Height to inside edge of TF coil
    #  Roughly equal to average of (inboard build from TF coil to plasma
    #  centre) and (outboard build from plasma centre to TF coil)

    build_variables.z_tf_inside_half = 0.5e0 * (
        (
            build_variables.dr_shld_vv_gap_inboard
            + build_variables.dr_vv_inboard
            + build_variables.dr_shld_inboard
            + build_variables.dr_blkt_inboard
            + build_variables.dr_fw_inboard
            + build_variables.dr_fw_plasma_gap_inboard
            + physics_variables.rminor
        )
        + (
            physics_variables.rminor
            + build_variables.dr_fw_plasma_gap_outboard
            + build_variables.dr_fw_outboard
            + build_variables.dr_blkt_outboard
            + build_variables.dr_shld_outboard
            + build_variables.dr_vv_outboard
            + build_variables.dr_shld_vv_gap_outboard
        )
    )

    #  Outer divertor strike point radius, set equal to major radius

    build_variables.rspo = physics_variables.rmajor

    #  First wall area: scales with minor radius

    # Average minor radius of the first wall
    awall = physics_variables.rminor + 0.5e0 * (
        build_variables.dr_fw_plasma_gap_inboard
        + build_variables.dr_fw_plasma_gap_outboard
    )
    first_wall_variables.a_fw_total = (
        physics_variables.a_plasma_surface * awall / physics_variables.rminor
    )

    if heat_transport_variables.ipowerflow == 0:
        first_wall_variables.a_fw_total = (
            1.0e0 - fwbs_variables.fhole
        ) * first_wall_variables.a_fw_total
    else:
        first_wall_variables.a_fw_total = (
            1.0e0
            - fwbs_variables.fhole
            - fwbs_variables.f_ster_div_single
            - fwbs_variables.f_a_fw_outboard_hcd
        ) * first_wall_variables.a_fw_total

    if f_output:
        #  Print out device build
        output(stellarator)


def output(stellarator):
    po.oheadr(stellarator.outfile, "Radial Build")

    po.ovarre(
        stellarator.outfile,
        "Avail. Space (m)",
        "(available_radial_space)",
        build_variables.available_radial_space,
    )
    po.ovarre(
        stellarator.outfile,
        "Req. Space (m)",
        "(required_radial_space)",
        build_variables.required_radial_space,
    )

    radius = 0.0e0
    po.obuild(stellarator.outfile, "Device centreline", 0.0e0, radius)

    drbild = (
        build_variables.dr_bore + build_variables.dr_cs + build_variables.dr_cs_tf_gap
    )
    radius = radius + drbild
    po.obuild(stellarator.outfile, "Machine dr_bore", drbild, radius, "(dr_bore)")
    po.ovarre(
        stellarator.outfile, "Machine build_variables.dr_bore (m)", "(dr_bore)", drbild
    )

    radius = radius + build_variables.dr_tf_inboard
    po.obuild(
        stellarator.outfile,
        "Coil inboard leg",
        build_variables.dr_tf_inboard,
        radius,
        "(dr_tf_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Coil inboard leg (m)",
        "(deltf)",
        build_variables.dr_tf_inboard,
    )

    radius = radius + build_variables.dr_shld_vv_gap_inboard
    po.obuild(
        stellarator.outfile,
        "Gap",
        build_variables.dr_shld_vv_gap_inboard,
        radius,
        "(dr_shld_vv_gap_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Gap (m)",
        "(dr_shld_vv_gap_inboard)",
        build_variables.dr_shld_vv_gap_inboard,
    )

    radius = radius + build_variables.dr_vv_inboard
    po.obuild(
        stellarator.outfile,
        "Vacuum vessel",
        build_variables.dr_vv_inboard,
        radius,
        "(dr_vv_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Vacuum vessel radial thickness (m)",
        "(dr_vv_inboard)",
        build_variables.dr_vv_inboard,
    )

    radius = radius + build_variables.dr_shld_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard shield",
        build_variables.dr_shld_inboard,
        radius,
        "(dr_shld_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inner radiation shield radial thickness (m)",
        "(dr_shld_inboard)",
        build_variables.dr_shld_inboard,
    )

    radius = radius + build_variables.dr_blkt_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard blanket",
        build_variables.dr_blkt_inboard,
        radius,
        "(dr_blkt_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inboard blanket radial thickness (m)",
        "(dr_blkt_inboard)",
        build_variables.dr_blkt_inboard,
    )

    radius = radius + build_variables.dr_fw_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard first wall",
        build_variables.dr_fw_inboard,
        radius,
        "(dr_fw_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inboard first wall radial thickness (m)",
        "(dr_fw_inboard)",
        build_variables.dr_fw_inboard,
    )

    radius = radius + build_variables.dr_fw_plasma_gap_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard scrape-off",
        build_variables.dr_fw_plasma_gap_inboard,
        radius,
        "(dr_fw_plasma_gap_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inboard scrape-off radial thickness (m)",
        "(dr_fw_plasma_gap_inboard)",
        build_variables.dr_fw_plasma_gap_inboard,
    )

    radius = radius + physics_variables.rminor
    po.obuild(
        stellarator.outfile,
        "Plasma geometric centre",
        physics_variables.rminor,
        radius,
        "(rminor)",
    )

    radius = radius + physics_variables.rminor
    po.obuild(
        stellarator.outfile,
        "Plasma outboard edge",
        physics_variables.rminor,
        radius,
        "(rminor)",
    )

    radius = radius + build_variables.dr_fw_plasma_gap_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard scrape-off",
        build_variables.dr_fw_plasma_gap_outboard,
        radius,
        "(dr_fw_plasma_gap_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outboard scrape-off radial thickness (m)",
        "(dr_fw_plasma_gap_outboard)",
        build_variables.dr_fw_plasma_gap_outboard,
    )

    radius = radius + build_variables.dr_fw_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard first wall",
        build_variables.dr_fw_outboard,
        radius,
        "(dr_fw_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outboard first wall radial thickness (m)",
        "(dr_fw_outboard)",
        build_variables.dr_fw_outboard,
    )

    radius = radius + build_variables.dr_blkt_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard blanket",
        build_variables.dr_blkt_outboard,
        radius,
        "(dr_blkt_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outboard blanket radial thickness (m)",
        "(dr_blkt_outboard)",
        build_variables.dr_blkt_outboard,
    )

    radius = radius + build_variables.dr_shld_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard shield",
        build_variables.dr_shld_outboard,
        radius,
        "(dr_shld_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outer radiation shield radial thickness (m)",
        "(dr_shld_outboard)",
        build_variables.dr_shld_outboard,
    )

    radius = radius + build_variables.dr_vv_outboard
    po.obuild(
        stellarator.outfile,
        "Vacuum vessel",
        build_variables.dr_vv_outboard,
        radius,
        "(dr_vv_outboard)",
    )

    radius = radius + build_variables.dr_shld_vv_gap_outboard
    po.obuild(
        stellarator.outfile,
        "Gap",
        build_variables.dr_shld_vv_gap_outboard,
        radius,
        "(dr_shld_vv_gap_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Gap (m)",
        "(dr_shld_vv_gap_outboard)",
        build_variables.dr_shld_vv_gap_outboard,
    )

    radius = radius + build_variables.dr_tf_outboard
    po.obuild(
        stellarator.outfile,
        "Coil outboard leg",
        build_variables.dr_tf_outboard,
        radius,
        "(dr_tf_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Coil outboard leg radial thickness (m)",
        "(dr_tf_outboard)",
        build_variables.dr_tf_outboard,
    )
