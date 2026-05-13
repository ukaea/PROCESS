from process.core import process_output as po
from process.core.model import DataStructure
from process.data_structure import (
    heat_transport_variables,
    physics_variables,
    stellarator_configuration,
)
from process.data_structure import (
    stellarator_variables as st,
)


def st_build(stellarator, f_output: bool, data: DataStructure):
    """Routine to determine the build of a stellarator machine

    This routine determines the build of the stellarator machine.
    The values calculated are based on the mean minor radius, etc.,
    as the actual radial and vertical build thicknesses vary with
    toroidal angle.

    Parameters
    ----------
    stellarator :

    f_output:

    data: DataStructure
        data structure object to provide model data

    """
    if data.fwbs.blktmodel > 0:
        data.build.dr_blkt_inboard = (
            data.build.blbuith + data.build.blbmith + data.build.blbpith
        )
        data.build.dr_blkt_outboard = (
            data.build.blbuoth + data.build.blbmoth + data.build.blbpoth
        )
        data.build.dz_shld_upper = 0.5e0 * (
            data.build.dr_shld_inboard + data.build.dr_shld_outboard
        )

    #  Top/bottom blanket thickness

    data.build.dz_blkt_upper = 0.5e0 * (
        data.build.dr_blkt_inboard + data.build.dr_blkt_outboard
    )

    # First Wall
    data.build.dr_fw_inboard = (
        2.0e0 * data.fwbs.radius_fw_channel + 2.0e0 * data.fwbs.dr_fw_wall
    )
    data.build.dr_fw_outboard = data.build.dr_fw_inboard

    data.build.dr_bore = physics_variables.rmajor - (
        data.build.dr_cs
        + data.build.dr_cs_tf_gap
        + data.build.dr_tf_inboard
        + data.build.dr_shld_vv_gap_inboard
        + data.build.dr_vv_inboard
        + data.build.dr_shld_inboard
        + data.build.dr_blkt_inboard
        + data.build.dr_fw_inboard
        + data.build.dr_fw_plasma_gap_inboard
        + physics_variables.rminor
    )

    #  Radial build to centre of plasma (should be equal to physics_variables.rmajor)
    data.build.rbld = (
        data.build.dr_bore
        + data.build.dr_cs
        + data.build.dr_cs_tf_gap
        + data.build.dr_tf_inboard
        + data.build.dr_shld_vv_gap_inboard
        + data.build.dr_vv_inboard
        + data.build.dr_shld_inboard
        + data.build.dr_blkt_inboard
        + data.build.dr_fw_inboard
        + data.build.dr_fw_plasma_gap_inboard
        + physics_variables.rminor
    )

    # Bc stellarators cannot scale physics_variables.rminor reasonably well an additional constraint equation is required,
    # that ensures that there is enough space between coils and plasma.
    data.build.required_radial_space = (
        data.build.dr_tf_inboard / 2.0e0
        + data.build.dr_shld_vv_gap_inboard
        + data.build.dr_vv_inboard
        + data.build.dr_shld_inboard
        + data.build.dr_blkt_inboard
        + data.build.dr_fw_inboard
        + data.build.dr_fw_plasma_gap_inboard
    )

    # derivative_min_LCFS_coils_dist  for how strong the stellarator shape changes wrt to aspect ratio
    data.build.available_radial_space = (
        st.r_coil_minor * st.f_coil_shape - physics_variables.rminor
    )
    +stellarator_configuration.stella_config_derivative_min_lcfs_coils_dist * (
        physics_variables.rminor
        - st.f_st_rmajor * stellarator_configuration.stella_config_rminor_ref
    )
    #  Radius to inner edge of inboard shield
    data.build.r_shld_inboard_inner = (
        physics_variables.rmajor
        - physics_variables.rminor
        - data.build.dr_fw_plasma_gap_inboard
        - data.build.dr_fw_inboard
        - data.build.dr_blkt_inboard
        - data.build.dr_shld_inboard
    )

    #  Radius to outer edge of outboard shield
    data.build.r_shld_outboard_outer = (
        physics_variables.rmajor
        + physics_variables.rminor
        + data.build.dr_fw_plasma_gap_outboard
        + data.build.dr_fw_outboard
        + data.build.dr_blkt_outboard
        + data.build.dr_shld_outboard
    )

    #  Thickness of outboard TF coil legs
    data.build.dr_tf_outboard = data.build.dr_tf_inboard

    #  Radius to centre of outboard TF coil legs

    data.build.dr_shld_vv_gap_outboard = data.build.gapomin
    data.build.r_tf_outboard_mid = (
        data.build.r_shld_outboard_outer
        + data.build.dr_vv_outboard
        + data.build.dr_shld_vv_gap_outboard
        + 0.5e0 * data.build.dr_tf_outboard
    )

    #  Height to inside edge of TF coil
    #  Roughly equal to average of (inboard build from TF coil to plasma
    #  centre) and (outboard build from plasma centre to TF coil)

    data.build.z_tf_inside_half = 0.5e0 * (
        (
            data.build.dr_shld_vv_gap_inboard
            + data.build.dr_vv_inboard
            + data.build.dr_shld_inboard
            + data.build.dr_blkt_inboard
            + data.build.dr_fw_inboard
            + data.build.dr_fw_plasma_gap_inboard
            + physics_variables.rminor
        )
        + (
            physics_variables.rminor
            + data.build.dr_fw_plasma_gap_outboard
            + data.build.dr_fw_outboard
            + data.build.dr_blkt_outboard
            + data.build.dr_shld_outboard
            + data.build.dr_vv_outboard
            + data.build.dr_shld_vv_gap_outboard
        )
    )

    #  Outer divertor strike point radius, set equal to major radius

    data.build.rspo = physics_variables.rmajor

    #  First wall area: scales with minor radius

    # Average minor radius of the first wall
    awall = physics_variables.rminor + 0.5e0 * (
        data.build.dr_fw_plasma_gap_inboard + data.build.dr_fw_plasma_gap_outboard
    )
    data.first_wall.a_fw_total = (
        physics_variables.a_plasma_surface * awall / physics_variables.rminor
    )

    if heat_transport_variables.ipowerflow == 0:
        data.first_wall.a_fw_total = (
            1.0e0 - data.fwbs.fhole
        ) * data.first_wall.a_fw_total
    else:
        data.first_wall.a_fw_total = (
            1.0e0
            - data.fwbs.fhole
            - data.fwbs.f_ster_div_single
            - data.fwbs.f_a_fw_outboard_hcd
        ) * data.first_wall.a_fw_total

    if f_output:
        #  Print out device build
        output(stellarator, data)


def output(stellarator, data):
    po.oheadr(stellarator.outfile, "Radial Build")

    po.ovarre(
        stellarator.outfile,
        "Avail. Space (m)",
        "(available_radial_space)",
        data.build.available_radial_space,
    )
    po.ovarre(
        stellarator.outfile,
        "Req. Space (m)",
        "(required_radial_space)",
        data.build.required_radial_space,
    )

    radius = 0.0e0
    po.obuild(stellarator.outfile, "Device centreline", 0.0e0, radius)

    drbild = data.build.dr_bore + data.build.dr_cs + data.build.dr_cs_tf_gap
    radius += drbild
    po.obuild(stellarator.outfile, "Machine dr_bore", drbild, radius, "(dr_bore)")
    po.ovarre(stellarator.outfile, "Machine data.build.dr_bore (m)", "(dr_bore)", drbild)

    radius += data.build.dr_tf_inboard
    po.obuild(
        stellarator.outfile,
        "Coil inboard leg",
        data.build.dr_tf_inboard,
        radius,
        "(dr_tf_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Coil inboard leg (m)",
        "(deltf)",
        data.build.dr_tf_inboard,
    )

    radius += data.build.dr_shld_vv_gap_inboard
    po.obuild(
        stellarator.outfile,
        "Gap",
        data.build.dr_shld_vv_gap_inboard,
        radius,
        "(dr_shld_vv_gap_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Gap (m)",
        "(dr_shld_vv_gap_inboard)",
        data.build.dr_shld_vv_gap_inboard,
    )

    radius += data.build.dr_vv_inboard
    po.obuild(
        stellarator.outfile,
        "Vacuum vessel",
        data.build.dr_vv_inboard,
        radius,
        "(dr_vv_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Vacuum vessel radial thickness (m)",
        "(dr_vv_inboard)",
        data.build.dr_vv_inboard,
    )

    radius += data.build.dr_shld_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard shield",
        data.build.dr_shld_inboard,
        radius,
        "(dr_shld_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inner radiation shield radial thickness (m)",
        "(dr_shld_inboard)",
        data.build.dr_shld_inboard,
    )

    radius += data.build.dr_blkt_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard blanket",
        data.build.dr_blkt_inboard,
        radius,
        "(dr_blkt_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inboard blanket radial thickness (m)",
        "(dr_blkt_inboard)",
        data.build.dr_blkt_inboard,
    )

    radius += data.build.dr_fw_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard first wall",
        data.build.dr_fw_inboard,
        radius,
        "(dr_fw_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inboard first wall radial thickness (m)",
        "(dr_fw_inboard)",
        data.build.dr_fw_inboard,
    )

    radius += data.build.dr_fw_plasma_gap_inboard
    po.obuild(
        stellarator.outfile,
        "Inboard scrape-off",
        data.build.dr_fw_plasma_gap_inboard,
        radius,
        "(dr_fw_plasma_gap_inboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Inboard scrape-off radial thickness (m)",
        "(dr_fw_plasma_gap_inboard)",
        data.build.dr_fw_plasma_gap_inboard,
    )

    radius += physics_variables.rminor
    po.obuild(
        stellarator.outfile,
        "Plasma geometric centre",
        physics_variables.rminor,
        radius,
        "(rminor)",
    )

    radius += physics_variables.rminor
    po.obuild(
        stellarator.outfile,
        "Plasma outboard edge",
        physics_variables.rminor,
        radius,
        "(rminor)",
    )

    radius += data.build.dr_fw_plasma_gap_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard scrape-off",
        data.build.dr_fw_plasma_gap_outboard,
        radius,
        "(dr_fw_plasma_gap_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outboard scrape-off radial thickness (m)",
        "(dr_fw_plasma_gap_outboard)",
        data.build.dr_fw_plasma_gap_outboard,
    )

    radius += data.build.dr_fw_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard first wall",
        data.build.dr_fw_outboard,
        radius,
        "(dr_fw_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outboard first wall radial thickness (m)",
        "(dr_fw_outboard)",
        data.build.dr_fw_outboard,
    )

    radius += data.build.dr_blkt_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard blanket",
        data.build.dr_blkt_outboard,
        radius,
        "(dr_blkt_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outboard blanket radial thickness (m)",
        "(dr_blkt_outboard)",
        data.build.dr_blkt_outboard,
    )

    radius += data.build.dr_shld_outboard
    po.obuild(
        stellarator.outfile,
        "Outboard shield",
        data.build.dr_shld_outboard,
        radius,
        "(dr_shld_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Outer radiation shield radial thickness (m)",
        "(dr_shld_outboard)",
        data.build.dr_shld_outboard,
    )

    radius += data.build.dr_vv_outboard
    po.obuild(
        stellarator.outfile,
        "Vacuum vessel",
        data.build.dr_vv_outboard,
        radius,
        "(dr_vv_outboard)",
    )

    radius += data.build.dr_shld_vv_gap_outboard
    po.obuild(
        stellarator.outfile,
        "Gap",
        data.build.dr_shld_vv_gap_outboard,
        radius,
        "(dr_shld_vv_gap_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Gap (m)",
        "(dr_shld_vv_gap_outboard)",
        data.build.dr_shld_vv_gap_outboard,
    )

    radius += data.build.dr_tf_outboard
    po.obuild(
        stellarator.outfile,
        "Coil outboard leg",
        data.build.dr_tf_outboard,
        radius,
        "(dr_tf_outboard)",
    )
    po.ovarre(
        stellarator.outfile,
        "Coil outboard leg radial thickness (m)",
        "(dr_tf_outboard)",
        data.build.dr_tf_outboard,
    )
