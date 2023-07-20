#! /usr/bin/env python

"""

  Program to output an MCNP compatible file from MFILE.DAT

  James Morris

  13/06.2014

  Notes:
    + 13/06/2014 Original version created.
    + 01/07/2014 Added option for CTF
    + 09/09/2014 Added extra TF coil information
"""

import argparse
from process.io.mfile import MFile


RADIAL_BUILD = [
    "bore",
    "ohcth",
    "gapoh",
    "tfcth",
    "gapds",
    "d_vv_in",
    "shldith",
    "blnkith",
    "fwith",
    "scrapli",
    "rminor",
    "scraplo",
    "fwoth",
    "blnkoth",
    "shldoth",
    "gapsto",
    "tfthko",
]


class ProcessEllipse(object):
    """A class to store ellipse data"""

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        return


class ProcessCylinder(object):
    """A class to store ellipse data"""

    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy
        return


class ProcessPlane(object):
    """A class to store ellipse data"""

    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy
        return


class ProcessTFArcs(object):
    """A class to store TF arc data"""

    def __init__(self, c_x, c_y, x_1, y_1, x_2, y_2):
        self.c_x = c_x
        self.c_y = c_y
        self.x_1 = x_1
        self.y_1 = y_1
        self.x_2 = x_2
        self.y_2 = y_2
        return


class InfoHolder(object):
    """A Class to store information"""

    def __init__(self, name, value):
        self.name = name
        self.value = value
        return


def main(cl_args):
    """Main function for converting MFILE data into MCNP format

    Arguments:
      cl_args --> command line arguments

    """

    mfile_data = MFile(cl_args.f)
    shapes = dict()

    if cl_args.ctf:
        populate_ctf_ellipse_data(shapes, mfile_data)
    else:
        populate_tok_ellipse_data(shapes, mfile_data)

    write_shapes_to_file(shapes, cl_args.o)


def get_xyz(mfile_obj):
    """Function to get x,y,z data for all ellipses"""
    x = mfile_obj.data["rmajor"].get_scan(-1) - mfile_obj.data["rminor"].get_scan(
        -1
    ) * mfile_obj.data["triang95"].get_scan(-1)
    y = 0.0
    z = 0.0
    return x, y, z


def populate_ctf_ellipse_data(shape_objs, mf_data):
    """Function to populate the ellipse objects with data for CTF case"""

    r_major = mf_data.data["rmajor"].get_scan(-1)
    r_minor = mf_data.data["rminor"].get_scan(-1)
    kappa_95 = mf_data.data["kappa95"].get_scan(-1)
    delta_95 = mf_data.data["triang95"].get_scan(-1)
    h_plasma = kappa_95 * r_minor
    t_tf = mf_data.data["tfcth"].get_scan(-1)
    t_scrape_o = mf_data.data["scraplo"].get_scan(-1)
    t_fw_i = mf_data.data["fwith"].get_scan(-1)
    t_fw_o = mf_data.data["fwoth"].get_scan(-1)
    t_shield_i = mf_data.data["shldith"].get_scan(-1)
    t_shield_o = mf_data.data["shldoth"].get_scan(-1)
    t_blanket_o = mf_data.data["blnkoth"].get_scan(-1)
    r_c = r_major - (r_minor * delta_95)

    # TF centre column
    thickness_tf = t_tf
    height_tf = h_plasma + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    shape_objs["c1"] = ProcessCylinder(thickness_tf, height_tf)

    # Inboard shield
    thickness_shield = t_shield_i
    height_shield = h_plasma + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    shape_objs["c2"] = ProcessCylinder(thickness_shield + thickness_tf, height_shield)

    # Inboard first wall
    thickness_fw = t_fw_i
    height_fw = h_plasma + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    shape_objs["c3"] = ProcessCylinder(
        thickness_fw + thickness_tf + thickness_shield, height_fw
    )

    # Upper shield mid-section
    horizontal_thickness_shield_u = r_c - (t_tf + t_shield_i + t_fw_i)
    z_upper_shield_u = h_plasma + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    r_outer_shield_u = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p1"] = ProcessPlane(r_outer_shield_u, z_upper_shield_u)

    # Upper blanket mid-section
    z_upper_blanket_u = h_plasma + t_scrape_o + t_fw_o + t_blanket_o
    r_outer_blanket_u = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p2"] = ProcessPlane(r_outer_blanket_u, z_upper_blanket_u)

    # Upper first wall mid-section
    z_upper_fw_u = h_plasma + t_scrape_o + t_fw_o
    z_lower_fw_u = h_plasma + t_scrape_o
    r_outer_fw_u = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p3"] = ProcessPlane(r_outer_fw_u, z_upper_fw_u)

    # Upper first wall mid-section lower edge
    z_upper_fw_u = h_plasma + t_scrape_o + t_fw_o
    z_lower_fw_u = h_plasma + t_scrape_o
    r_outer_fw_u = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p4"] = ProcessPlane(r_outer_fw_u, z_lower_fw_u)

    # Lower first wall mid-section
    z_upper_fw_l = -h_plasma - t_scrape_o
    z_lower_fw_l = -h_plasma - t_scrape_o - t_fw_o
    r_outer_fw_l = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p4"] = ProcessPlane(r_outer_fw_l, z_upper_fw_l)

    # Lower first wall mid-section
    z_upper_fw_l = -h_plasma - t_scrape_o
    z_lower_fw_l = -h_plasma - t_scrape_o - t_fw_o
    r_outer_fw_l = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p5"] = ProcessPlane(r_outer_fw_l, z_lower_fw_l)

    # Lower blanket mid-section
    z_lower_blanket_l = -h_plasma - t_scrape_o - t_fw_o - t_blanket_o
    r_outer_blanket_l = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p6"] = ProcessPlane(r_outer_blanket_l, z_lower_blanket_l)

    # Lower shield mid-section
    z_lower_shield_l = -h_plasma - t_scrape_o - t_fw_o - t_blanket_o - t_shield_o
    r_outer_shield_l = t_tf + t_shield_i + t_fw_i + horizontal_thickness_shield_u
    shape_objs["p7"] = ProcessPlane(r_outer_shield_l, z_lower_shield_l)

    # Ellipse info
    x, y, z = get_xyz(mf_data)

    # Ellipse 1: outboard first wall inner edge
    shape_objs["e1"] = ProcessEllipse(x, y, z)
    shape_objs["e1"].a = (r_major - r_c) + r_minor + t_scrape_o
    shape_objs["e1"].b = h_plasma + t_scrape_o
    shape_objs["e1"].c = r_c

    # Ellipse 2: outboard blanket inner edge
    shape_objs["e2"] = ProcessEllipse(x, y, z)
    shape_objs["e2"].a = (r_major - r_c) + r_minor + t_scrape_o + t_fw_o
    shape_objs["e2"].b = h_plasma + t_scrape_o + t_fw_o
    shape_objs["e2"].c = r_c

    # Ellipse 3: outboard shield inner edge
    shape_objs["e3"] = ProcessEllipse(x, y, z)
    shape_objs["e3"].a = (r_major - r_c) + r_minor + t_scrape_o + t_fw_o + t_blanket_o
    shape_objs["e3"].b = h_plasma + t_scrape_o + t_fw_o + t_blanket_o
    shape_objs["e3"].c = r_c

    # Ellipse 4: outboard shield outer edge
    shape_objs["e4"] = ProcessEllipse(x, y, z)
    shape_objs["e4"].a = (
        (r_major - r_c) + r_minor + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    )
    shape_objs["e4"].b = h_plasma + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    shape_objs["e4"].c = r_c

    # Extra information
    # TF coil inner leg thickness
    name = "TF Coil inner leg thickness (m)"
    value = t_tf
    shape_objs["i1"] = InfoHolder(name, value)

    # TF coil inner leg thickness
    name = "TF Coil outer leg thickness (m)"
    value = mf_data.data["tfthko"].get_scan(-1)
    shape_objs["i2"] = InfoHolder(name, value)

    # TF coil inner plasma facing case thickness
    name = "TF Coil inner leg plasma facing case thickness (m)"
    value = mf_data.data["casthi"].get_scan(-1)
    shape_objs["i3"] = InfoHolder(name, value)

    # TF coil inner case thickness
    name = "TF Coil inner case thickness (m)"
    value = mf_data.data["thkcas"].get_scan(-1)
    shape_objs["i4"] = InfoHolder(name, value)

    # TF coil winding pack thickness
    name = "TF Coil winding pack thickness (m)"
    value = mf_data.data["dr_tf_wp"].get_scan(-1)
    shape_objs["i5"] = InfoHolder(name, value)

    return shape_objs


def populate_tok_ellipse_data(shape_objs, mf_data):
    """Function to populate the ellipse objects with data for CTF case"""

    r_major = mf_data.data["rmajor"].get_scan(-1)
    r_minor = mf_data.data["rminor"].get_scan(-1)
    kappa_95 = mf_data.data["kappa95"].get_scan(-1)
    delta_95 = mf_data.data["triang95"].get_scan(-1)
    h_plasma = kappa_95 * r_minor
    t_tf = mf_data.data["tfcth"].get_scan(-1)
    t_scrape_i = mf_data.data["scrapli"].get_scan(-1)
    t_scrape_o = mf_data.data["scraplo"].get_scan(-1)
    t_fw_i = mf_data.data["fwith"].get_scan(-1)
    t_fw_o = mf_data.data["fwoth"].get_scan(-1)
    t_shield_i = mf_data.data["shldith"].get_scan(-1)
    t_shield_o = mf_data.data["shldoth"].get_scan(-1)
    t_blanket_i = mf_data.data["blnkith"].get_scan(-1)
    t_blanket_o = mf_data.data["blnkoth"].get_scan(-1)
    r_c = r_major - (r_minor * delta_95)

    # Ellipse info
    x, y, z = get_xyz(mf_data)

    # Ellipse 1: inboard shield outer edge
    shape_objs["e1"] = ProcessEllipse(x, y, z)
    shape_objs["e1"].a = (
        r_minor - (r_major - r_c) + t_scrape_i + t_fw_i + t_blanket_i + t_shield_i
    )
    shape_objs["e1"].b = h_plasma + t_scrape_o + t_fw_i + t_blanket_i + t_shield_i
    shape_objs["e1"].c = r_c

    # Ellipse 2: inboard blanket outer edge
    shape_objs["e2"] = ProcessEllipse(x, y, z)
    shape_objs["e2"].a = r_minor - (r_major - r_c) + t_scrape_i + t_fw_i + t_blanket_i
    shape_objs["e2"].b = h_plasma + t_scrape_i + t_fw_i + t_blanket_i
    shape_objs["e2"].c = r_c

    # Ellipse 3: inboard first wall outer edge
    shape_objs["e3"] = ProcessEllipse(x, y, z)
    shape_objs["e3"].a = r_minor - (r_major - r_c) + t_scrape_i + t_fw_i + t_blanket_o
    shape_objs["e3"].b = h_plasma + t_scrape_i + t_fw_i
    shape_objs["e3"].c = r_c

    # Ellipse 4: inboard first wall inner edge
    shape_objs["e4"] = ProcessEllipse(x, y, z)
    shape_objs["e4"].a = r_minor - (r_major - r_c) + t_scrape_i
    shape_objs["e4"].b = h_plasma + t_scrape_i
    shape_objs["e4"].c = r_c

    # Ellipse 5: outboard first wall inner edge
    shape_objs["e5"] = ProcessEllipse(x, y, z)
    shape_objs["e5"].a = (r_major - r_c) + r_minor + t_scrape_o
    shape_objs["e5"].b = h_plasma + t_scrape_o
    shape_objs["e5"].c = r_c

    # Ellipse 6: outboard blanket inner edge
    shape_objs["e6"] = ProcessEllipse(x, y, z)
    shape_objs["e6"].a = (r_major - r_c) + r_minor + t_scrape_o + t_fw_o
    shape_objs["e6"].b = h_plasma + t_scrape_o + t_fw_o
    shape_objs["e6"].c = r_c

    # Ellipse 7: outboard shield inner edge
    shape_objs["e7"] = ProcessEllipse(x, y, z)
    shape_objs["e7"].a = (r_major - r_c) + r_minor + t_scrape_o + t_fw_o + t_blanket_o
    shape_objs["e7"].b = h_plasma + t_scrape_o + t_fw_o + t_blanket_o
    shape_objs["e7"].c = r_c

    # Ellipse 8: outboard shield outer edge
    shape_objs["e8"] = ProcessEllipse(x, y, z)
    shape_objs["e8"].a = (
        (r_major - r_c) + r_minor + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    )
    shape_objs["e8"].b = h_plasma + t_scrape_o + t_fw_o + t_blanket_o + t_shield_o
    shape_objs["e8"].c = r_c

    # TF arc 1
    c_x = mf_data.data["xctfc(1)"].get_scan(-1)
    c_y = mf_data.data["yctfc(1)"].get_scan(-1)
    x_1 = mf_data.data["xarc(1)"].get_scan(-1)
    y_1 = mf_data.data["yarc(1)"].get_scan(-1)
    x_2 = mf_data.data["xarc(2)"].get_scan(-1)
    y_2 = mf_data.data["yarc(2)"].get_scan(-1)
    shape_objs["a1"] = ProcessTFArcs(c_x, c_y, x_1, y_1, x_2, y_2)

    # TF arc 2
    c_x = mf_data.data["xctfc(2)"].get_scan(-1)
    c_y = mf_data.data["yctfc(2)"].get_scan(-1)
    x_1 = mf_data.data["xarc(2)"].get_scan(-1)
    y_1 = mf_data.data["yarc(2)"].get_scan(-1)
    x_2 = mf_data.data["xarc(3)"].get_scan(-1)
    y_2 = mf_data.data["yarc(3)"].get_scan(-1)
    shape_objs["a2"] = ProcessTFArcs(c_x, c_y, x_1, y_1, x_2, y_2)

    # TF arc 3
    c_x = mf_data.data["xctfc(3)"].get_scan(-1)
    c_y = mf_data.data["yctfc(3)"].get_scan(-1)
    x_1 = mf_data.data["xarc(3)"].get_scan(-1)
    y_1 = mf_data.data["yarc(3)"].get_scan(-1)
    x_2 = mf_data.data["xarc(4)"].get_scan(-1)
    y_2 = mf_data.data["yarc(4)"].get_scan(-1)
    shape_objs["a3"] = ProcessTFArcs(c_x, c_y, x_1, y_1, x_2, y_2)

    # TF arc 4
    c_x = mf_data.data["xctfc(4)"].get_scan(-1)
    c_y = mf_data.data["yctfc(4)"].get_scan(-1)
    x_1 = mf_data.data["xarc(4)"].get_scan(-1)
    y_1 = mf_data.data["yarc(4)"].get_scan(-1)
    x_2 = mf_data.data["xarc(5)"].get_scan(-1)
    y_2 = mf_data.data["yarc(5)"].get_scan(-1)
    shape_objs["a4"] = ProcessTFArcs(c_x, c_y, x_1, y_1, x_2, y_2)

    # Extra information
    # TF coil inner leg thickness
    name = "TF Coil inner leg thickness (m)"
    value = t_tf
    shape_objs["i1"] = InfoHolder(name, value)

    # TF coil inner leg thickness
    name = "TF Coil outer leg thickness (m)"
    value = mf_data.data["tfthko"].get_scan(-1)
    shape_objs["i2"] = InfoHolder(name, value)

    # TF coil inner plasma facing case thickness
    name = "TF Coil inner leg plasma facing case thickness (m)"
    value = mf_data.data["casthi"].get_scan(-1)
    shape_objs["i3"] = InfoHolder(name, value)

    # TF coil inner case thickness
    name = "TF Coil inner case thickness (m)"
    value = mf_data.data["thkcas"].get_scan(-1)
    shape_objs["i4"] = InfoHolder(name, value)

    # TF coil winding pack thickness
    name = "TF Coil winding pack thickness (m)"
    value = mf_data.data["dr_tf_wp"].get_scan(-1)
    shape_objs["i5"] = InfoHolder(name, value)

    return shape_objs


def write_shapes_to_file(shape_data, filename):
    """Function to write the ellipse data to a MCNP file."""
    output_file = open(filename, "w")
    output_file.write("TZ x y z a b c\n")
    for shape in shape_data:
        if "e" in shape:
            dat = shape_data[shape]
            x = dat.x
            y = dat.y
            z = dat.z
            a = dat.a
            b = dat.b
            c = dat.c
            line = "TZ %.3f %.3f %.3f %.3f %.3f %.3f\n" % (x, y, z, a, b, c)
            output_file.write(line)
    output_file.write("\n")

    output_file.write("CZ x\n")
    for shape in shape_data:
        if "c" in shape:
            dat = shape_data[shape]
            dx = dat.dx
            line = "CZ %.3f\n" % dx
            output_file.write(line)
    output_file.write("\n")

    output_file.write("PZ y\n")
    for shape in shape_data:
        if "p" in shape:
            dat = shape_data[shape]
            dy = dat.dy
            line = "PZ %.3f\n" % dy
            output_file.write(line)

    output_file.write("\n")

    output_file.write("AZ c_x c_y x_1 y_1 x_2 y_2\n")
    for shape in shape_data:
        if "a" in shape:
            dat = shape_data[shape]
            c_x = dat.c_x
            c_y = dat.c_y
            x_1 = dat.x_1
            y_1 = dat.y_1
            x_2 = dat.x_2
            y_2 = dat.y_2
            line = "AZ %.3f %.3f %.3f %.3f %.3f %.3f\n" % (c_x, c_y, x_1, y_1, x_2, y_2)
            output_file.write(line)
    output_file.write("\n")

    output_file.write("Info y\n")
    for shape in shape_data:
        if "i" in shape:
            dat = shape_data[shape]
            name = dat.name
            value = dat.value
            line = "%s \t=\t %.3f\n" % (name, value)
            output_file.write(line)

    output_file.close()


if __name__ == "__main__":

    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Process MFILE.DAT into " "PROCESS.MCNP file."
    )

    parser.add_argument(
        "-f",
        metavar="f",
        type=str,
        default="MFILE.DAT",
        help="File to read as MFILE.DAT",
    )

    parser.add_argument(
        "-o",
        metavar="o",
        type=str,
        default="PROCESS.MCNP",
        help="File to write as PROCESS.MCNP",
    )

    parser.add_argument("--ctf", help="True/False flag for CTF", action="store_true")

    args = parser.parse_args()

    main(args)
