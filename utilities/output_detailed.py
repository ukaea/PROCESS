"""

  Create PROCESS output summary in html

  James Morris 06/06/2017
  CCFE

"""

# Third party libraries
import os
import sys
import json
import argparse
from collections import OrderedDict
from grip import export

# PROCESS libraries
from process.io.in_dat import InDat
from process.io.mfile import MFile
from process_io_lib.process_plots import (
    plot_pulse_timings,
    radial_bar_plot,
    plasma_profiles_plot,
)

try:
    import process_io_lib.process_dicts as proc_dict
except ImportError:
    print(
        "The Python dictionaries have not yet been created. Please run",
        " 'make dicts'!",
    )
    exit()

# Constants
EXCLUSIONS = ["ixc", "icc"]
MODULES = dict()
ICC_FULL = proc_dict.DICT_ICC_FULL
ICC_VARS = proc_dict.DICT_ICC_VARS
IXC_FULL = proc_dict.DICT_IXC_FULL
IXC_DEF = proc_dict.DICT_IXC_DEFAULT
DES_FIMP = proc_dict.DICT_FIMP
MODULES_FULL = proc_dict.DICT_MODULE
DEFAULTS = proc_dict.DICT_DEFAULT
DESCRIPTIONS = proc_dict.DICT_DESCRIPTIONS
FOM = proc_dict.DICT_OPTIMISATION_VARS


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Markdown


def bold(x):
    return "**" + str(x) + "**"


def unbold(x):
    return x.replace("*", "")


def code(x):
    return "`" + str(x) + "`"


def bold_italics(x):
    return "**_" + str(x) + "_**"


def output_line(x):
    OUTFILE.write(x + "\n\n")


def heading(x, y):
    OUTFILE.write("{0} {1}\n\n".format("#" * x, y))


def output_below_table_comment(x, y, z):
    OUTFILE.write("> **{0}** : **{1}** : {2}\n\n".format(x, y, z))


def output_newline():
    OUTFILE.write("\n\n")


def heading_comment(x):
    if x in COMMENTS:
        OUTFILE.write("> {0}\n\n".format(" ".join(COMMENTS[x])))


def return_to_top():
    OUTFILE.write("[:arrow_up:](#contents)\n\n")


def output_pic(x, y, z):
    OUTFILE.write("![{0}]({1} '{2}')\n\n".format(x, y, z))


def bold_info(x, y):
    OUTFILE.write("{0}: {1}\n\n".format(bold(x), y))


def content_heading(x):
    tag = x.lower().replace(" ", "-").replace(":", "")
    x = x.replace("-", " ")
    OUTFILE.write("[{0}](#{1})\n\n".format(x, tag))


def content_subheading(x):
    tag = x.lower().replace(" ", "-").replace(":", "")
    OUTFILE.write("{0}[{1}](#{2})\n\n".format("&nbsp;" * 8, x, tag))


def content_subsubheading(x):
    tag = x.lower().replace(" ", "-").replace(":", "")
    OUTFILE.write("{0}[{1}](#{2})\n\n".format("&nbsp;" * 16, x, tag))


def section_ro(section_name):
    """
    print mailto: for section RO
    """
    ro_name = DATA[section_name]["ro"]["name"]
    ro_email = DATA[section_name]["ro"]["email"]
    ro_lab = DATA[section_name]["ro"]["institute"]
    OUTFILE.write(
        "PMU RO: [{0}](mailto:{1}) ({2})\n\n".format(ro_name, ro_email, ro_lab)
    )


def table_heading(headings):
    header_line = "|"
    under_line = "|"
    for item in headings:
        header_line += " {0} |".format(item)
        under_line += " --- |"
    OUTFILE.write(header_line + "\n")
    OUTFILE.write(under_line + "\n")


def table_line(items, form):
    item_line = "|"
    item_format = "{0:" + "{0}".format(form) + "}"
    for item in items:
        if type(item) == float:
            item_line += item_format.format(item) + " |"
        elif type(item) == int:
            item_line += " {0} |".format(item)
        else:
            try:
                eval(item)
                item_line += item_format.format(item) + " |"
            except (NameError, ValueError, SyntaxError, TypeError):
                item_line += " {0} |".format(item)
    OUTFILE.write(item_line + "\n")


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def get_itvar_values():
    """
    Get iteration variable values from MFILE
    """
    itvar_vals = dict()

    for item in MFILE_LIST:
        if item[:5] == "itvar":
            name = MFILE.data[item].var_description
            value = MFILE.data[item].get_scan(-1)
            itvar_vals[name] = value
    return itvar_vals


def get_indat_comments():
    """
    Get IN.DAT #tag comments
    """
    lines = open(COMMAND_ARGS.f).readlines()

    for line in lines:
        if line[:2] == "*#":
            try:
                split_line = line.split(":", 1)
                key = split_line[0][2:].replace(" ", "")
                value = split_line[1].strip(" ").strip("\n")
                COMMENTS.setdefault(key, []).append(value)
            except Exception:
                print("Error with comment for line in IN.DAT :: {0}".format(line))


def constraint_comments(ky, iccs):
    """
    Display constraint comments for this module
    """

    for icc in iccs:
        con_key = "constraint-{0}".format(icc)
        if con_key in COMMENTS.keys():
            con_name = ICC_FULL[str(icc)]["name"]
            comment = " ".join(COMMENTS[con_key])
            if any(c.isalpha() for c in comment):
                output_below_table_comment(icc, con_name, comment)


def iteration_comments(ky, ixcs):
    """
    Display iteration comments for this module
    """

    for ixc in ixcs:
        itv_key = "iteration-variable-{0}".format(ixc)
        if itv_key in COMMENTS.keys():
            itv_name = IXC_FULL[str(ixc)]["name"]
            comment = " ".join(COMMENTS[itv_key])
            if len(comment) > 1:
                output_below_table_comment(ixc, itv_name, comment)


def output_contents():
    """
    Print contents page
    """

    bold_info("Username", MFILE.data["username"].get_scan(-1))
    bold_info("Date", MFILE.data["date"].get_scan(-1))
    bold_info("Time", MFILE.data["time"].get_scan(-1))
    bold_info("PROCESS version", MFILE.data["procver"].get_scan(-1))
    bold_info("PROCESS Tag Number", MFILE.data["tagno"].get_scan(-1))
    bold_info("Run description", MFILE.data["runtitle"].get_scan(-1))

    figure_of_merit = FOM[MFILE.data["minmax"].get_scan(-1)]
    if MFILE.data["minmax"].get_scan(-1) < 0:
        fom_type = "maximise"
    else:
        fom_type = "minimise"
    bold_info("Figure of merit", "{0} {1}".format(fom_type, figure_of_merit.lower()))

    bold_info(
        "PROCESS references",
        "[website](http://www.ccfe.ac.uk/powerplants.aspx), [physics paper](http://www.sciencedirect.com/science/article/pii/S0920379614005961), [engineering paper](http://www.sciencedirect.com/science/article/pii/S0920379616300072)",
    )
    bold_info(
        "Contacts",
        "[James Morris](mailto:james.morris2@ukaea.uk), [Hanni Lux](mailto:hanni.lux@ukaea.uk), [Michael Kovari](mailto:michael.kovari@ukaea.uk)",
    )

    # Top level comment
    output_line(bold("IN.DAT Comment"))
    heading_comment("header-title")

    heading(1, "Contents")

    content_heading("Output PDF")

    content_heading("Builds")
    content_subheading("Radial Build")
    content_subheading("Vertical Build")

    for k, v in DATA.items():
        if v["project"]:
            content_heading(k)

            if len(v["children"]) != 0:
                for child in v["children"]:
                    content_subheading(child)

                    if len(DATA[child]["children"]) != 0:
                        for child2 in DATA[child]["children"]:
                            content_subsubheading(child2)


def include_diagram(filename, caption):
    """
    Put requested image in output
    """
    global CAPTIONNUMBER

    try:
        output_pic(caption, filename, caption)
        output_line("Figure {0}: {1}".format(CAPTIONNUMBER, caption))
        CAPTIONNUMBER += 1
    except Exception:
        print("Cannot find image {0}".format(filename))
        pass


def output_constraints(k, numbers=False):
    """
    Output constraints list of constraint mods
    """

    MODULE_ICC_LIST = list()
    if numbers:
        cons = DATA[k]["constraints"]
        for c in cons:
            if c in CONSTRAINTS:
                MODULE_ICC_LIST.append(c)
            # else:
            #     print("Requested constraint {0} not a constraint for this MFILE".format(c))
    else:
        MODULE_ICC_LIST = sum(
            [MODULES[key] for key in MODULES.keys() if key in DATA[k]["constraints"]],
            [],
        )

    if len(MODULE_ICC_LIST) == 0:
        return

    heading(4, "Constraints")

    table_heading(
        [
            "Constraint",
            "Description",
            "Limit Name",
            "Limit",
            "Value",
            "F-Value Name",
            "F-Value Value",
        ]
    )

    for item in MODULE_ICC_LIST:

        con_name = ICC_FULL[str(item)]["name"]

        if ICC_VARS[str(item)] != "consistency" and ICC_VARS[str(item)] != "empty":
            if "f" in ICC_VARS[str(item)].keys():
                con_f_val_name = ICC_VARS[str(item)]["f"]
                if con_f_val_name not in MFILE_LIST:
                    if con_f_val_name in IT_VAR_VALUES.keys():
                        con_f_val = IT_VAR_VALUES[con_f_val_name]
                    else:
                        if con_f_val_name in INPUTS_LIST:
                            con_f_val = INFILE.data[con_f_val_name].value
                        else:
                            con_f_val = IXC_DEF[con_f_val_name]
                            con_f_val = bold(con_f_val)
                else:
                    con_f_val = MFILE.data[con_f_val_name].get_scan(-1)

                con_lim_name = ICC_VARS[str(item)]["v"]

                if con_lim_name in MFILE.data.keys():
                    con_lim = MFILE.data[con_lim_name].get_scan(-1)
                else:
                    if con_lim_name in DEFAULTS.keys():
                        con_lim = DEFAULTS[con_lim_name]
                    else:
                        con_lim = bold("calculated")
        else:
            con_f_val_name = "-"
            con_lim_name = "consistency"
            con_f_val = "-"
            con_lim = "-"

        if con_f_val != "-" and con_lim != "-" and con_lim != bold("calculated"):
            if type(con_lim) == str:
                con_lim = unbold(con_lim)
            if type(con_f_val) == str:
                con_f_val = unbold(con_f_val)
            actual_value = "{0:.2e}".format(float(con_lim) * float(con_f_val))
            con_lim = "{0:.2e}".format(float(con_lim))
            con_f_val = "{0:.2e}".format(float(con_f_val))
        else:
            actual_value = "-"

        table_line(
            [
                str(item),
                con_name,
                con_lim_name,
                con_lim,
                actual_value,
                con_f_val_name,
                con_f_val,
            ],
            ".4g",
        )

    constraint_comments(k, MODULE_ICC_LIST)


def output_itvars(k, numbers=False):
    """
    Output iteration variables for section k
    """

    if numbers:
        MODULES_IXC_NUMBERS_LIST = list()
        itvars_list = DATA[k]["iteration-variables"]
        for it in itvars_list:
            if it in IT_VARS:
                MODULES_IXC_NUMBERS_LIST.append(it)
            # else:
            #     print("Requested iteration variable {0} not iteration variable for this MFILE".format(it))
        MODULES_IXC_NAMES_LIST = [
            IXC_FULL[str(ixc)]["name"] for ixc in MODULES_IXC_NUMBERS_LIST
        ]

    else:
        mod_names = DATA[k]["iteration-variables"]
        MODULES_IXC_NAMES_LIST = [
            IXC_FULL[str(ixc)]["name"]
            for ixc in IT_VARS
            if any(
                IXC_FULL[str(ixc)]["name"] in MODULES_FULL[mod_name]
                for mod_name in mod_names
            )
        ]

        MODULES_IXC_NUMBERS_LIST = [
            ixc
            for ixc in IT_VARS
            if any(
                IXC_FULL[str(ixc)]["name"] in MODULES_FULL[mod_name]
                for mod_name in mod_names
            )
        ]

    if len(MODULES_IXC_NAMES_LIST) == 0:
        return

    heading(4, "Iteration Variables")
    output_line("* Values in **bold** are **not default** but user inputs.")
    table_heading(
        [
            "No.",
            "Name",
            "Final Value",
            "Description",
            "Starting Value",
            "Lower Bound",
            "Upper Bound",
        ]
    )

    for item in IT_VARS:

        item_name = IXC_FULL[str(item)]["name"]
        if "fimp(" in item_name:
            item_description = DES_FIMP[item_name]
        else:
            item_description = DESCRIPTIONS[item_name].split("\n")[0]

        if item_name in IT_VAR_VALUES.keys():
            item_value = IT_VAR_VALUES[item_name]
        else:
            item_value = MFILE.data[item_name].get_scan(-1)

        if item_name in INPUTS_LIST:
            starting_value = INFILE.data[item_name].value
            starting_value = bold(starting_value)
        else:
            if "fimp(" in item_name:
                starting_value = item_value
                fimp_index = int(item_name.split("(")[-1].split(")")[0])
                starting_value = INFILE.data["fimp"].value[fimp_index - 1]
            else:
                starting_value = IXC_DEF[item_name]

        if item_name in MODULES_IXC_NAMES_LIST:

            if str(item) in BOUNDS.keys():
                if "l" in BOUNDS[str(item)].keys():
                    low_bound = BOUNDS[str(item)]["l"]
                    low_bound = bold(low_bound)
                else:
                    low_bound = IXC_FULL[str(item)]["lb"]

                if "u" in BOUNDS[str(item)].keys():
                    up_bound = BOUNDS[str(item)]["u"]
                    up_bound = bold(up_bound)
                else:
                    up_bound = IXC_FULL[str(item)]["ub"]

            else:
                low_bound = IXC_FULL[str(item)]["lb"]
                up_bound = IXC_FULL[str(item)]["ub"]

            table_line(
                [
                    str(item),
                    code(item_name),
                    item_value,
                    item_description,
                    starting_value,
                    low_bound,
                    up_bound,
                ],
                ".4g",
            )

    iteration_comments(k, MODULES_IXC_NUMBERS_LIST)


def output_inputs(k, manual=False):
    """
    Output inputs for section k
    """

    mfile_values = list()
    if manual:
        MODULE_INPUT_LIST = list()
        inputs = DATA[k]["inputs"]
        for inp in inputs:
            if inp in INPUTS_LIST:
                MODULE_INPUT_LIST.append(inp)
            elif inp in MFILE_LIST:
                if (
                    MFILE.data[inp].var_flag != "OP"
                    and MFILE.data[inp].var_flag != "ITV"
                ):
                    mfile_values.append(inp)
                    MODULE_INPUT_LIST.append(inp)
            # else:
            #     print("The input you requested {0} is not in the input file".format(inp))
    else:
        mod_names = DATA[k]["input-sections"]
        MODULE_INPUT_LIST = [
            item
            for item in INPUTS_LIST
            if any(item in MODULES_FULL[mod_name] for mod_name in mod_names)
        ]

    if len(MODULE_INPUT_LIST) == 0:
        return

    heading(4, "Inputs")
    output_line(
        "* Values in **bold** are **default** input values. Output in MFILE but **NOT** specified in IN.DAT."
    )
    table_heading(["Input", "Value", "Description", "Comment"])

    for item in MODULE_INPUT_LIST:
        if (
            item not in IT_VAR_LIST
            and item not in EXCLUSIONS
            and item not in DATA[k]["exclusions"]
        ):

            if "fimp(" in item:
                item_name = "fimp"
            else:
                item_name = item

            if item_name not in INPUTS_LIST:
                item_value = str(MFILE.data[item].get_scan(-1))
            else:
                if "fimp(" not in item:
                    item_value = INFILE.data[item].value
                else:
                    fimp_index = int(item.split("(")[-1].split(")")[0])
                    item_value = str(INFILE.data["fimp"].value[fimp_index - 1])

            if "fimp(" not in item:
                item_des = DESCRIPTIONS[item].split("\n")[0]

            if "in-" + item in COMMENTS.keys():
                comment = " ".join(COMMENTS["in-" + item])
                if comment == "":
                    comment = ""
            else:
                comment = ""

            if type(item_value) == list:
                table_line([code(item), bold("array"), item_des], ".4g")
                for i in range(len(item_value)):
                    if "fimp" in item:
                        item_des = DES_FIMP["fimp({0})".format(i + 1)]
                    else:
                        item_des = "-"
                    table_line(
                        [
                            code(item) + "[{0}]".format(i),
                            item_value[i],
                            item_des,
                            comment,
                        ],
                        ".4g",
                    )
            elif "," in item_value:
                table_line([code(item), bold("array"), item_des], ".4g")
                item_value = item_value.split(",")
                for i in range(len(item_value)):
                    if item_value[i] != "":
                        if "fimp" in item:
                            item_des = DES_FIMP["fimp({0})".format(i + 1)]
                        else:
                            item_des = "-"
                        table_line(
                            [
                                code(item) + "[{0}]".format(i),
                                item_value[i],
                                item_des,
                                comment,
                            ],
                            ".4g",
                        )
            else:
                if item_name not in INPUTS_LIST:
                    table_line([code(item), bold(item_value), item_des, comment], ".4g")
                else:
                    table_line([code(item), item_value, item_des, comment], ".4g")


def output_outputs(k, mod_out_list=[]):
    """
    Output outputs
    """

    MODULE_OUTPUTS_LIST = DATA[k]["outputs"] + mod_out_list
    heading(4, "Outputs")
    table_heading(["Output", "Value", "Description"])

    for item in MODULE_OUTPUTS_LIST:
        item_value = MFILE.data[item].get_scan(-1)
        if item not in INPUTS_LIST and item not in DATA[k]["exclusions"]:

            try:
                item_des = MFILE.data[item].var_description.replace("_", " ")
                if item_des.lower() == item.lower():
                    item_name = ""
                else:
                    item_name = code(item)
                table_line([item_name.replace(" ", "_"), item_value, item_des], ".4g")
            except AttributeError:
                print("Warning: Skipping item:  {0}".format(item))
                pass


def output_output_sections(k):
    """
    Output the outputs from another sections
    """

    outputs_list = list()

    for item in DATA[k]["output-sections"]:
        outputs_list += DATA[item]["outputs"]

    output_outputs(k, mod_out_list=outputs_list)


def output_section(k):
    """
    Output various items from section depending on JSON config file
    """

    if len(DATA[k]["constraints"]) != 0:
        if DATA[k]["project"]:
            if type(DATA[k]["constraints"][0]) == str:
                output_constraints(k)
            else:
                output_constraints(k, numbers=True)
        else:
            output_constraints(k, numbers=True)

    if len(DATA[k]["iteration-variables"]) != 0:
        if DATA[k]["project"]:
            if type(DATA[k]["constraints"][0]) == str:
                output_itvars(k)
            else:
                output_itvars(k, numbers=True)
        else:
            output_itvars(k, numbers=True)

    if len(DATA[k]["input-sections"]) != 0:
        output_inputs(k)

    if len(DATA[k]["inputs"]) != 0:
        output_inputs(k, manual=True)

    if len(DATA[k]["output-sections"]) != 0:
        output_output_sections(k)

    if len(DATA[k]["outputs"]) != 0:
        output_outputs(k)


def output_section_topper(level, section_name, value):
    """
    outputs the top level information for the section
    """

    heading(level, section_name)

    if level == 1:
        section_ro(section_name)
    heading_name = "header-{0}".format(section_name.lower().replace(" ", "-"))
    heading_comment(heading_name)
    return_to_top()
    if value["diagrams"] != []:
        for item in value["diagrams"]:
            include_diagram(item["filename"], item["caption"])


def output_projects():
    """
    Output the projects
    """

    for k, v in DATA.items():
        if v["project"]:
            output_section_topper(1, k, v)
            output_section(k)

            if len(v["children"]) != 0:
                for child in v["children"]:
                    output_section_topper(2, child, DATA[child])
                    output_section(child)

                    if len(DATA[child]["children"]) != 0:
                        for child2 in DATA[child]["children"]:
                            output_section_topper(3, child2, DATA[child2])
                            output_section(child2)


def output_vertical_build():
    """
    Output for vertical build
    """

    heading(2, "Vertical Build")
    table_heading(["Name", "Thickness [m]", "Height [m]", "Description"])

    if MFILE.data["i_single_null"].get_scan(-1):
        vert_build = DATA["Builds"]["Vertical Build"]["single"]
    else:
        vert_build = DATA["Builds"]["Vertical Build"]["double"]

    tot_height = (
        MFILE.data["hmax"].get_scan(-1)
        - MFILE.data["shldlth"].get_scan(-1)
        - MFILE.data["divfix"].get_scan(-1)
        - MFILE.data["vgap"].get_scan(-1)
        + MFILE.data["vgaptop"].get_scan(-1)
        + MFILE.data["shldtth"].get_scan(-1)
        + MFILE.data["fwtth"].get_scan(-1)
        + MFILE.data["blnktth"].get_scan(-1)
        + MFILE.data["vvblgap"].get_scan(-1)
        + MFILE.data["tfcth"].get_scan(-1)
    )

    v_build_sum = tot_height

    for k, v in vert_build.items():

        item_name = vert_build[k]["name"]
        item_des = vert_build[k]["des"]
        item_combined = vert_build[k]["combination"]

        if item_combined:
            combo_type = vert_build[k]["combo_type"]
            items = item_name.replace(" ", "").split(combo_type)
            item_1_value = MFILE.data[items[0]].get_scan(-1)
            item_2_value = MFILE.data[items[1]].get_scan(-1)
            if combo_type == "+":
                item_value = item_1_value + item_2_value
            elif combo_type == "*":
                item_value = item_1_value * item_2_value
        else:
            item_value = MFILE.data[item_name].get_scan(-1)

        if int(k) < 10:
            table_line([code(item_name), item_value, v_build_sum, item_des], ".3f")
            v_build_sum -= item_value
        elif int(k) == 10:
            table_line([code(item_name), item_value, v_build_sum, item_des], ".3f")
            v_build_sum = 0
            table_line([code("Midplane"), 0.0, v_build_sum, "Device midplane"], ".3f")
        elif int(k) > 10:
            v_build_sum -= item_value
            table_line([code(item_name), item_value, v_build_sum, item_des], ".3f")


def output_radial_build():
    """
    Output radial build
    """

    heading(2, "Radial Build")

    try:
        radial_plot_file = DATA["Builds"]["radial_diagram"]["filename"]
        radial_plot_caption = DATA["Builds"]["radial_diagram"]["caption"]
        include_diagram(radial_plot_file, radial_plot_caption)
    except Exception:
        print("Cannot find image {0}".format(radial_plot_file["filename"]))
        pass

    table_heading(["Name", "Thickness [m]", "Radial Position [m]", "Description"])

    r_build_sum = 0
    rad_build = DATA["Builds"]["Radial Build"]

    for k, v in rad_build.items():

        item_name = rad_build[k]["name"]
        item_des = rad_build[k]["des"]
        item_combined = rad_build[k]["combination"]

        if item_combined:
            if rad_build[k]["combo_type"] == "+":
                item_1 = item_name.split("+")[0].replace(" ", "")
                item_2 = item_name.split("+")[1].replace(" ", "")
                item_1_value = MFILE.data[item_1].get_scan(-1)
                item_2_value = MFILE.data[item_2].get_scan(-1)
                item_value = item_1_value + item_2_value
            elif rad_build[k]["combo_type"] == "*":
                item_1 = item_name.split("*")[0].replace(" ", "")
                item_2 = item_name.split("*")[1].replace(" ", "")
                item_1_value = MFILE.data[item_1].get_scan(-1)
                item_2_value = MFILE.data[item_2].get_scan(-1)
                item_value = item_1_value * item_2_value
        else:
            item_value = MFILE.data[item_name].get_scan(-1)

        r_build_sum += item_value

        table_line([code(item_name), item_value, r_build_sum, item_des], ".3f")

    return


def output_modules():
    """
    Output the modules to markdown
    """

    icc = INFILE.data["icc"].value

    for item in icc:
        module = proc_dict.DICT_ICC_MODULE[str(item)]
        MODULES.setdefault(module, []).append(item)

    get_indat_comments()

    output_contents()

    heading(1, "Output PDF")
    include_diagram("process_diagram.png", "PROCESS output diagrams")

    output_section_topper(1, "Builds", DATA["Builds"])
    output_radial_build()
    output_vertical_build()

    output_projects()

    OUTFILE.close()


def create_plots():
    """
    Create all required plots for output
    """

    if not os.path.exists("documentation/figures/"):
        save_p = ""
    else:
        save_p = "documentation/figures/"

    # create pulse timings plot
    plot_pulse_timings(MFILE, save=True, show=False, save_path=save_p)

    # create radial build plot
    # TODO: non-hard coded list needed
    radial_build = [
        "bore",
        "ohcth",
        "precomp",
        "gapoh",
        "tfcth",
        "deltf",
        "thshield_ib",
        "gapds",
        "d_vv_in",
        "shldith",
        "vvblgap",
        "blnkith",
        "fwith",
        "scrapli",
        "rminor",
    ]
    radial_bar_plot(radial_build, MFILE, show=False, save=True, save_path=save_p)

    # create plasma profiles plot
    plasma_profiles_plot(MFILE, show=False, save=True, save_path=save_p)


def main(cargs):
    """
    Main
    """

    create_plots()

    output_modules()

    export(
        path=cargs.o + ".md",
        password="e35a21bfce5462bebbecc2e43d12bf4ec2ba469d",
        render_wide=True,
        render_inline=True,
        out_filename=cargs.o + ".html",
        title="PROCESS Output",
    )

    print("Over...")

    return


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        description="Create PROCESS output document."
        "For info contact james.morris2@ukaea.uk"
    )

    PARSER.add_argument(
        "-f", metavar="MFILENAME", type=str, default="", help="specify PROCESS MFILE"
    )

    PARSER.add_argument(
        "-j",
        metavar="JSONFILE",
        type=str,
        default="output_detailed.json",
        help="specify JSON file location and name",
    )

    PARSER.add_argument(
        "-o",
        metavar="OUTFILENAME",
        type=str,
        default="output_detailed",
        help="specify output file",
    )

    COMMAND_ARGS = PARSER.parse_args()

    # read json
    try:
        json_filename = COMMAND_ARGS.j
        DATA = json.load(open(json_filename), object_pairs_hook=OrderedDict)
    except ValueError as e:
        print("Error in JSON config file. Please correct error: {0}".format(e))
        sys.exit()
    except FileNotFoundError:
        print(
            "Cannot find JSON config file. Expecting file in: {0}".format(
                COMMAND_ARGS.j
            )
        )
        sys.exit()

    if COMMAND_ARGS.f:

        # read mfile
        MFILE = MFile(filename=COMMAND_ARGS.f)
        MFILE_LIST = MFILE.data.keys()

        # read input file
        INFILE = InDat(filename=COMMAND_ARGS.f, start_line=MFILE.mfile_end)
        INPUTS_LIST = INFILE.data.keys()

        # initialise COMMENTS dictionary
        COMMENTS = dict()

        CONSTRAINTS = INFILE.data["icc"].value
        IT_VARS = INFILE.data["ixc"].value
        IT_VAR_LIST = [IXC_FULL[str(itvar)]["name"] for itvar in IT_VARS]
        BOUNDS = INFILE.data["bounds"].value
        IT_VAR_VALUES = get_itvar_values()

        OUTFILE = open(COMMAND_ARGS.o + ".md", "w")

        # caption numbering
        CAPTIONNUMBER = 1

        main(COMMAND_ARGS)

    else:
        print("Please enter a reference MFILE with -f!")
