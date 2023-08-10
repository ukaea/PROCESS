"""

  Create PROCESS output document

  James Morris 11/05/2017
  CCFE

"""

# Third party libraries
import sys
import json
import argparse
from collections import OrderedDict
from grip import export

# PROCESS libraries
from process.io.in_dat import InDat
from process.io.mfile import MFile
from dict_tools import proc_dict

# Constants
MODULES = dict()
ICC_FULL = proc_dict.DICT_ICC_FULL
ICC_VARS = proc_dict.DICT_ICC_VARS
IXC_FULL = proc_dict.DICT_IXC_FULL
IXC_DEF = proc_dict.DICT_IXC_DEFAULT
DES_FIMP = proc_dict.DICT_FIMP
MODULES_FULL = proc_dict.DICT_MODULE
DEFAULTS = proc_dict.DICT_DEFAULT
DESCRIPTIONS = proc_dict.DICT_DESCRIPTIONS

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


def output_png(x, y, z):
    OUTFILE.write("![{0}]({1} '{2}')\n\n".format(x, y, z))


def bold_info(x, y):
    OUTFILE.write("{0}: {1}\n\n".format(bold(x), y))


def content_heading(x):
    tag = x.lower().replace(" ", "-").replace(":", "")
    x = x.replace("-", " ")
    OUTFILE.write("[{0}](#{1})\n\n".format(x, tag))


def content_subheading(x):
    tag = x.lower().replace(" ", "-").replace(":", "")
    OUTFILE.write("&nbsp;&nbsp;&nbsp;&nbsp;[{0}](#{1})\n\n".format(x, tag))


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
        if isinstance(item, float):
            item_line += item_format.format(item) + " |"
        elif isinstance(item, int):
            item_line += " {0} |".format(item)
        else:
            try:
                eval(item)
                item_line += item_format.format(item) + " |"
                # item_line += " {0:.4g} |".format(item)
            except (NameError, ValueError, SyntaxError, TypeError):
                item_line += " {0} |".format(item)
    OUTFILE.write(item_line + "\n")


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def check_empty(key):
    """
    Check not an empty section
    """
    con_mods = DATA["projects"][key]["icc_modules"]
    mod_names = DATA["projects"][key]["in_dat_modules"]
    tot_count = 0
    con_count = 0
    it_count = 0
    in_count = 0
    out_count = 0
    result = dict()

    CONSTRAINTS = INFILE.data["icc"].value
    IT_VARS = INFILE.data["ixc"].value
    IT_VAR_LIST = [IXC_FULL[str(itvar)]["name"] for itvar in IT_VARS]
    MODULE_VAR_LIST = sum([MODULES_FULL[mod_name] for mod_name in mod_names], [])
    MODULE_ICC_LIST = sum(
        [MODULES[key] for key in MODULES.keys() if key in con_mods], []
    )

    for item in CONSTRAINTS:
        if item in MODULE_ICC_LIST:
            tot_count += 1
            con_count += 1

    for item in IT_VARS:
        item_name = IXC_FULL[str(item)]["name"]
        if item_name in MODULE_VAR_LIST:
            tot_count += 1
            it_count += 1

    for item in INPUTS_LIST:
        if item in MODULE_VAR_LIST:
            if item not in IT_VAR_LIST:
                tot_count += 1
                in_count += 1

    for item in MFILE_LIST:
        if item in MODULE_VAR_LIST:
            if item not in INPUTS_LIST:
                tot_count += 1
                out_count += 1

    result["tot_count"] = tot_count
    result["con_count"] = con_count
    result["it_count"] = it_count
    result["in_count"] = in_count
    result["out_count"] = out_count

    return result


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
            split_line = line.split(":")
            key = split_line[0][2:].replace(" ", "")
            value = split_line[1].strip(" ")
            COMMENTS.setdefault(key, []).append(value)


def output_contents():
    """
    Print contents page
    """

    bold_info("Username", MFILE.data["username"].get_scan(-1))
    bold_info("Date", MFILE.data["date"].get_scan(-1))
    bold_info("Time", MFILE.data["time"].get_scan(-1))
    bold_info("PROCESS version", MFILE.data["procver"].get_scan(-1))
    bold_info("Run description", MFILE.data["runtitle"].get_scan(-1))

    # Top level comment
    output_line(bold("IN.DAT Comment"))
    heading_comment("header-title")

    heading(1, "Contents")

    content_heading("Diagrams")

    for k, v in DATA["projects"].items():

        result = check_empty(k)

        if result["tot_count"] != 0:

            content_heading(k)

            if k == "General":
                content_subheading("Radial Build")
                content_subheading("Vertical Build")

            # for item in v["mfile_sections"]:
            #     content_subheading(item)


def constraint_comments(ky, iccs):
    """
    Display constraint comments for this module
    """

    for icc in iccs:
        con_key = "constraint-{0}".format(icc)
        if con_key in COMMENTS.keys():
            con_name = ICC_FULL[str(icc)]["name"]
            comment = " ".join(COMMENTS[con_key])
            output_below_table_comment(icc, con_name, comment)


def output_constraints(k):
    """
    Output constraints for section k
    """

    con_mods = DATA["projects"][k]["icc_modules"]
    MODULE_ICC_LIST = sum(
        [MODULES[key] for key in MODULES.keys() if key in con_mods], []
    )

    heading(3, "Constraints")

    table_heading(
        [
            "Constraint",
            "Description",
            "F-Value Name",
            "F-Value Value",
            "Limit Name",
            "Limit",
            "Value",
        ]
    )

    for item in CONSTRAINTS:
        if item in MODULE_ICC_LIST:

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
                if isinstance(str):
                    con_lim = unbold(con_lim)
                if isinstance(con_f_val, str):
                    con_f_val = unbold(con_f_val)
                actual_value = "{0:.2e}".format(float(con_lim) * float(con_f_val))
                con_lim = "{0:.2e}".format(float(con_lim))
                con_f_val = "{0:.2e}".format(float(con_f_val))
            else:
                actual_value = "-"

            # if item in MODULES[con_mod]:
            table_line(
                [
                    str(item),
                    con_name,
                    con_f_val_name,
                    con_f_val,
                    con_lim_name,
                    con_lim,
                    actual_value,
                ],
                ".4g",
            )

    constraint_comments(k, MODULE_ICC_LIST)


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


def output_itvars(k):
    """
    Output iteration variables for section k
    """

    mod_names = DATA["projects"][k]["in_dat_modules"]
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

    heading(3, "Iteration Variables")
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


def output_inputs(k):
    """
    Output inputs for section k
    """

    mod_names = DATA["projects"][k]["in_dat_modules"]
    MODULE_INPUT_LIST = [
        item
        for item in INPUTS_LIST
        if any(item in MODULES_FULL[mod_name] for mod_name in mod_names)
    ]

    heading(3, "Inputs")
    table_heading(["Input", "Value", "Description", "Comment"])

    for item in INPUTS_LIST:
        if item in MODULE_INPUT_LIST:
            if item not in IT_VAR_LIST and item not in DATA["exclusions"]["inputs"]:

                item_value = INFILE.data[item].value
                item_des = DESCRIPTIONS[item].split("\n")[0]

                if "in-" + item in COMMENTS.keys():
                    comment = COMMENTS["in-" + item][0].replace("\n", "")
                    if comment == "":
                        comment = ""
                else:
                    comment = ""

                if "fimp(" in item:
                    i_des = DES_FIMP[item]

                if isinstance(item_value, list):
                    table_line([code(item), bold("array"), item_des], ".4g")
                    for i in range(len(item_value)):
                        if "fimp(" in item:
                            item_des = i_des[i]
                        else:
                            item_des = "-"
                        table_line(
                            [
                                code(item) + "[{0}]".format(i),
                                item_value[i],
                                item_des[4:],
                                comment,
                            ],
                            ".4g",
                        )
                elif "," in item_value:
                    table_line([code(item), bold("array"), item_des], ".4g")
                    item_value = item_value.split(",")
                    for i in range(len(item_value)):
                        if item_value[i] != "":
                            table_line(
                                [
                                    code(item) + "[{0}]".format(i),
                                    item_value[i],
                                    "-",
                                    comment,
                                ],
                                ".4g",
                            )
                else:
                    table_line([code(item), item_value, item_des, comment], ".4g")


def output_outputs(k):
    """
    Output outputs for section k
    """

    mod_names = DATA["projects"][k]["mfile_sections"]
    MODULE_OUTPUTS_LIST = [
        item
        for item in MFILE_LIST
        if any(MFILE.data[item].var_mod == mod_name for mod_name in mod_names)
    ]

    heading(3, "Outputs")
    table_heading(["Output", "Value", "Description"])

    for item in MFILE_LIST:
        if item in MODULE_OUTPUTS_LIST:

            item_value = MFILE.data[item].get_scan(-1)
            if item not in INPUTS_LIST:

                try:
                    item_des = MFILE.data[item].var_description.replace("_", " ")
                    if item_des.lower() == item.lower():
                        item_name = ""
                    else:
                        item_name = code(item)
                    table_line(
                        [item_name.replace(" ", "_"), item_value, item_des], ".4g"
                    )
                except AttributeError:
                    print("Warning: Skipping item:  {0}".format(item))
                    pass


def vertical_build():
    """
    Output for vertical build
    """

    heading(2, "Vertical Build")
    table_heading(["Name", "Thickness [m]", "Height [m]", "Description"])

    if MFILE.data["i_single_null"].get_scan(-1):
        vert_build = DATA["machine build"]["vertical"]["single"]
    else:
        vert_build = DATA["machine build"]["vertical"]["double"]

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


def radial_build():
    """
    Output the machine radial build
    """

    heading(2, "Radial Build")
    table_heading(["Name", "Thickness [m]", "Radial Position [m]", "Description"])

    r_build_sum = 0
    rad_build = DATA["machine build"]["radial"]

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


def output_project(key, result):
    """
    Output sections
    """

    heading(1, key)
    return_to_top()
    heading_comment("header-{0}".format(key))

    if key == "General":
        radial_build()
        vertical_build()

    # Constraints
    if result["con_count"] != 0:
        output_constraints(key)

    if result["it_count"] != 0:
        output_itvars(key)

    if result["in_count"] != 0:
        output_inputs(key)

    if result["out_count"] != 0:
        output_outputs(key)


def output_diagrams():
    """
    Embed plot_proc output
    """

    heading(1, "Diagrams")
    try:
        output_png("alt text", "process_diagram.png", "PROCESS output diagrams")
    except Exception:
        pass


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

    output_diagrams()

    for k, v in DATA["projects"].items():

        result = check_empty(k)

        if result["tot_count"] != 0:

            output_project(k, result)

    OUTFILE.close()


def main(cargs):
    """
    Main
    """

    output_modules()

    export(
        path="output_summary.md",
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
        "-o",
        metavar="OUTFILENAME",
        type=str,
        default="output_summary",
        help="specify output file",
    )

    COMMAND_ARGS = PARSER.parse_args()

    # read json
    try:
        DATA = json.load(open("output_summary.json"), object_pairs_hook=OrderedDict)
    except ValueError as e:
        print("Error in JSON config file. Please correct error: {0}".format(e))
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

        main(COMMAND_ARGS)

    else:
        print("Please enter a reference MFILE with -f!")
