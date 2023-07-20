"""Useful functions for the GUI

"""

from dicts.gui_dicts import GUI_MODULE, GUI_LABLXC, GUI_LABLCC

from lib.guiindat import GuiInDat


def dict_to_in_dat(di):
    """Converts a dictionary recieved from the client to a
    GuiInDat object. The ixc and icc arrays are created
    from the itervar_#, constraint_# style checkboxes

    """
    in_dat = GuiInDat()
    for dummy, tuplist in GUI_MODULE.items():
        # all of the usual variables that appear in GUI_MODULE
        # can just be assigned
        for tu in tuplist:
            varname = tu[1]
            in_dat[varname] = di[varname]

    # do iteration variables
    ixc = []
    for num in GUI_LABLXC.keys():
        itername = "itervar_" + str(num)
        boundlname = "boundl(" + str(num) + ")"
        bounduname = "boundu(" + str(num) + ")"
        # add to ixc if variable is set as itervar
        if itername in di and di[itername] == "on":
            ixc.append(num)
        # copy the bounds
        in_dat[boundlname] = di[boundlname]
        in_dat[bounduname] = di[bounduname]

    in_dat["ixc"] = ixc
    in_dat["nvar"] = len(ixc)

    # do constraint equations
    icc = []
    # f = open('test_text.txt', 'w')
    # for k in di:
    #     f.write('\n' + k)
    # f.write('\n\nNext is GUI_LABCC keys')
    for num in GUI_LABLCC.keys():
        # f.write('\n' + str(num))
        constname = "constraint_" + str(num)
        if constname in di and di[constname] == "on":
            icc.append(num)
    #     f.write('\n' + constname)
    # f.close()

    in_dat["icc"] = icc
    in_dat["neqns"] = len(icc)

    # do description
    if "Run_Description" in di:
        in_dat.Run_Description = di["Run_Description"]

    return in_dat


def split_dicts(request):
    """Converts a submission from the client to two dictionaries.
    Reference values start with 'ref_'. These are split into
    their own dictionary so they can be dealt with seperately

    """
    assert request.method == "POST"
    in_dict = {}
    ref_dict = {}
    for key, value in request.POST.items():
        if "ref_" in key:
            ref_dict[key[4:]] = value
        else:
            in_dict[key] = value

    return in_dict, ref_dict
