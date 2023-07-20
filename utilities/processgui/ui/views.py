"""Views for the GUI.

"""
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext

from process_io_lib.process_dicts import DICTIONARY_VERSION

from dicts.gui_dicts import GUI_MODULE, GUI_LABLXC, GUI_LABLCC
from lib.mylib import split_dicts, dict_to_in_dat
from lib.guiindat import GuiInDat


def index(request):
    """The view for the main page"""
    if request.method == "GET":
        # show the page with default values
        context = create_context()
        req_cont = RequestContext(request)
        return render_to_response("ui.html", context, context_instance=req_cont)

    assert request.method == "POST"
    # split the submitted data into two dictionaries. One for reference values
    # one for input values
    in_dict, ref_dict = split_dicts(request)

    # f = open('test_views.txt', 'w')
    # f.write('\n\nNext is in_dict')
    # for k in in_dict:
    #     f.write('\n' + k)
    # f.write('\n\nNext is ref_dict')
    # for k in ref_dict:
    #     f.write('\n' + k)

    if "submit_download" in request.POST:
        # if the user wants to download the file
        in_dat_string = str(dict_to_in_dat(in_dict))
        response = HttpResponse(in_dat_string, content_type="text/DAT")
        response["Content-Disposition"] = 'attachment; filename="IN.DAT"'
        return response

    elif "inputfile" in request.FILES:
        # if the user has uploaded an input file
        text = request.FILES["inputfile"].read().decode()
        in_dat = GuiInDat()
        in_dat.readlines(text.split("\n"))
        # read in the new in_dat values from the file and the
        # reference values from the form
        ref_dat = dict_to_in_dat(ref_dict)

    elif "reffile" in request.FILES:
        # if the user has uploaded a reference file
        text = request.FILES["reffile"].read().decode()
        ref_dat = GuiInDat()
        ref_dat.readlines(text.split("\n"))
        # read in the ref_dat values from the file and the input
        # values from the form
        in_dat = dict_to_in_dat(in_dict)

    context = create_context(in_dat, ref_dat)
    req_cont = RequestContext(request)
    return render_to_response("ui.html", context, context_instance=req_cont)


def echo(request):
    """Echos every submitted value back"""
    assert request.method == "POST"
    ret = ""
    for key, value in sorted(request.POST.items()):
        ret += key + " : " + value + "<br>"
    return HttpResponse(ret)


def create_context(in_dat=None, ref_dat=None):
    """Given two GuiInDat objects, creates a context used to expand
    the template

    """
    context = {}
    if not in_dat:
        in_dat = GuiInDat()
    if not ref_dat:
        ref_dat = GuiInDat()

    context["dictionaryversion"] = DICTIONARY_VERSION
    context["moduledata"] = GUI_MODULE
    context["inputvalues"] = in_dat
    context["refvalues"] = ref_dat
    context["lablxc"] = GUI_LABLXC
    context["lablcc"] = GUI_LABLCC

    return context
