'''
Generate directory structure and input files for selected case. 
06.2025 Walkowiak
'''

from pathlib import Path
import os, shutil
import numpy as np
import sys

from process.io.in_dat import InDat

def main( B, R, prefix = 'squid', main_name=None, clean_start=True):
    """
    Generate input files in the directory defined by main_name
    prefix is used to fine stella_conf file and input template 
    (if no input match the prefix, the default input.IN.DAT is used)
    TODO for now it works only for bt scan, a general version can be useful
    """
    default_dir = 'stellarator_test/MC_start'
    templates_dir = 'stellarator_test/templates'

    if os.path.isfile(templates_dir+'/'+prefix+'.IN.DAT'):
        input_file_path = templates_dir+'/'+prefix+'.IN.DAT'
    else:
        input_file_path = templates_dir+'/input.IN.DAT'
    
    create_directory(Path(default_dir+'/'+main_name), clean_start)


    i = InDat(input_file_path)
    i.remove_iteration_variable(2)  # remove bt from iteration variables
    i.remove_iteration_variable(3)  # remove Rmajor from iteration variables

    i.add_parameter("bt", B)
    i.add_parameter("rmajor", R)

    print(i.data['bounds'])

    ixc_list = i.data["ixc"].value
    for ixc in ixc_list:
        print(ixc, i.data[str(ixc)].value)

    sys.exit()

    iteriation_variables = []


    
    cases = B_list
    for case in cases:
        i.add_parameter("bt", case)
        case_path = default_dir+'/'+main_name+'/B_'+str(case)
        os.mkdir(case_path)
        i.write_in_dat()
        i.write_in_dat(output_filename=case_path+'/'+prefix+'.IN.DAT')
        shutil.copyfile(templates_dir+'/run_me.py', case_path+'/run_me.py')
        shutil.copyfile('stellarator_test/config_files/'+prefix+'.stella_conf.json', case_path+'/'+prefix+'.stella_conf.json',)


def create_directory(dirpath, clean_start=True):
    if dirpath.exists() and dirpath.is_dir() and clean_start:
        try:
            shutil.rmtree(dirpath)
        except Exception as e: print(e)

    try:
        os.mkdir(dirpath )
    except Exception as e: print(e)
    else: print(f'Fresh {dirpath} directory created') 
    

if __name__ == "__main__":
    # main('helias', prefix = 'helias')
    main(B=7.0, R=21.0, prefix='Helias_7T', main_name='test')

