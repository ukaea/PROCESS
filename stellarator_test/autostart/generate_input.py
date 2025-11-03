'''
Generate directory structure and input files for selected case. 
06.2025 Walkowiak
'''

from pathlib import Path
import os, shutil
import numpy as np
import warnings

from process.io.in_dat import InDat

def main(main_name, prefix, B_min = 5.0, B_max = 6.0, create_scan=False, clean_start=True):
    """
    Generate input files in the directory defined by main_name
    prefix is used to fine stella_conf file and input template 
    (if no input match the prefix, the default input.IN.DAT is used)
    TODO for now it works only for bt scan, a general version can be useful
    """
    default_dir = 'stellarator_test/autostart'
    templates_dir = 'stellarator_test/templates'

    input_file_path = os.path.join(templates_dir+'/'+prefix+'.IN.DAT')
    if not os.path.isfile(input_file_path):
        warnings.warn('\nNo input file found with given prefix, using default\n', stacklevel=2)
        input_file_path = templates_dir+'/input.IN.DAT'
        
    
    create_directory(Path(default_dir+'/'+main_name), clean_start)

    B_list = list(np.linspace(B_min,B_max,round((B_max-B_min)*10+1)))
    print('B list: ', B_list)

    if create_scan:
        i = InDat(input_file_path)

        i.remove_iteration_variable(2)  # remove bt from iteration variables
        i.add_parameter('nsweep', 28)   # variable selection: 28 -> bt
        i.add_parameter('isweep', len(B_list))    # number of scan points
        i.add_parameter("sweep", B_list)    # scan points

        case = 'scan'
        case_path = default_dir+'/'+main_name+'/'+str(case)
        os.mkdir(case_path)
        i.write_in_dat()
        i.write_in_dat(output_filename=case_path+'/'+prefix+'.IN.DAT')
        shutil.copyfile(templates_dir+'/run_me.py', case_path+'/run_me.py')
        shutil.copyfile('stellarator_test/config_files/'+prefix+'.stella_conf.json', case_path+'/'+prefix+'.stella_conf.json',)

    i = InDat(input_file_path)
    i.remove_iteration_variable(2)  # remove bt from iteration variables

    
    cases = B_list
    for case in cases:
        case_path = default_dir+'/'+main_name+'/B_'+f'{case:.1f}'
        if not os.path.isdir(case_path):
            os.mkdir(case_path)

        i.add_parameter("bt", case)
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
    main('low_blanket', prefix = 'squid')

