from pathlib import Path
import os, shutil
import numpy as np

from process.io.in_dat import InDat

def main(main_name='cases', prefix = 'squid', create_scan=False):

    default_dir = 'stellarator_test/autostart'
    templates_dir = 'stellarator_test/templates'
    
    create_directory(Path(default_dir+'/'+main_name))

    B_min = 5
    B_max = 6
    B_list = list(np.linspace(B_min,B_max,(B_max-B_min)*10+1))
    print('B list: ', B_list)

    if create_scan:
        i = InDat(templates_dir+'/input.IN.DAT')

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

    i = InDat(templates_dir+'/input.IN.DAT')
    i.remove_iteration_variable(2)  # remove bt from iteration variables

    cases = B_list
    for case in cases:
        i.add_parameter("bt", case)
        case_path = default_dir+'/'+main_name+'/B_'+str(case)
        os.mkdir(case_path)
        i.write_in_dat()
        i.write_in_dat(output_filename=case_path+'/'+prefix+'.IN.DAT')
        shutil.copyfile(templates_dir+'/run_me.py', case_path+'/run_me.py')
        shutil.copyfile('stellarator_test/config_files/'+prefix+'.stella_conf.json', case_path+'/'+prefix+'.stella_conf.json',)


def create_directory(dirpath):
    if dirpath.exists() and dirpath.is_dir():
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

