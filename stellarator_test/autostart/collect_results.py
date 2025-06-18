from process.io.mfile import MFile

from pathlib import Path
import os, shutil
import matplotlib.pyplot as plt


def main(main_name='cases', prefix = None, param = 'rmajor'):
    """
    Collect and plot output from MFILE.DAT in main_name directory
    prefix is a name of the MFILE.DAT file
    param is PROCESS parameter name loaded form the input
    """
    default_dir = 'stellarator_test/autostart'

    case_name = []
    results = []
    
    for case in os.listdir(default_dir+'/'+main_name):   
        mfile_path = os.path.join(default_dir+'/'+main_name, case, prefix+'.MFILE.DAT')
        m = MFile(filename=mfile_path)

        if m.data[param].get_number_of_scans() == 1:
            case_name.append(float(case[-3:]))
            results.append(m.data[param].get_scan(-1))

    
    print(case_name)
    print(results)
    plot_results(case_name, results, param)


def plot_results(case_name, results, param):
    plt.plot(case_name, results)
    plt.xlabel('bt')
    plt.ylabel(param)
    plt.show()


if __name__ == "__main__":
    main('updated_beta5', prefix = 'updated')