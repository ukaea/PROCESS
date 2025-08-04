from process.io.mfile import MFile

from pathlib import Path
import os, shutil
import matplotlib.pyplot as plt


def main(main_name, prefix, param_x='bt', param_y = 'rmajor'):
    """
    Collect and plot output from MFILE.DAT in main_name directory
    prefix is a name of the MFILE.DAT file
    param is PROCESS parameter name loaded form the input
    """
    default_dir = 'stellarator_test/autostart'

    case_name = []
    results = []
    output = {}
    
    for case in os.listdir(default_dir+'/'+main_name):   
        mfile_path = os.path.join(default_dir+'/'+main_name, case, prefix+'.MFILE.DAT')
        if os.path.isfile(mfile_path):
            m = MFile(filename=mfile_path)

            if m.data[param_y].get_number_of_scans() == 1 and m.data['ifail'].get_scan(-1) == 1:
                # case_name.append(float(case[-3:]))
                case_name.append(m.data[param_x].get_scan(-1))
                results.append(m.data[param_y].get_scan(-1))
                output[m.data[param_x].get_scan(-1)] = m.data[param_y].get_scan(-1)

    output = dict(sorted(output.items()))
    print(output.keys())

    for key, value in output.items():
        print(key, value)
    plot_results(output.keys(), output.values(), param_y)


def plot_results(case_name, results, param):
    plt.plot(case_name, results, marker='o', linestyle='-')
    plt.xlabel('bt')
    plt.ylabel(param)
    plt.ylim((16.0, 25.0))
    plt.show()


if __name__ == "__main__":

    case_name = 'helias5_7T'
    prefix = 'helias5_7T'
    main(case_name, prefix = prefix)