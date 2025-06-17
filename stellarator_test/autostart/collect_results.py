from process.io.mfile import MFile

from pathlib import Path
import os, shutil
import matplotlib.pyplot as plt


def main(main_name='cases', prefix = 'squid'):

    default_dir = 'stellarator_test/autostart'

    case_name = []
    results = []
    param = 'rmajor'

    for case in os.listdir(default_dir+'/'+main_name):   
        mfile_path = os.path.join(default_dir+'/'+main_name, case, prefix+'.MFILE.DAT')
        m = MFile(filename=mfile_path)

        if m.data[param].get_number_of_scans() == 1:
            case_name.append(float(case[-3:]))
            results.append(m.data[param].get_scan(-1))

    
    print(case_name)
    print(results)
    plot_results(case_name, results)


def plot_results(case_name, results):
    plt.plot(case_name, results)
    plt.xlabel('bt')
    plt.ylabel('rmajor')
    plt.show()


if __name__ == "__main__":
    main('low_blanket', prefix = 'squid')