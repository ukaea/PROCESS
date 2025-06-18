from pathlib import Path
import os, shutil
import subprocess

def main(main_name, prefix, skip_calculated=True):
    """
    Run cases in given main_name directory
    prefix is used to find input and config and then to name the output files
    cases which contain output files are skipped by default, to change that switch skip_calculated=False
    """
    default_dir = 'stellarator_test/autostart'

    for case in os.listdir(default_dir+'/'+main_name):
        if not (os.path.isfile(os.path.join(default_dir, main_name, case, prefix+'IN.DAT')) and skip_calculated):
            runpath = os.path.join(default_dir, main_name, case, 'run_me.py')
            subprocess.run(["python", runpath, '-n'+prefix])

if __name__ == "__main__":
    main('low_blanket', prefix = 'squid')