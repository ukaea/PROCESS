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
        # check if there is input file
        if os.path.isfile(os.path.join(default_dir, main_name, case, prefix+'.IN.DAT')):
            # check if recalculation should take place
            if not (os.path.isfile(os.path.join(default_dir, main_name, case, prefix+'.OUT.DAT'))) or not skip_calculated:
                print(f'\nCalculating case: {case}\n')
                runpath = os.path.join(default_dir, main_name, case, 'run_me.py')
                subprocess.run(["python", runpath, '-n'+prefix])
            else:
                print(f"Found results for {case}\n  Skipping...")
        else:
            print(f"Input at {os.path.join(default_dir, main_name, case)} not found")


if __name__ == "__main__":
    main('helias5_7T_2', prefix = 'helias5_7T')