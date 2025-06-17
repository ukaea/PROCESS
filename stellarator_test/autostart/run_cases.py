from pathlib import Path
import os, shutil
import subprocess

def main(main_name='cases', prefix = 'squid'):

    default_dir = 'stellarator_test/autostart'

    for case in os.listdir(default_dir+'/'+main_name):
        runpath = os.path.join(default_dir+'/'+main_name, case, 'run_me.py')
        subprocess.run(["python", runpath])

if __name__ == "__main__":
    main('low_blanket', prefix = 'squid')