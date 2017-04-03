import sys
from subprocess import Popen, PIPE, STDOUT
import os
from bandup_python_wrapper.environ import BandUp_exe

def run_bandup(args):
    # Running BandUP
    os.chdir(args.results_dir)
    bandup_run_options = [BandUp_exe] + args.argv
    with open("out_BandUP.dat", 'w') as f:
        bandup_run = Popen(bandup_run_options, stdout=PIPE, stderr=STDOUT)
        for line in iter(bandup_run.stdout.readline, ''):
            sys.stdout.write(line)
            f.write(line)
