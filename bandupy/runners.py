from subprocess import Popen, PIPE, STDOUT
import os
import sys
# Imports from within the package
from .environ import BandUp_exe
from .files import (
    mkdir,
    create_bandup_input,
    create_bandup_plot_input,
)
from .plot import (
    print_opening_message,
    BandUpPlot, 
    produce_figure,
)

def run_bandup(args):
    start_dir = os.getcwd()
    # Running BandUP
    os.chdir(args.results_dir)
    bandup_run_options = [BandUp_exe] + args.argv
    with open("out_BandUP.dat", 'w') as f:
        bandup_run = Popen(bandup_run_options, stdout=PIPE, stderr=STDOUT)
        for line in iter(bandup_run.stdout.readline, ''):
            sys.stdout.write(line)
            f.write(line)
    os.chdir(start_dir)

def make_plot(args):
    print_opening_message()
    os.chdir(args.plotdir)
    plot = BandUpPlot(args)
    produce_figure(plot)

def run_requested_task(args):
    if(args.main_task=='unfold'):
        mkdir(args.results_dir, ignore_existing=True)
        create_bandup_input(args)
        run_bandup(args)
    elif(args.main_task=='plot'):
        mkdir(args.plotdir, ignore_existing=True)
        create_bandup_plot_input(args)
        make_plot(args)
    elif(args.main_task=='pre-unfold'):
       pass
    else:
        print 'Task "%s" not available.'%(args.main_task)
