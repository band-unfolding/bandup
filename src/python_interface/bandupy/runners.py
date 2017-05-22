# Copyright (C) 2013-2017 Paulo V. C. Medeiros
# This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
#
# BandUP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
#  BandUP is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with BandUP.  If not, see <http://www.gnu.org/licenses/>.
from subprocess import Popen, PIPE, STDOUT
import os
import sys
# Imports from within the package
from .constants import BANDUP_BIN, BANDUP_PRE_UNFOLDING_BIN, WORKING_DIR
from .files import (
    mkdir,
    create_bandup_input,
    create_bandup_plot_input,
)
from .plot import make_plot
from .vasp import procar2bandup
from .orbital_contributions import get_unfolded_orb_projs

def run_bandup(args):
    #start_dir = os.getcwd()
    start_dir = WORKING_DIR
    # Running BandUP
    os.chdir(args.results_dir)
    bandup_run_options = [BANDUP_BIN] + args.argv
    with open("out_BandUP.dat", 'w') as f:
        bandup_run = Popen(bandup_run_options, stdout=PIPE, stderr=STDOUT)
        for line in iter(bandup_run.stdout.readline, ''):
            sys.stdout.write(line)
            f.write(line)
    if(args.orbitals):
        get_orbital_projections_and_duals(args)
    os.chdir(start_dir)

def run_pre_bandup_tool(args):
    start_dir = WORKING_DIR
    # Running BandUP pre-unfolding tool
    os.chdir(args.inputs_dir)
    bandup_pre_unf_run_options = [BANDUP_PRE_UNFOLDING_BIN] + args.argv
    with open("out_BandUP_get_SCKPTS_pre_unfolding.dat", 'w') as f:
        bandup_pre_unf_run = Popen(bandup_pre_unf_run_options, 
                                   stdout=PIPE, stderr=STDOUT)
        for line in iter(bandup_pre_unf_run.stdout.readline, ''):
            sys.stdout.write(line)
            f.write(line)
    os.chdir(start_dir)


def run_requested_task(args):
    if(args.main_task=='unfold'):
        mkdir(args.results_dir, ignore_existing=True)
        create_bandup_input(args)
        run_bandup(args)
    elif(args.main_task=='plot'):
        if(args.gui):
            from .plot_gui.main_window import open_plot_gui
            open_plot_gui()
        else:
            mkdir(args.plotdir, ignore_existing=True)
            create_bandup_plot_input(args)
            make_plot(args)
    elif(args.main_task=='kpts-sc-get'):
       run_pre_bandup_tool(args)
    elif(args.main_task=='projected-unfold'):
        get_unfolded_orb_projs(args, clip_contributions=True, verbose=True)
    else:
        print('Task "%s" not available.'%(args.main_task))

def get_orbital_projections_and_duals(args):
    if(args.qe or args.castep or args.abinit):
        raise ValueError('Orbital projections not yet implemented for current PW code!')
    else:
        procar2bandup(fpath=os.path.join(args.wavefunc_calc_dir, 'PROCAR'))
