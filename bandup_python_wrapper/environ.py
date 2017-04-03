# Copyright (C) 2017 Paulo V. C. Medeiros
# A python wrapper to BandUP and its plotting tool
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
import os
import sys
import glob

try:
    bandup_dir = os.environ['BANDUPDIR']
except(KeyError):
    msg = 'TEST'
    msg = ('The variable "BANDUPDIR" is not defined in your environment.\n' +
           18*' '+"Please set it yo BandUP's main directory.\n"+
           18*' '+'    > This is the directory where "build.sh" is located.\n'+
           18*' '+'Example:\n'+
           28*' '+'export BANDUPDIR=$HOME/codes/BandUP\n'+
           18*' '+'To have this automatically set next time you open a terminal, \n'+
           18*' '+'you can add the same command to, e.g., the file "~/.bash_profile".'
          )
    raise EnvironmentError(msg)

user_home = os.path.expanduser('~')
plot_path = os.path.join(bandup_dir, "utils", "post_unfolding", "plot")
plot_script = os.path.join(plot_path, "plot_unfolded_EBS_BandUP.py")
sys.path.insert(0, os.path.dirname(plot_script))

working_dir = os.getcwd()
BandUp_exe = os.path.join(bandup_dir, "BandUP_bin", "BandUP.x")

bandup_raw = {'symmetry_avgd':'unfolded_EBS_symmetry-averaged.dat',
              'not_symmetry_avgd':'unfolded_EBS_not-symmetry_averaged.dat'}

try:
    default_self_consist_calc_dir = glob.glob(
                                        os.path.join(os.path.dirname(working_dir),
                                                     'step_1*'))[0]
except(IndexError):
    default_self_consist_calc_dir = None

try:
    default_pre_unfolding_out_dir = glob.glob(
                                              os.path.join(os.path.dirname(working_dir),
                                                           'step_2*'))[0]
except(IndexError):
    default_pre_unfolding_out_dir = None
default_pre_unfolding_inputs_dir = None
if(default_pre_unfolding_out_dir is not None):
    default_pre_unfolding_inputs_dir = os.path.join(default_pre_unfolding_out_dir,
                                                    'input_files')


try:
    default_wavefunc_calc_dir = glob.glob(os.path.join(os.path.dirname(working_dir),
                                                       'step_3*'))[0]
except(IndexError):
    default_wavefunc_calc_dir = None

default_results_dir = os.path.join(working_dir, 'band_unfolding')
default_plot_dir = os.path.join(default_results_dir, 'plot')

default_pckpts_file = 'KPOINTS_prim_cell.in'
default_out_sckpts_file = 'KPOINTS_supercell.out'
default_pc_file = 'prim_cell_lattice.in'
default_sc_file = 'supercell_lattice.in'
