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

try:
    bandup_dir = os.environ['BANDUPDIR']
except(KeyError):
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

