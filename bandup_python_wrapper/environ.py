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
import matplotlib as mpl
from .warnings_wrapper import warnings

package_dir = os.path.dirname(os.path.realpath(__file__))
package_parent_dir = os.path.dirname(package_dir)
ideal_bandup_src_dir = os.path.dirname(package_parent_dir)
ideal_bandup_dir = os.path.dirname(ideal_bandup_src_dir)
ideal_bandup_src_dir_is_valid = os.path.isfile(
                                    os.path.join(ideal_bandup_src_dir,'main_BandUP.f90')
                                )
ideal_bandup_dir_is_valid = ideal_bandup_src_dir_is_valid
print ideal_bandup_dir
try:
    bandup_dir = os.environ['BANDUPDIR']
    bandup_dir_accessible = os.access(bandup_dir, os.R_OK)
    if(not bandup_dir_accessible):
        raise EnvironmentError
except(EnvironmentError):
    msg = ('The variable "BANDUPDIR" is defined in your environment, but\n '+
           17*' '+'it points to an invalid path.\n'+
           18*' '+'Please either:\n'+
           18*' '+"(i)  set it to a valid path (BandUP's main directory), or \n"+
           18*' '+'(ii) delete it and try again (the code will try a workaround)')
    raise EnvironmentError(msg)
except(KeyError):
    if(ideal_bandup_dir_is_valid):
        msg = ('\n'+
               4*' '+'The variable "BANDUPDIR" is not defined in your environment.\n' +
               4*' '+"However, we have found what *seems to be* a BandUP dir at\n"+
               4*' '+"%s\n"%(ideal_bandup_dir) +
               4*' '+"This directory will be used. If you don't want this to happen,\n"+
               4*' '+"please set 'BANDUPDIR' to your preferred BandUP main dir.\n"+
               5*' '+'   > This should be a valid dir (normally where "build.sh" is)\n'+
               4*' '+'Example:\n'+
               14*' '+'export BANDUPDIR=$HOME/codes/BandUP\n'+
               4*' '+'To have this automatically set next time you open a terminal, \n'+
               4*' '+'add the same command to, e.g., the "~/.bash_profile" file.'
              )
        bandup_dir = ideal_bandup_dir
        warnings.warn(msg)
    else:
        msg = ('The variable "BANDUPDIR" is not defined in your environment.\n' +
               18*' '+"Please set it to BandUP's main directory.\n"+
               18*' '+'    > This is the directory where "build.sh" is located.\n'+
               18*' '+'Example:\n'+
               28*' '+'export BANDUPDIR=$HOME/codes/BandUP\n'+
               18*' '+'To have this automatically set next time you open a terminal, \n'+
               18*' '+'add the same command to, e.g., the "~/.bash_profile" file.'
              )
        raise EnvironmentError(msg)

original_matplotlib_backend = mpl.get_backend()
user_home = os.path.expanduser('~')
plot_path = os.path.join(bandup_dir, "utils", "post_unfolding", "plot")
plot_script = os.path.join(plot_path, "plot_unfolded_EBS_BandUP.py")
sys.path.insert(0, os.path.dirname(plot_script))


working_dir = os.getcwd()
BandUp_exe = os.path.join(bandup_dir, "BandUP_bin", "BandUP.x")

