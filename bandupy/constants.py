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
import matplotlib as mpl

def __get_parent_dir(path, level=1):
    parent = path
    for ilevel in range(level):
        parent = os.path.dirname(parent)
    return parent


PACKAGE_DIR = os.path.dirname(os.path.realpath(__file__))
INTERFACE_MAIN_SOURCE_DIR = __get_parent_dir(PACKAGE_DIR, 1)
BANDUP_DIR = __get_parent_dir(PACKAGE_DIR, 3)

BANDUP_BIN_DIR = os.path.join(BANDUP_DIR, 'BandUP_bin')
BANDUP_BIN = os.path.join(BANDUP_BIN_DIR, 'BandUP.x')
BANDUP_PRE_UNFOLDING_BIN = os.path.join(BANDUP_DIR, 'utils', 'pre_unfolding', 
                                        'get_SCKPTS_pre_BandUP', 
                                        'get_SCKPTS_pre_BandUP.x')

PACKAGE_VERSION = 'UNKNOWN'
with open(os.path.join(BANDUP_DIR, 'src', 'general_io_mod.f90'), 'r') as f:
    for line in f:
        if('package_version' in line and '::' in line):
            PACKAGE_VERSION = line.split('=')[-1].replace('"','').strip()
            break

ORIGINAL_MATPLOTLIB_BACKEND = mpl.get_backend()

USER_HOME = os.path.expanduser('~')
WORKING_DIR = os.getcwd()

