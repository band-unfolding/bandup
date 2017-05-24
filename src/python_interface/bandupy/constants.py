# Copyright (C) 2017 Paulo V. C. Medeiros
# Constants used by the Python interface to BandUP
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
BANDUP_SRC_DIR = __get_parent_dir(PACKAGE_DIR, 2)
BANDUP_DIR = __get_parent_dir(PACKAGE_DIR, 3)
BANDUP_CONFIG_DIR = os.path.join(os.path.expanduser('~'), '.bandup')
BANDUP_CONFIG_FILE = os.path.join(BANDUP_CONFIG_DIR, 'config')

BANDUP_BIN = os.path.join(BANDUP_SRC_DIR, 'BandUP.x')
BANDUP_PRE_UNFOLDING_BIN = os.path.join(BANDUP_SRC_DIR, 'get_SCKPTS_pre_BandUP.x')

ORIGINAL_MATPLOTLIB_BACKEND = mpl.get_backend()

USER_HOME = os.path.expanduser('~')
WORKING_DIR = os.getcwd()

