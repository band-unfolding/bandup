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
from .warnings_wrapper import warnings


_package_dir = os.path.dirname(os.path.realpath(__file__))
_package_parent_dir = os.path.dirname(package_dir)
_bandup_src_dir = os.path.dirname(package_parent_dir)
BANDUP_DIR = os.path.dirname(ideal_bandup_src_dir)

_bandup_bin_dir = os.path.join(BANDUP_DIR, 'BandUP_bin')
BANDUP_BIN = os.path.join(_bandup_bin_dir, 'BandUP.x')

ORIGINAL_MATPLOTLIB_BACKEND = mpl.get_backend()

USER_HOME = os.path.expanduser('~')
WORKING_DIR = os.getcwd()

