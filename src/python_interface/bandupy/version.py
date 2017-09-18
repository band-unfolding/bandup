# Copyright (C) 2017 Paulo V. C. Medeiros
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
from subprocess import Popen, PIPE
import os
from .constants import BANDUP_SRC_DIR

def get_tag_from_source():
    source_file = os.path.join(BANDUP_SRC_DIR, 'constants_and_types_mod.f90')
    source_tag = None
    with open(source_file, 'r') as f:
        for line in f:
            if('tag_for_push' in line and '::' in line):
                source_tag = (
                    line.split('::')[1].split('=')[1].replace('"','').replace("'","")
                )
                source_tag = 'v'+source_tag.strip()
                break
    return source_tag

def get_latest_git_tag():
    tag = None
    try:
        get_tag_run = Popen(['git', 'describe', '--tags', '--dirty'],
                               stdout=PIPE, stderr=PIPE)
        stdout, stderr = get_tag_run.communicate()
        if(not stderr.strip()):
            tag = stdout.strip()
    except(OSError):
        # In case git is not in the user's PATH
        pass
    return tag

def get_package_version():
    pv = get_latest_git_tag()
    if(pv is None):
        pv = get_tag_from_source()
    return pv

__version__ = get_package_version
