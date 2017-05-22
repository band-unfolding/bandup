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
import sys

def arg_passed(arg_name):
    # Takes abbreviations into account
    possible_matches = []
    for passed_arg in sys.argv[1:]:
        if(passed_arg.strip()==arg_name[:len(passed_arg.strip())]):
            possible_matches.append(passed_arg)
    argname_passed = len(possible_matches) == 1
    return argname_passed
