# Copyright (C) 2017 Paulo V. C. Medeiros
# Module used by the Python interface to BandUP to aid with list-related tasks
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
import numpy as np

def str_list_to_int_range(lst, convert_fortran2python_indexing=False, 
                          include_last_from_subranges=True):
    result_range = []
    for item in lst:
        try:
            int_item = int(item)
            if(convert_fortran2python_indexing): 
                int_item -= 1
            result_range.append(int_item)
        except(ValueError):
            start, end = map(int, item.split('-'))
            if(start>end): start, end = end, start
            if(convert_fortran2python_indexing):
                start -= 1
                end -= 1
            if(include_last_from_subranges):
                end += 1
            result_range += range(start, end)
    result_range = sorted(list(set(result_range)))
    return result_range

def int_list_to_str_range(lst):
    lst = sorted(list(set(lst)))
    diffs = np.diff(lst)
    intervals = []
    interval = [lst[0], lst[0]]
    for idiff, diff in enumerate(diffs):
        if(diff!=1):
            interval[1] = lst[idiff]
            intervals.append(tuple(interval))
            interval[0] = lst[idiff+1] 
        interval[1] = lst[idiff+1] 
    intervals.append(tuple(interval))

    to_join = []
    for interval in intervals:
        if(interval[0]==interval[1]):
            to_join.append(str(interval[0]))
        elif(interval[1]==interval[0]+1):
            to_join.append(str(interval[0]))
            to_join.append(str(interval[1]))
        else:
            to_join.append('%d-%d'%(interval[0], interval[1]))
    str_intervals = ', '.join(to_join)
    return str_intervals
