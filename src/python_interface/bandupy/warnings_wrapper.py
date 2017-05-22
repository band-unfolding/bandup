# Copyright (C) 2017 Paulo V. C. Medeiros
# A convenient wrapper for the warnings method for use with BandUP
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
import warnings
import os
import sys

class WarningError(UserWarning):
    pass

def _new_showwarning(message, category = UserWarning, filename = '', lineno = -1):
    if(category==WarningError):
        warn_msg = 'ERROR '
    else:
        warn_msg = 'WARNING '
    warn_msg += '(%s, line %d): %s'%(os.path.basename(filename), lineno, 
                                           message)
    if(category==WarningError):
        warn_msg += '\nCannot continue. Stopping now.'
    print(warn_msg)
    if(category==WarningError):
        sys.exit(1)
warnings.showwarning = _new_showwarning

    
