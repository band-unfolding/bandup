# Copyright (C) 2017 Paulo V. C. Medeiros
# Routines used by the build script
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
from __future__ import print_function
from subprocess import Popen, PIPE, STDOUT
import shutil
import os
from .files import mkdir, rmdir
from .constants import BANDUP_SRC_DIR, BANDUP_CONFIG_FILE
from .warnings_wrapper import warnings, WarningError

def system_has_compiler(compiler):
    try:
        compiler_call = Popen([compiler], stdout=PIPE, stderr=STDOUT)
        return True
    except(OSError):
        return False
def get_gcc():
    comp = 'gcc'
    for v in range(2,7)[::-1]:
        new_comp = 'gcc-%d'%(v)
        if(system_has_compiler(new_comp)):
            comp = new_comp
            break
    return comp
def compatible_c(fortcomp):
    if('ifort' in fortcomp): return 'icc'
    elif('gfortran' in fortcomp or 
         'nagfor' in fortcomp): 
        return get_gcc()
    else: return None
def compatible_cpp(fortcomp):
    if('ifort' in fortcomp or
       'nagfor' in fortcomp): 
        return 'fpp'
    elif('gfortran' in fortcomp):
        return 'cpp'
    else: return None
def compatible_omp_flags(fortcomp):
    fcflag, cflag = "", ""
    if('ifort' in fortcomp):
        fcflag, cflag = '-openmp', '-openmp'
    elif('gfortran' in fortcomp):
        fcflag, cflag = '-fopenmp', '-fopenmp'
    elif('nagfor' in fortcomp):
        fcflag, cflag = '-openmp', '-fopenmp'
    return fcflag, cflag
def assert_valid_compiler(args, supported_fortran_compilers):
    user_chose_compiler = args.compiler is not None
    if(user_chose_compiler):
        if(not system_has_compiler(args.compiler)):
            msg = 'The compiler you chose (%s) is not available.'%(args.compiler)
            warnings.warn(msg, category=WarningError)
        print('BandUP will be compiled using %s (compiler chosen by you)'%(
              args.compiler))
    else:
        for comp in supported_fortran_compilers:
            if(system_has_compiler(comp)):
                args.compiler = comp
                break
        if(args.compiler is None):
            msg = 'Could not find a supported Fortran compiler.'
            warnings.warn(msg, category=WarningError)
        print('BandUP will be compiled using %s (compiler determined automatically)'%(
              args.compiler))
    return args

def castep_interface_available(calling_from_build_script=False):
    ret = False
    if(calling_from_build_script):
        check2xsf_dir = os.path.join(BANDUP_SRC_DIR, 'external', 
                                     'check2xsf2_modules_for_BandUP')
        castep_modules_dir = os.path.join(BANDUP_SRC_DIR, 'castep_related_modules')
        ret = os.path.isdir(check2xsf_dir) and os.path.isdir(castep_modules_dir)
    else:
        with open(BANDUP_CONFIG_FILE, 'r') as f:
            for line in f:
                if('CASTEP_SUPPORT' in line):
                    support = line.split('=')[1].strip()
                    if(support.lower() == 'true'):
                        ret = True
    return ret

def qe_interface_available(calling_from_build_script=False):
    ret = False
    if(calling_from_build_script):
        qe_needed_lib_file = os.path.join(BANDUP_SRC_DIR, 'external', 
                                          'espresso-5.1_modules_for_BandUP',
                                          'flib', 'flib.a')
        ret = os.path.isfile(qe_needed_lib_file)
    else:
        with open(BANDUP_CONFIG_FILE, 'r') as f:
            for line in f:
                if('QE_SUPPORT' in line):
                    support = line.split('=')[1].strip()
                    if(support.lower() == 'true'):
                        ret = True
    return ret
